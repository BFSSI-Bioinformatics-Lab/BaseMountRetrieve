import os
import click
import shutil
import logging
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from dataclasses import dataclass
from typing import Optional
from tabulate import tabulate
from BaseMountRetrieve.__init__ import __version__, __author__, __email__

logger = logging.getLogger()
logging.basicConfig(
    format=f'\033[92m \033[1m [%(levelname)s]\t\033[0m %(message)s',
    level=logging.INFO)


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    logging.info(f"Version: {__version__}")
    logging.info(f"Author: {__author__}")
    logging.info(f"Email: {__email__}")
    quit()


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    return Path(value)


@dataclass
class BaseMountSample:
    """
    Dataclass to store data for a Sample from BaseMount (sample_id, R1, R2)
    """
    run_dir: Path
    sample_dir: Path  # ~/basemount/Projects/PRJ1/Samples/SAMPLE1
    sample_id: str
    sample_name: str

    project_dir: Path = None
    r1: Path = None
    r2: Path = None

    def __post_init__(self):
        self.fastq_dir = self.sample_dir / 'Files'
        try:
            self.validate_fastq_in_sample_dir(self.fastq_dir)
        except AssertionError:
            logging.error(
                f"ERROR: Could not validate FASTQ files for {self.sample_id} ({self.sample_name}) at {self.sample_dir}")
            return

        r1, r2 = self.get_reads()
        self.r1 = r1
        self.r2 = r2

    def get_reads(self) -> tuple:
        """
        Grabs fwd and rev reads from the Files/ directory
        """
        r1 = list(self.fastq_dir.glob("*_R1_*"))[0]
        r2 = list(self.fastq_dir.glob("*_R2_*"))[0]
        return r1, r2

    @staticmethod
    def validate_fastq_in_sample_dir(fastq_dir: Path):
        """
        Validates fwd and rev reads are present in the Files/ directory
        """
        r1 = list(fastq_dir.glob("*_R1_*"))
        r2 = list(fastq_dir.glob("*_R2_*"))
        assert len(r1) == 1
        assert len(r2) == 1
        assert r1[0].is_file()
        assert r2[0].is_file()


@dataclass
class BaseMountRun:
    """
    Dataclass to store data for a Run from BaseMount
    """
    run_dir: Path

    project_dir: Path = None
    runinfoxml: Path = None
    stats_json: Path = None
    runparametersxml: Path = None
    interop_dir: Path = None
    interop_files: [Path] = None
    experiment_name: str = None

    def __post_init__(self):
        self.properties_dir = self.run_dir / 'Properties'

        self.run_id = self.run_dir.name
        self.log_dir = self.get_log_dir()

        # Get samplesheet, read into df
        self.samplesheet = self.get_samplesheet()
        if self.experiment_name is None:
            self.experiment_name = self.extract_experiment_name(samplesheet=self.samplesheet)
            logger.debug(f'Set experiment name to {self.experiment_name}')
        self.samplesheet_df = self.parse_samplesheet(samplesheet=self.samplesheet)

        # Store this to verify that that FASTQ samples match up with this - requeued runs will break this matchup and
        # we need to log a warning
        samplesheet_samples = self.samplesheet_df['Sample_ID'].tolist()

        # Print out the samplesheet in the console; useful for debugging
        print(tabulate(self.samplesheet_df, headers='keys', tablefmt='psql'))

        # Get InterOp dir and file contents
        self.interop_dir = self.get_interop_dir()

        if self.interop_dir is not None:
            self.interop_files = self.get_interop_files()

        # Get samples and set up BasemountSample objects
        self.sample_dirs = self.get_sample_dirs()
        self.sample_id_dict = self.get_sample_id_dict()
        self.sample_objects = self.generate_sample_objects()

        # Do cross validation between sample_objects and samplesheet_samples
        sample_object_ids = [s.sample_id for s in self.sample_objects]
        difference_list = list(set(sample_object_ids) - set(samplesheet_samples))
        if len(difference_list) > 0:
            logger.warning(f'Found discrepanices between listed samples in the samplesheet and samples found in FASTQ '
                           f'directories. This is likely due to an analysis requeue in BaseSpace. Discrepant samples:')
            for s in difference_list:
                print(f'\t\t{s}')

        # Get runinfo, runparameters, stats
        self.runinfoxml = self.get_runinfoxml()
        self.runparametersxml = self.get_runparametersxml()

        # Get all log files
        if self.log_dir is not None:
            self.stats_json = self.get_stats_json()
            self.logfiles = self.get_log_files()

    def get_log_dir(self) -> Optional[Path]:
        log_dir = self.run_dir / 'Logs'
        if log_dir.exists():
            logger.debug(f'Detected log directory at {log_dir}')
            return log_dir
        else:
            # NOTE: The BaseMount API is a mess, so the entire Output.Samples folder must be searched to return the
            # proper log directory. If the directory can't be found, will return None.
            samples_dir = self.properties_dir / 'Output.Samples'
            sample_folders = list(samples_dir.glob("*"))
            logger.warning(f'Could not detect standard log directory at {log_dir}, digging a bit deeper...')

            # Iterate over the sample folders til we find one with the proper Logs directory, then return
            for sample_folder in sample_folders:
                app_session_log_dir = sample_folder / 'ParentAppSession' / 'Logs'
                if app_session_log_dir.is_dir():
                    logger.debug(f'Detected alternate log directory at {app_session_log_dir}')
                    return app_session_log_dir
        logger.error('ERROR: Could not find any log directory after deep search!')
        return None

    def get_interop_files(self) -> [Path]:
        """
        Collects all files in the InterOp/ directory and filters out any junk (directories, hidden files)
        """
        interop_files = list(self.interop_dir.glob("*"))
        interop_files = [f for f in interop_files if f.is_file()]
        interop_files = [f for f in interop_files if not f.name.startswith(".")]
        return interop_files

    def get_log_files(self) -> [Path]:
        """
        Collects all files in the Logs/ directory and filters out hidden files
        :return: List containing paths to all log files
        """
        if self.log_dir is None:
            return []
        logfiles = list(self.log_dir.glob("*"))
        logfiles = [f for f in logfiles if not f.name.startswith(".")]
        return logfiles

    def get_stats_json(self) -> Optional[Path]:
        """
        Method to retrieve the Stats.json file from self.log_dir
        :return: Path to Stats.json or None if nothing can be found
        """
        stats_json = self.log_dir / 'Stats.json'
        if stats_json.exists():
            return stats_json
        return None

    def get_interop_dir(self) -> Path:
        """
        Collects and verifies the InterOp directory for a particular Run
        """
        interop_dir = self.run_dir / 'Files' / 'InterOp'
        if not interop_dir.is_dir():
            logging.warning(
                f"Could not locate InterOp data for {self.experiment_name} at {interop_dir}. "
                f"Confirm researcher shared Run on BaseSpace.")
            interop_dir = None
        return interop_dir

    def get_sample_id_dict(self) -> dict:
        """
        Extracts the Sample IDs from a BaseMount run directory, creates a dict containing Sample_ID:Sample_dir links
        """
        sample_id_dict = {}
        for sample_dir in self.sample_dirs:
            sample_properties = sample_dir / 'SampleProperties'
            sample_id = None
            sample_name = None
            if sample_properties.exists():
                with open(str(sample_properties), 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if 'Name:' in line:
                            sample_name = line.split(" ", 1)[1].strip()
                        elif 'SampleId:' in line:
                            sample_id = line.split(" ", 1)[1].strip()
                        else:
                            continue
            sample_id_dict[sample_id] = {'sample_dir': sample_dir, 'sample_name': sample_name}

        # Debugging
        logger.debug(f"Extracted the following sample details for {self.run_id}:")
        logger.debug(f"sample_id\t\tsample_name")
        for key, val in sample_id_dict.items():
            logger.debug(f"{key}\t{val['sample_name']}")

        return sample_id_dict

    def get_samplesheet(self) -> Path:
        """
        Grabs and verifies the SampleSheet for a Run. Checks several known locations on BaseMount.
        """
        possible_samplesheet_locations = [
            self.properties_dir / "Input.sample-sheet",
            self.properties_dir.parent / 'Files' / 'SampleSheet.csv',
            self.properties_dir / 'Input.Libraries' / '0' / 'Properties' / 'Output.Runs' / '0' / 'Files' / 'SampleSheet.csv',
        ]

        for samplesheet in possible_samplesheet_locations:
            if samplesheet.is_file():
                logger.debug(f'Found SampleSheet at {samplesheet}')
                return samplesheet
            else:
                continue
        raise FileNotFoundError(f"Could not find SampleSheet in any of the expected locations!")

    def generate_sample_objects(self) -> [BaseMountSample]:
        """
        Generates a list of BaseMountSample data objects belonging to the respective Run
        """
        logging.debug(f"Generating BaseMountSample objects for {self.run_id} ({self.experiment_name})")
        sample_object_list = []
        for sample_id, metadata in self.sample_id_dict.items():
            sample_object = BaseMountSample(project_dir=self.project_dir,
                                            run_dir=self.run_dir,
                                            sample_dir=metadata['sample_dir'],
                                            sample_id=sample_id,
                                            sample_name=metadata['sample_name'])
            sample_object_list.append(sample_object)
        return sample_object_list

    def get_sample_dirs(self) -> list:
        """
        Grabs all Samples directories for the Run, filters out junk.
        Searches in ../Run/Properties/Output.Samples
        """
        sample_dirs = list(self.run_dir.glob("Sample.*"))
        # Try another location
        if len(sample_dirs) == 0:
            sample_dirs = []
            sample_dirs_unfiltered = list((self.run_dir / 'Properties' / 'Output.Samples').glob("*"))
            # Sometimes these can be empty, so filter them out
            for sample_dir in sample_dirs_unfiltered:
                if (sample_dir / 'SampleProperties').exists():
                    sample_dirs.append(sample_dir)

        # Filter out junk
        sample_dirs = [sample_dir for sample_dir in sample_dirs if 'Undetermined' not in sample_dir.name]
        sample_dirs = [sample_dir for sample_dir in sample_dirs if sample_dir.is_dir()]

        return sample_dirs

    def get_runparametersxml(self) -> Optional[Path]:
        """
        # TODO: Check if this file exists anywhere else on BaseMount
        Tries to grab the RunParameters.xml file if its present in the expected location
        """
        runparametersxml_1 = self.properties_dir / 'Input.Runs' / '0' / 'Files' / 'RunParameters.xml'
        runparametersxml_2 = self.run_dir / 'Files' / 'RunParameters.xml'
        if runparametersxml_1.is_file():
            return runparametersxml_1
        elif runparametersxml_2.is_file():
            return runparametersxml_2
        else:
            logging.warning(f"Could not locate RunParameters.xml for {self.experiment_name}")
            return None

    def get_runinfoxml(self) -> Optional[Path]:
        """
        Grabs the RunInfo.xml file
        """

        # Try to find RunInfo.xml
        runinfoxml_1 = self.properties_dir / 'Input.Runs' / '0' / 'Files' / 'RunInfo.xml'
        runinfoxml_2 = self.run_dir / 'Logs' / 'RunInfo.xml'
        runinfoxml_3 = self.run_dir / 'Files' / 'RunInfo.xml'
        if runinfoxml_1.is_file():
            return runinfoxml_1
        elif runinfoxml_2.is_file():
            return runinfoxml_2
        elif runinfoxml_3.is_file():
            return runinfoxml_3
        else:
            logging.warning(f"Could not locate RunInfo.xml for {self.experiment_name}")
            return None

    @staticmethod
    def parse_samplesheet(samplesheet: Path) -> pd.DataFrame:
        """
        Reads SampleSheet.csv and returns dataframe (all header information will be stripped)
        :param samplesheet: Path to SampleSheet.csv
        :return: pandas df of SampleSheet.csv with head section stripped away
        """
        counter = 1
        with open(str(samplesheet)) as f:
            for line in f:
                if '[Data]' in line:
                    break
                else:
                    counter += 1
        df = pd.read_csv(samplesheet, sep=",", index_col=False, skiprows=counter)
        return df

    @staticmethod
    def extract_experiment_name(samplesheet: Path) -> str:
        """
        Retrieves the 'Experiment Name' from SampleSheet.csv
        :param samplesheet: Path to SampleSheet.csv
        :return: value of 'Experiment Name'
        """
        with open(str(samplesheet)) as f:
            for line in f:
                if 'Experiment Name' in line:
                    experiment_name = line.split(',')[1].strip()
                    return experiment_name
                elif 'Description' in line:
                    experiment_name = line.split(',')[1].strip()
                    return experiment_name
            else:
                raise Exception(f"Could not find 'Experiment Name' in {samplesheet}")


@dataclass
class BaseMountProject:
    """
    Dataclass to store data for a Project from BaseMount
    """
    project_dir: Path

    def __post_init__(self):
        self.sample_dir = self.project_dir / 'Samples'
        self.appsessions = self.project_dir / 'AppSessions.v1'
        self.run_dirs = self.get_runs_dirs()
        self.run_objects = self.generate_run_objects()

    def generate_run_objects(self) -> [BaseMountRun]:
        """
        Generates a list of BaseMountRun data objects belonging to this Project
        """
        run_object_list = []
        for run_dir in self.run_dirs:
            run_object = BaseMountRun(run_dir=run_dir, project_dir=self.project_dir)
            run_object_list.append(run_object)
        return run_object_list

    def get_runs_dirs(self) -> list:
        """
        Grabs all valid Run directories
        """
        run_dirs = list(self.appsessions.glob('FASTQ*'))
        return run_dirs


def retrieve_experiment_contents_from_basemount(run_dir: Path, out_dir: Path, rename: bool):
    logging.info(f"Started retrieving contents of {run_dir.name}")

    # Create output directory if it doesn't already exist
    out_dir.mkdir(exist_ok=True, parents=True)

    # Validate directory
    if not run_dir.exists():
        raise FileNotFoundError(f"Directory {run_dir} does not exist! Quitting.")
    experiment_name = run_dir.name

    # Establish Run object
    run_obj = BaseMountRun(run_dir=run_dir, experiment_name=experiment_name)

    # Copy reads and run data to out_dir
    create_run_folder_skeleton(out_dir=out_dir)
    copy_metadata_files(run_obj=run_obj, out_dir=out_dir)
    copy_log_files(run_obj=run_obj, out_dir=out_dir)
    copy_interop_files(run_obj=run_obj, out_dir=out_dir)
    copy_reads(run_obj=run_obj, out_dir=out_dir, rename=rename)


def retrieve_project_contents_from_basemount(project_dir: Path, out_dir: Path, rename: bool):
    """
    Main method to analyze BaseMount folder contents, establish dataclasses (Project, Run, Sample), and copy to out_dir
    """
    logging.info(f"Started retrieving contents {project_dir.name}")
    # Create output directory if it doesn't already exist
    out_dir.mkdir(exist_ok=True)

    logging.info(f"Analyzing contents of {project_dir} ...")
    project = BaseMountProject(project_dir=project_dir)
    for run_obj in project.run_objects:
        logging.info(f"Processing {run_obj.run_id}...")

        # Setup output directory
        run_dir_out = out_dir / run_obj.run_id

        # Check if outdir already exists, skip if it does
        if run_dir_out.exists():
            logging.info(f"Run directory for {run_obj.run_id} already exists, skipping")
            continue

        create_run_folder_skeleton(out_dir=run_dir_out)
        copy_metadata_files(run_obj=run_obj, out_dir=run_dir_out)
        copy_log_files(run_obj=run_obj, out_dir=run_dir_out)
        copy_interop_files(run_obj=run_obj, out_dir=run_dir_out)
        copy_reads(run_obj=run_obj, out_dir=run_dir_out, rename=rename)


def shutil_if_exists(src: Path, dst: Path):
    """
    Copies a file to dst, sets generous permissions, and fails quietly if src is None
    """
    try:
        shutil.copy(str(src), str(dst))
        os.chmod(str(dst), 0o666)  # 666: read/write but no execute
    except FileNotFoundError:
        return


def copy_metadata_files(run_obj: BaseMountRun, out_dir: Path):
    logging.debug(f"Copying metadata for {run_obj.run_id}")
    metadata_files = [
        (run_obj.samplesheet, out_dir / 'SampleSheet.csv'),
        (run_obj.runinfoxml, out_dir / 'RunInfo.xml'),
        (run_obj.stats_json, out_dir / 'Stats.json')
    ]
    for src, dst in metadata_files:
        shutil_if_exists(src, dst)


def copy_log_files(run_obj: BaseMountRun, out_dir: Path):
    logging.debug(f"Copying log file contents for {run_obj.run_id}")
    for logfile in run_obj.logfiles:
        shutil.copy(str(logfile), str(out_dir / 'Logs' / logfile.name))
        os.chmod(str(out_dir / 'Logs' / logfile.name), 0o666)


def copy_interop_files(run_obj: BaseMountRun, out_dir: Path):
    logging.debug(f"Searching for InterOp files...")
    if run_obj.interop_files is not None:
        logging.debug(f"Copying InterOp contents for {run_obj.run_id}")
        for interop_file in run_obj.interop_files:
            interop_file_out = out_dir / 'InterOp' / interop_file.name
            shutil.copy(str(interop_file), str(interop_file_out))
            os.chmod(str(interop_file_out), 0o666)
    else:
        logging.debug(f"InterOp files not available for {run_obj.run_id}, skipping")


def copy_reads(run_obj: BaseMountRun, out_dir: Path, rename: bool):
    logging.info(f"Copying reads for {run_obj.run_id}...")
    for sample_obj in tqdm(iterable=run_obj.sample_objects, miniters=1):
        # Rename samples to {SampleID_(R1/R2).fastq.gz} if flag is set, otherwise keep as-is
        if rename:
            r1_out = out_dir / 'Data' / 'Intensities' / 'BaseCalls' / (sample_obj.sample_id + "_R1.fastq.gz")
            r2_out = out_dir / 'Data' / 'Intensities' / 'BaseCalls' / (sample_obj.sample_id + "_R2.fastq.gz")
        else:
            r1_out = out_dir / 'Data' / 'Intensities' / 'BaseCalls' / sample_obj.r1.name
            r2_out = out_dir / 'Data' / 'Intensities' / 'BaseCalls' / sample_obj.r2.name
        shutil.copy(str(sample_obj.r1), str(r1_out))
        shutil.copy(str(sample_obj.r2), str(r2_out))
        os.chmod(str(r1_out), 0o666)
        os.chmod(str(r2_out), 0o666)


def create_run_folder_skeleton(out_dir: Path):
    """
    Creates the skeleton structure for a mock local MiSeq run
    :param out_dir: This should be the path to .../project_name/experiment_name
    """
    base_folders = ['Config',
                    (Path('Data') / Path('Intensities') / Path('BaseCalls')),
                    'Images',
                    'InterOp',
                    'Logs',
                    'Recipes',
                    'Thumbnail_Images']
    out_dir.mkdir(exist_ok=True)
    for f in base_folders:
        (out_dir / f).mkdir(parents=True, exist_ok=True)


@click.command(help="BaseMountRetrieve will tap into the mounted BaseMount filesystem and retrieve an entire Project "
                    "(or single Run) in the output style of a local MiSeq run.")
@click.option('-p', '--project-dir',
              type=click.Path(exists=True),
              required=False,
              default=None,
              help='Path to the directory on BaseMount for a particular Project. e.g. ../basemount/Projects/[project]. '
                   'Cannot be used at the same time as the --run-dir parameter.',
              callback=convert_to_path)
@click.option('-r', '--run-dir',
              type=click.Path(exists=True),
              required=False,
              default=None,
              help='Path to the directory on BaseMount for a particular Run. e.g. ../basemount/Runs/[your run]. '
                   'Cannot be used at the same time as the --project-dir parameter.',
              callback=convert_to_path)
@click.option('-o', '--out-dir',
              type=click.Path(exists=False),
              required=True,
              default=None,
              help='Directory to dump all runs for project.',
              callback=convert_to_path)
@click.option('--rename',
              help='Use this flag to automatically re-name the R1 and R2 files to just include the Sample ID.',
              is_flag=True,
              default=False)
@click.option('-v', '--verbose',
              help='Use this flag to enable more verbose output.',
              is_flag=True,
              default=False)
@click.option('--version',
              help='Use this flag to print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              expose_value=False)
def cli(project_dir, run_dir, out_dir, rename, verbose):
    logging.info(f"Started BaseMountRetrieve (v{__version__})")

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug("Enabled VERBOSE mode")

    # Input validation
    if project_dir is not None and run_dir is not None:
        logging.error(f"ERROR: Please only provide one of --project-dir or --run-dir, not both!")
        quit()
    if project_dir is None and run_dir is None:
        logging.error(f"ERROR: Please provide one of --project-dir or --run-dir!")
        quit()
    if project_dir and project_dir.parent.name != 'Projects':
        logging.error(f"ERROR: Parent folder name is not as expected (expected: 'Projects', "
                      f"actual: '{project_dir.parent.name}'). "
                      f"Ensure you have provided the proper path to your target Project directory from the mounted "
                      f"BaseMount file system.")
        quit()
    if run_dir and run_dir.parent.name != 'Runs':
        logging.error(f"ERROR: Parent folder name is not as expected (expected: 'Runs', "
                      f"actual: '{run_dir.parent.name}'). "
                      f"Ensure you have provided the proper path to your target Run directory from the mounted "
                      f"BaseMount file system.")
        quit()

    logging.debug(f"Rename samples: {rename}")
    logging.debug(f"out_dir:\t{out_dir}")

    if project_dir:
        logging.debug(f"project_dir:\t{project_dir}")
        retrieve_project_contents_from_basemount(project_dir=project_dir, out_dir=out_dir, rename=rename)
    elif run_dir:
        logging.debug(f"run_dir:\t{run_dir}")
        retrieve_experiment_contents_from_basemount(run_dir=run_dir, out_dir=out_dir, rename=rename)
    logging.info("Done!")


if __name__ == "__main__":
    cli()

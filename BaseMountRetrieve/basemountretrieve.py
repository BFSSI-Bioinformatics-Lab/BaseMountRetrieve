import os
import click
import shutil
import logging
import pandas as pd
from tqdm import tqdm
from typing import Union
from pathlib import Path
from dataclasses import dataclass
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
    runparametersxml: Path = None
    interop_dir: Path = None
    interop_files: [Path] = None
    experiment_name: str = None

    def __post_init__(self):
        self.log_dir = self.run_dir / 'Logs'
        self.properties_dir = self.run_dir / 'Properties'
        self.run_id = self.run_dir.name

        # Get samplesheet, read into df
        self.samplesheet = self.get_samplesheet()
        if self.experiment_name is None:
            self.experiment_name = self.extract_experiment_name(samplesheet=self.samplesheet)
        self.samplesheet_df = self.parse_samplesheet(samplesheet=self.samplesheet)

        # Get InterOp dir and file contents
        self.interop_dir = self.get_interop_dir()

        if self.interop_dir is not None:
            self.interop_files = self.get_interop_files()

        # Get samples and set up BasemountSample objects
        self.sample_dirs = self.get_sample_dirs()
        self.sample_id_dict = self.get_sample_id_dict()
        self.sample_objects = self.generate_sample_objects()

        # Get runinfo and runparameters
        self.runinfoxml = self.get_runinfoxml()
        self.runparametersxml = self.get_runparametersxml()

        # Get all logfiles
        self.logfiles = self.get_log_files()

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
        """
        logfiles = list(self.log_dir.glob("*"))
        logfiles = [f for f in logfiles if not f.name.startswith(".")]
        return logfiles

    def get_interop_dir(self) -> Path:
        """
        Collects and verifies the InterOp directory for a particular Run
        """
        interop_dir = self.run_dir / 'Files' / 'InterOp'
        if not interop_dir.is_dir():
            logging.warning(
                f"Could not locate InterOp data for {self.experiment_name}. "
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
        return sample_id_dict

    def get_samplesheet(self) -> Path:
        """
        Grabs and verifies the SampleSheet for a Run
        """
        samplesheet = self.properties_dir / "Input.sample-sheet"
        if samplesheet.is_file():
            return samplesheet
        else:
            samplesheet = self.properties_dir.parent / 'Files' / 'SampleSheet.csv'
            try:
                assert (samplesheet.exists())
            except FileNotFoundError:
                raise FileNotFoundError(f"Could not find SampleSheet at expected location: {samplesheet}")
            return samplesheet

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
        Grabs all Samples directories for the Run, filters out junk
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

    def get_runparametersxml(self) -> Union[Path, None]:
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

    def get_runinfoxml(self) -> Union[Path, None]:
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
            logging.error(f"Could not locate RunInfo.xml for {self.experiment_name}")
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

    def validate_samplesheet(self):
        pass


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


def retrieve_experiment_contents_from_basemount(experiment_name: str, basemount_dir: Path, out_dir: Path, rename: bool):
    # Create output directory if it doesn't already exist
    out_dir.mkdir(exist_ok=True, parents=True)

    # Validate /Runs directory
    runs_dir = basemount_dir / 'Runs'
    if not runs_dir.exists():
        logging.error(f"ERROR: Directory {runs_dir} does not exist!")
        quit()

    experiment_names = [x.name for x in list(runs_dir.glob("*"))]
    experiment_dir = None
    if experiment_name in experiment_names:
        experiment_dir = runs_dir / experiment_name
        assert experiment_dir.exists()
    else:
        logging.error(f"ERROR: Could not find experiment {experiment_name} in {runs_dir}")
        logging.debug(f"Detected the following experiment names: {experiment_names}")
        quit()

    # Establish Run object
    run_obj = BaseMountRun(run_dir=experiment_dir, experiment_name=experiment_name)

    create_run_folder_skeleton(out_dir=out_dir)
    copy_metadata_files(run_obj=run_obj, out_dir=out_dir)
    copy_log_files(run_obj=run_obj, out_dir=out_dir)
    copy_interop_files(run_obj=run_obj, out_dir=out_dir)
    copy_reads(run_obj=run_obj, out_dir=out_dir, rename=rename)


def retrieve_project_contents_from_basemount(project_dir: Path, out_dir: Path, rename: bool):
    """
    Main method to analyze BaseMount folder contents, establish dataclasses (Project, Run, Sample), and copy to out_dir
    """

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


def copy_metadata_files(run_obj: BaseMountRun, out_dir: Path):
    logging.debug(f"Copying metadata for {run_obj.run_id}")
    shutil.copy(str(run_obj.samplesheet), str(out_dir / 'SampleSheet.csv'))
    shutil.copy(str(run_obj.runinfoxml), str(out_dir / 'RunInfo.xml'))
    os.chmod(str(out_dir / 'RunInfo.xml'), 0o775)
    os.chmod(str(out_dir / 'SampleSheet.csv'), 0o775)
    if run_obj.runparametersxml is not None:
        shutil.copy(str(run_obj.runparametersxml), str(out_dir / 'RunParameters.xml'))
        os.chmod(str(out_dir / 'RunParameters.xml'), 0o775)


def copy_log_files(run_obj: BaseMountRun, out_dir: Path):
    logging.debug(f"Copying log file contents for {run_obj.run_id}")
    for logfile in run_obj.logfiles:
        shutil.copy(str(logfile), str(out_dir / 'Logs' / logfile.name))
        os.chmod(str(out_dir / 'Logs' / logfile.name), 0o775)


def copy_interop_files(run_obj: BaseMountRun, out_dir: Path):
    logging.debug(f"Searching for InterOp files...")
    if run_obj.interop_files is not None:
        logging.debug(f"Copying InterOp contents for {run_obj.run_id}")
        for interop_file in run_obj.interop_files:
            interop_file_out = out_dir / 'InterOp' / interop_file.name
            shutil.copy(str(interop_file), str(interop_file_out))
            os.chmod(str(interop_file_out), 0o775)
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
        os.chmod(str(r1_out), 0o775)
        os.chmod(str(r2_out), 0o775)


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


@click.command(help="BaseMountRetrieve will tap into the mounted BaseMount filesystem and retrieve all of the runs for"
                    " a given project in the output style of a local MiSeq run.")
@click.option('-p', '--project-dir',
              type=click.Path(exists=True),
              required=False,
              default=None,
              help='Path to the directory on BaseMount for a particular project. e.g. '
                   'basemount/Projects/[your project].',
              callback=convert_to_path)
@click.option('-e', '--experiment-name',
              type=click.STRING,
              required=False,
              default=None,
              help='Searches BaseMount for the name of a run/experiment and attempts to retrieve its contents. '
                   'Must be used along with the --basemount-dir parameter. '
                   'Cannot be used alongside the --project-dir flag.')
@click.option('-b', '--basemount-dir',
              type=click.Path(exists=True),
              required=False,
              default=None,
              help='Path to root directory for BaseMount, e.g. ~/basemount/. '
                   'Must be supplied alongside --experiment-name parameter.',
              callback=convert_to_path)
@click.option('-o', '--out-dir',
              type=click.Path(exists=False),
              required=True,
              default=None,
              help='Directory to dump all runs for project.',
              callback=convert_to_path)
@click.option('-r', '--rename',
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
def cli(project_dir, experiment_name, basemount_dir, out_dir, rename, verbose):
    logging.info(f"Started BaseMountRetrieve (v{__version__})")

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug("Enabled VERBOSE mode")

    # Input validation
    if project_dir is not None and experiment_name is not None:
        logging.error(f"ERROR: Please only provide one of --project-dir or --experiment-name, not both!")
        quit()
    if project_dir is None and experiment_name is None:
        logging.error(f"ERROR: Must provide a value for --project-dir or --experiment-name")
        quit()
    if not (experiment_name is not None and basemount_dir is not None):
        logging.error("ERROR: Must provide both --experiment-name and --basemount-dir values")
        quit()

    logging.debug(f"out_dir:\t\t{out_dir}")
    logging.debug(f"Rename samples: {rename}")

    # TODO: Conduct check to see if runs already exist at an earlier stage to save waiting time
    if project_dir:
        logging.debug(f"project_dir:\t\t{project_dir}")
        retrieve_project_contents_from_basemount(project_dir=project_dir, out_dir=out_dir, rename=rename)
    elif experiment_name:
        logging.debug(f"experiment_name:\t\t{experiment_name}")
        logging.debug(f"basemount_dir:\t\t{basemount_dir}")
        retrieve_experiment_contents_from_basemount(experiment_name=experiment_name, basemount_dir=basemount_dir,
                                                    out_dir=out_dir, rename=rename)
    logging.info("Done!")


if __name__ == "__main__":
    cli()

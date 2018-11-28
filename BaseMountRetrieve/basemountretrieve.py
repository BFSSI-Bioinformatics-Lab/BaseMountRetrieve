import os
import click
import shutil
import logging
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from dataclasses import dataclass
from BaseMountRetrieve.__init__ import __version__, __author__, __email__

logger = logging.getLogger()
logging.basicConfig(
    format=f'\033[92m \033[1m [%(levelname)s]\tbasemountretrieve:\033[0m %(message)s',
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
    project_dir: Path
    run_dir: Path
    sample_dir: Path  # ~/basemount/Projects/PRJ1/Samples/SAMPLE1
    sample_id: str

    r1: Path = None
    r2: Path = None

    def __post_init__(self):
        self.fastq_dir = self.sample_dir / 'Files'
        try:
            self.validate_fastq_in_sample_dir(self.fastq_dir)
        except AssertionError:
            logging.error(f"ERROR: Could not validate FASTQ files for {self.sample_id} at {self.sample_dir}. Skipping")
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
    project_dir: Path

    runinfoxml: Path = None
    runparametersxml: Path = None
    interop_dir: Path = None
    interop_files: [Path] = None

    def __post_init__(self):
        self.log_dir = self.run_dir / 'Logs'
        self.properties_dir = self.run_dir / 'Properties'
        self.run_id = self.run_dir.name

        # Get samplesheet, read into df
        self.samplesheet = self.get_samplesheet()
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
        interop_dir = self.project_dir.parents[1] / 'Runs' / self.experiment_name / 'Files' / 'InterOp'
        if not interop_dir.is_dir():
            logging.warning(
                f"Could not locate InterOp data for {self.experiment_name} - confirm researcher shared Run on BaseSpace")
            interop_dir = None
        return interop_dir

    def get_sample_id_dict(self) -> dict:
        """
        Extracts the Sample IDs from a BaseMount run directory, creates a dict containing Sample_ID:Sample_dir links
        """
        sample_id_dict = {}
        for sample_dir in self.sample_dirs:
            sample_id = sample_dir.name.split(".", 2)[2]
            sample_id_dict[sample_id] = sample_dir
        return sample_id_dict

    def get_samplesheet(self) -> Path:
        """
        Grabs and verifies the SampleSheet for a Run
        """
        samplesheet = self.properties_dir / "Input.sample-sheet"
        if samplesheet.is_file():
            return samplesheet
        else:
            raise FileNotFoundError(f"Could not find SampleSheet at expected location: {samplesheet}")

    def generate_sample_objects(self) -> [BaseMountSample]:
        """
        Generates a list of BaseMountSample data objects belonging to the respective Run
        """
        logging.debug(f"Generating BaseMountSample objects for {self.run_id} ({self.experiment_name})")
        sample_object_list = []
        for sample_id, sample_dir in self.sample_id_dict.items():
            sample_object = BaseMountSample(project_dir=self.project_dir,
                                            run_dir=self.run_dir,
                                            sample_dir=sample_dir,
                                            sample_id=sample_id)
            sample_object_list.append(sample_object)
        return sample_object_list

    def get_sample_dirs(self) -> list:
        """
        Grabs all Samples directories for the Run, filters out junk
        """
        sample_dirs = list(self.run_dir.glob("Sample.*"))

        # Filter out junk
        sample_dirs = [sample_dir for sample_dir in sample_dirs if 'Undetermined' not in sample_dir.name]
        sample_dirs = [sample_dir for sample_dir in sample_dirs if sample_dir.is_dir()]

        return sample_dirs

    def get_runparametersxml(self) -> Path:
        """
        # TODO: Check if this file exists anywhere else on BaseMount
        Tries to grab the RunParameters.xml file if its present in the expected location
        """
        runparametersxml = self.properties_dir / 'Input.Runs' / '0' / 'Files' / 'RunParameters.xml'
        if not runparametersxml.is_file():
            logging.warning(f"Could not locate RunParameters.xml for {self.experiment_name}")
            return None
        else:
            return runparametersxml

    def get_runinfoxml(self) -> Path:
        """
        Grabs the RunInfo.xml file
        """

        # Try to find RunInfo.xml in both known locations
        runinfoxml_1 = self.properties_dir / 'Input.Runs' / '0' / 'Files' / 'RunInfo.xml'
        runinfoxml_2 = self.run_dir / 'Logs' / 'RunInfo.xml'
        if runinfoxml_1.is_file():
            return runinfoxml_1
        elif runinfoxml_2.is_file():
            return runinfoxml_2
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


def retrieve_project_contents_from_basemount(project_dir: Path, out_dir: Path):
    """
    Main method to analyze BaseMount folder contents, establish dataclasses (Project, Run, Sample), and copy to out_dir
    """

    # Create output directory if it doesn't already exist
    out_dir.mkdir(exist_ok=True)

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

        # Copy run metadata over to outdir folder
        logging.debug(f"Copying metadata for {run_obj.run_id}")
        shutil.copy(str(run_obj.samplesheet), str(run_dir_out / 'SampleSheet.csv'))
        shutil.copy(str(run_obj.runinfoxml), str(run_dir_out / 'RunInfo.xml'))
        os.chmod(str(run_dir_out / 'RunInfo.xml'), 0o775)
        if run_obj.runparametersxml is not None:
            shutil.copy(str(run_obj.runparametersxml), str(run_dir_out / 'RunParameters.xml'))
            os.chmod(str(run_dir_out / 'RunParameters.xml'), 0o775)

        # Copy log files
        logging.debug(f"Copying log file contents for {run_obj.run_id}")
        for logfile in run_obj.logfiles:
            shutil.copy(str(logfile), str(run_dir_out / 'Logs' / logfile.name))
            os.chmod(str(run_dir_out / 'Logs' / logfile.name), 0o775)

        # Copy InterOp contents
        if run_obj.interop_files is not None:
            logging.debug(f"Copying InterOp contents for {run_obj.run_id}")
            for interop_file in run_obj.interop_files:
                interop_file_out = run_dir_out / 'InterOp' / interop_file.name
                shutil.copy(str(interop_file), str(interop_file_out))
                os.chmod(str(interop_file_out), 0o775)
        else:
            logging.debug(f"InterOp files not available for {run_obj.run_id}, skipping")

        # Copy reads over
        logging.info(f"Transferring reads for {run_obj.run_id} ({run_obj.experiment_name})...")
        for sample_obj in tqdm(iterable=run_obj.sample_objects, miniters=1):
            r1_out = run_dir_out / 'Data' / 'Intensities' / 'BaseCalls' / sample_obj.r1.name
            r2_out = run_dir_out / 'Data' / 'Intensities' / 'BaseCalls' / sample_obj.r2.name
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
@click.option('-p', '--project_dir',
              type=click.Path(exists=True),
              required=True,
              default=None,
              help='Path to the directory on BaseMount for a particular project. e.g. '
                   'basemount/Projects/[your project].',
              callback=convert_to_path)
@click.option('-o', '--out_dir',
              type=click.Path(exists=False),
              required=True,
              default=None,
              help='Directory to dump all runs for project.',
              callback=convert_to_path)
@click.option('-v', '--verbose',
              help='Specify this flag to enable more verbose output.',
              is_flag=True,
              default=False)
@click.option('--version',
              help='Specify this flag to print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              expose_value=False)
def cli(project_dir, out_dir, verbose):
    logging.info(f"Started BaseMountRetrieve (v{__version__})")

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug("Enabled VERBOSE mode")

    logging.debug(f"project_dir:\t\t{project_dir}")
    logging.debug(f"out_dir:\t\t{out_dir}")

    retrieve_project_contents_from_basemount(project_dir=project_dir, out_dir=out_dir)

    logging.info("Done!")


if __name__ == "__main__":
    cli()

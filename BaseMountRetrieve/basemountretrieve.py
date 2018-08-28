#!/usr/bin/env python3

__version__ = "0.2.4"
__author__ = "Forest Dussault"
__email__ = "forest.dussault@canada.ca"

import os
import click
import shutil
import logging
import pandas as pd
from pathlib import Path

script = os.path.basename(__file__)
logger = logging.getLogger()
logging.basicConfig(
    format=f'\033[92m \033[1m {script}:\033[0m %(message)s ',
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


@click.command()
@click.option('-p', '--projectdir',
              type=click.Path(exists=True),
              required=False,
              default=None,
              help='Path to the directory on BaseMount for a particular project. e.g. '
                   'basemount/Projects/[your project].',
              callback=convert_to_path)
@click.option('-o', '--outdir',
              type=click.Path(exists=False),
              required=True,
              default=None,
              help='Directory to dump all .fastq.gz files. Note that the Sample ID will be appended to the beginning '
                   'of the copied .fastq.gz file, which normally only contains the Sample Name.',
              callback=convert_to_path)
@click.option('--miseqsim',
              help='Specify this flag to simulate the MiSeq folder structure when retrieving from BaseSpace',
              is_flag=True,
              default=False)
@click.option('--version',
              help='Specify this flag to print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              expose_value=False)
def cli(projectdir, outdir, miseqsim):
    logging.info("Started BaseMountRetrieve")

    # Create output directory if it doesn't already exist
    os.makedirs(outdir, exist_ok=True)

    # Get samplesheets
    samplesheet_dict, run_translation_dict = retrieve_samplesheets(projectdir=projectdir, outdir=outdir)

    # Get list of samples to transfer, download files
    retrieve_samples(projectdir=projectdir, outdir=outdir, miseqsim=miseqsim)

    # Move everything around to simulate MiSeq folder structure
    if miseqsim:
        logfile_dict = retrieve_logfile_dict(projectdir)
        sample_dict = get_sample_dictionary(outdir)
        base_folders = ['Config',
                        'Data',
                        'Images',
                        'InterOp',
                        'Logs',
                        'Recipes',
                        'Thumbnail_Images']
        for run_id, samplesheet_path in samplesheet_dict.items():
            os.makedirs(outdir / run_id, exist_ok=True)
            shutil.copy(samplesheet_path, outdir / run_id / 'SampleSheet.csv')
            for f in base_folders:
                os.makedirs(outdir / run_id / f, exist_ok=True)

            read_folder = outdir / run_id / 'Data' / 'Intensities' / 'BaseCalls'
            os.makedirs(read_folder, exist_ok=True)

            df = read_samplesheet(samplesheet_path)
            sample_id_list = get_sample_id_list(df)
            for sample_id, reads in sample_dict.items():
                if sample_id in sample_id_list:
                    shutil.move(reads[0], read_folder / reads[0].name)
                    shutil.move(reads[1], read_folder / reads[1].name)

            # Copy all the log files over to the 'Logs' folder
            for verbose_run_name, log_list in logfile_dict.items():
                if run_translation_dict[verbose_run_name] == run_id:
                    for logfile in logfile_dict[verbose_run_name]:
                        outname = outdir / run_id / 'Logs' / logfile.name
                        shutil.copy(src=logfile, dst=outname)
                        os.chmod(str(outname), 0o775)  # Fix permissions

            retrieve_interop(run_id=run_id, projectdir=projectdir, outdir=outdir)

    # Delete remnant .csv files
    cleanup_csv = list(outdir.glob("*.csv"))
    for f in cleanup_csv:
        os.remove(f)

    logging.info(f"Process complete. Results available in {outdir}")


def retrieve_samples(projectdir: Path, outdir: Path, miseqsim: bool):
    # Gather all BaseMount file paths
    fastq_list = list(projectdir.glob("Samples/*/Files/*"))

    # Filter out hidden stuff
    fastq_list = [Path(x) for x in fastq_list if ".id." not in str(x)]

    # Prepare to copy files
    transfer_list = []

    if not miseqsim:
        outdir_files = list(outdir.glob('*'))
    else:
        outdir_files = list(outdir.glob('*/Data/Intensities/BaseCalls/*'))
    if len(outdir_files) == 0:
        transfer_list = fastq_list

    for i in fastq_list:
        # Get name components of sample
        sampleid = i.parents[1].name
        samplename = i.name

        # Prepare outfile names
        if sampleid not in samplename:
            outname = outdir / Path(sampleid + "_" + samplename)
        else:
            outname = outdir / Path(samplename)

        # Skip files that already exist
        for j in outdir_files:
            if outname.name not in str(j):
                transfer_list.append(i)
            else:
                logging.info(f"Skipping {i.name}")

    # Begin copying files
    for i in sorted(set(transfer_list)):
        # Get name components of sample
        sampleid = i.parents[1].name
        samplename = i.name

        # Copy to outdir
        if sampleid not in samplename:
            outname = outdir / Path(sampleid + "_" + samplename)
        else:
            outname = outdir / Path(samplename)

        if miseqsim:
            outdir_file_names = [x.name for x in outdir_files]
            tmp_name = sampleid + '_' + i.name
            if not tmp_name in outdir_file_names:
                logging.info(f"Copying {samplename}...")
                try:
                    shutil.copy(i, outname)
                    os.chmod(str(outname), 0o775)  # Fix permissions
                except IsADirectoryError:
                    logging.warning(f"WARNING: Could not copy {i} because it's a directory' ")
        else:
            if not outname.exists():
                logging.info(f"{outname.name} already exists. Skipping.")
            else:
                logging.info(f"Copying {samplename}...")
                shutil.copy(i, outname)  # shutil.copy is filesystem agnostic, unlike shutil.move, os.rename
                os.chmod(str(outname), 0o775)  # Fix permissions


def retrieve_logfile_dict(projectdir: Path) -> dict:
    run_folders = get_run_folders(projectdir)
    log_dict = dict()
    for verbose_run_name in run_folders:
        logfiles = list(verbose_run_name.glob('Logs/*'))
        logfiles = [Path(x) for x in logfiles if not Path(x).name.startswith(".")]
        log_dict[verbose_run_name.name] = logfiles
    return log_dict


def retrieve_samplesheets(projectdir: Path, outdir: Path) -> tuple:
    """
    Returns two dictionaries due to a poor design decision. TODO: Make this whole thing less sloppy
    :param projectdir:
    :param outdir:
    :return:
    """
    # Locate samplesheets
    samplesheets = list(projectdir.glob('AppSessions.v1/*/Properties/Input.sample-sheet'))
    samplesheets = [Path(x) for x in samplesheets if ".id." not in str(x)]

    if len(samplesheets) == 0:
        logging.error('ERROR: Could not find samplesheets for project. Quitting.')
        quit()

    # Copy samplesheets into outdir
    samplesheet_dict = dict()
    run_translation_dict = dict()
    for samplesheet in samplesheets:
        verbose_run_name = samplesheet.parents[1].name
        outname = outdir / (verbose_run_name + '.' + 'SampleSheet.csv')
        logging.info(f'Copying SampleSheet.csv for {samplesheet.parents[1].name} to {outname}')
        shutil.copy(str(samplesheet), str(outname))
        run_id = extract_run_name(samplesheet=samplesheet)
        run_translation_dict[verbose_run_name] = run_id
        samplesheet_dict[run_id] = outname
    return samplesheet_dict, run_translation_dict


def get_run_folders(projectdir: Path) -> list:
    runfolders = list(projectdir.glob('AppSessions.v1/*'))
    runfolders = [Path(x) for x in runfolders if not Path(x).name.startswith(".")]
    return runfolders


def read_samplesheet(samplesheet: Path) -> pd.DataFrame:
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


def validate_samplesheet_header(header: list) -> bool:
    """
    Validates that column names match expected values
    :param header: List of column names
    :return: True if header meets all expected values, False if not
    """
    expected_header = [
        'Sample_ID',
        'Sample_Name',
        'Sample_Plate',
        'Sample_Well',
        'I7_Index_ID',
        'index',
        'I5_Index_ID',
        'index2',
        'Sample_Project',
        'Description'
    ]
    if not set(header) == set(expected_header):
        raise Exception(f"Provided header {header} does not match expected header {expected_header}")
    else:
        return True


def get_sample_id_list(samplesheet_df: pd.DataFrame) -> list:
    """
    Returns list of all SampleIDs in SampleSheet dataframe
    :param samplesheet_df: df returned from read_samplesheet()
    :return: list of all Sample IDs
    """
    sample_id_list = list(samplesheet_df['Sample_ID'])
    return sample_id_list


def group_by_project(samplesheet_df: pd.DataFrame) -> dict:
    """
    Groups samples by project extracted from SampleSheet.csv.
    :param samplesheet_df: df returned from read_samplesheet()
    :return: project dictionary (Keys are project names, values are lists of associated samples)
    """
    project_list = list(samplesheet_df.groupby(['Sample_Project']).groups.keys())
    project_dict = {}
    for project in project_list:
        project_dict[project] = list(samplesheet_df[samplesheet_df['Sample_Project'] == project]['Sample_ID'])
    return project_dict


def retrieve_interop(run_id: str, projectdir: Path, outdir: Path):
    logging.info(f"Copying InterOp folder contents for {run_id}...")
    try:
        interop_folder = projectdir.parents[1] / 'Runs' / run_id / 'Files' / 'InterOp'
    except:
        logging.error("ERROR: Couldn't retrieve InterOp contents.")
        return

    interop_folder_contents = list(interop_folder.glob("*"))
    for f in interop_folder_contents:
        outname = outdir / run_id / 'InterOp' / f.name
        if outname.exists():
            logging.info(f"Skipping {outname.name} for {run_id}")
        else:
            logging.info(f"Copying {f}...")
            shutil.copy(src=f, dst=outname)
            try:
                os.chmod(str(outname), 0o775)
            except PermissionError:
                pass


def extract_run_name(samplesheet: Path) -> str:
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


def retrieve_fastqgz(directory: Path) -> [Path]:
    """
    :param directory: Path to folder containing output from MiSeq run
    :return: LIST of all .fastq.gz files in directory
    """
    fastq_file_list = list(directory.glob("*.f*q*"))
    return fastq_file_list


def retrieve_sampleids(fastq_file_list: [Path]) -> list:
    """
    :param fastq_file_list: List of fastq.gz filepaths generated by retrieve_fastqgz()
    :return: List of Sample IDs
    """
    # Iterate through all of the fastq files and grab the sampleID, append to list
    sample_id_list = list()
    for f in fastq_file_list:
        sample_id = f.name.split('_')[0]
        sample_id_list.append(sample_id)

    # Get unique sample IDs
    sample_id_list = list(set(sample_id_list))
    return sample_id_list


def get_readpair(sample_id: str, fastq_file_list: [Path], forward_id: str = "_R1",
                 reverse_id: str = "_R2") -> (list, None):
    """
    :param sample_id: String of sample ID
    :param fastq_file_list: List of fastq.gz file paths generated by retrieve_fastqgz()
    :param forward_id: ID indicating forward read in filename (e.g. _R1)
    :param reverse_id: ID indicating reverse read in filename (e.g. _R2)
    :return: the absolute filepaths of R1 and R2 for a given sample ID
    """

    r1, r2 = None, None
    for f in fastq_file_list:
        if sample_id in f.name:
            if forward_id in f.name:
                r1 = f
            elif reverse_id in f.name:
                r2 = f
    if r1 is not None and r2 is not None:
        return [r1, r2]
    else:
        logging.info('Could not pair {}'.format(sample_id))
        return None


def populate_sample_dictionary(sample_id_list: list, fastq_file_list: [Path]) -> dict:
    """
    :param sample_id_list: List of unique Sample IDs generated by retrieve_sampleids()
    :param fastq_file_list: List of fastq.gz file paths generated by retrieve_fastqgz()
    :return: dictionary with each Sample ID as a key and the read pairs as values
    """
    # Find file pairs for each unique sample ID
    sample_dictionary = {}
    for sample_id in sample_id_list:
        read_pair = get_readpair(sample_id, fastq_file_list)
        if read_pair is not None:
            sample_dictionary[sample_id] = read_pair
        else:
            pass
    return sample_dictionary


def get_sample_dictionary(directory: Path) -> dict:
    """
    Creates a sample dictionary with unique/valid sample IDs as keys and paths to forward and reverse reads as values
    :param directory: Path to a directory containing .fastq.gz files
    :return: Validated sample dictionary with sample_ID:R1,R2 structure
    """
    fastq_file_list = retrieve_fastqgz(directory)
    sample_id_list = retrieve_sampleids(fastq_file_list)
    sample_dictionary = populate_sample_dictionary(sample_id_list, fastq_file_list)
    if len(sample_dictionary) > 0:
        logging.info(f"Successfully paired {len(sample_dictionary)} of {len(sample_id_list)} samples")
    return sample_dictionary


if __name__ == "__main__":
    cli()

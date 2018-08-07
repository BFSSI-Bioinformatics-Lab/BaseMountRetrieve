#!/usr/bin/env python3

__version__ = "0.1.0"
__author__ = "Forest Dussault"
__email__ = "forest.dussault@canada.ca"

import os
import click
import shutil
import logging
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
@click.option('--version',
              help='Specify this flag to print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              expose_value=False)
def cli(projectdir, outdir):
    logging.info("Started BaseMountRetrieve")

    # Create output directory if it doesn't already exist
    os.makedirs(outdir, exist_ok=True)

    # Get samplesheets
    retrieve_samplesheets(projectdir=projectdir, outdir=outdir)

    # Get list of samples to transfer
    retrieve_samples(projectdir=projectdir, outdir=outdir)

    logging.info(f"Process complete. Results available in {outdir}")


def retrieve_samples(projectdir: Path, outdir: Path):
    # Gather all BaseMount file paths
    fastq_list = list(projectdir.glob("Samples/*/Files/*"))

    # Filter out hidden stuff
    fastq_list = [Path(x) for x in fastq_list if ".id." not in str(x)]

    # Avoid copying files that already exist
    transfer_list = []
    outdir_files = list(outdir.glob('*'))
    for j in fastq_list:
        # Get name components of sample
        sampleid = j.parents[1].name
        samplename = j.name

        # Copy to outdir and prepare file for chmod
        if sampleid not in samplename:
            outname = outdir / Path(sampleid + "_" + samplename)
        else:
            outname = outdir / Path(samplename)

        if outname.name not in outdir_files:
            transfer_list.append(j)
        else:
            logging.info(f"Skipping {j.name} (already present in {outdir})")

        modlist = []
        for i in transfer_list:
            logging.info(f"Copying {i.name}...")

            # Get name components of sample
            sampleid = i.parents[1].name
            samplename = i.name

            # Copy to outdir and prepare file for chmod
            if sampleid not in samplename:
                outname = outdir / Path(sampleid + "_" + samplename)
            else:
                outname = outdir / Path(samplename)
            shutil.copy(i,
                        outname)  # shutil.copy is filesystem agnostic, unlike shutil.move, os.rename, or Path.rename
            modlist.append(outname)

        # Fix permissions
        for f in modlist:
            os.chmod(str(f), 0o775)


def retrieve_samplesheets(projectdir: Path, outdir: Path):
    # Locate samplesheet
    try:
        samplesheets = list(projectdir.glob('AppSessions.v1/*/Properties/Input.sample-sheet'))
    except IndexError:
        print('ERROR: Could not find samplesheets for project')
        return None

    # Create copies
    for samplesheet in samplesheets:
        print(samplesheet)
        outname = outdir / Path(samplesheet.parents[1].name + '.' + 'SampleSheet.csv')
        shutil.copy(str(samplesheet), str(outname))


if __name__ == "__main__":
    cli()

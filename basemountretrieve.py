#!/usr/bin/env python3

__version__ = "0.0.1"
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


@click.command()
@click.option('-i', '--inputdir',
              type=click.Path(exists=True),
              required=False,
              default=None,
              help='Path to the Samples directory on BaseMount for a particular project. e.g. '
                   'basemount/Projects/[your project]/Samples. ')
@click.option('-o', '--outdir',
              type=click.Path(exists=False),
              required=False,
              default=None,
              help='Directory to dump all .fastq.gz files. Note that the Sample ID will be appended to the beginning '
                   'of the copied .fastq.gz file, which normally only contains the Sample Name.')
@click.option('--version',
              help='Specify this flag to print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              expose_value=False)
def cli(inputdir, outdir):
    inputdir = Path(inputdir)
    logging.info("Started BaseMountRetrieve")

    # Create output directory if it doesn't already exist
    os.makedirs(outdir, exist_ok=True)

    # Gather all BaseMount file paths
    fastq_list = inputdir.glob("*/Files/*")

    # Filter out hidden stuff
    fastq_list = [x for x in fastq_list if ".id." not in str(x)]

    modlist = []
    for i in fastq_list:
        logging.info(f"Copying {i.name}...")

        # Get name components of sample
        sampleid = i.parents[1].name
        samplename = i.name

        # Copy to outdir and prepare file for chmod
        outname = outdir / Path(sampleid + "_" + samplename)
        shutil.copy(i, outname)  # shutil.copy is filesystem agnostic, unlike shutil.move, os.rename, or Path.rename
        modlist.append(outname)

    # Fix permissions
    for f in modlist:
        os.chmod(str(f), 0o775)


if __name__ == "__main__":
    cli()

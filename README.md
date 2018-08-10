# BaseMountRetrieve
Internal BFSSI package for retrieving files from BaseMount.

### Requirements
- [BaseMount](https://basemount.basespace.illumina.com/)
- Python 3.6


### Usage
```bash
Usage: basemountretrieve.py [OPTIONS]

Options:
  -p, --projectdir PATH  Path to the directory on BaseMount for a particular
                         project. e.g. basemount/Projects/[your project].
  -o, --outdir PATH      Directory to dump all .fastq.gz files. Note that the
                         Sample ID will be appended to the beginning of the
                         copied .fastq.gz file, which normally only contains
                         the Sample Name.  [required]
  --miseqsim             Specify this flag to simulate the MiSeq folder
                         structure when retrieving from BaseSpace
  --version              Specify this flag to print the version and exit.
  --help                 Show this message and exit.

```
# BaseMountRetrieve
Package for retrieving files from BaseMount in the output style of a local MiSeq run.

### Output Directory Sample
```
Project
├── 20180714_WGS_M01308
│   ├── Config
│   ├── Data
│   │   └── Intensities
│   │       └── BaseCalls
│   ├── Images
│   ├── InterOp
│   ├── Logs
│   ├── Recipes
│   └── Thumbnail_Images
└── 20181102_WGS_M01308
    ├── Config
    ├── Data
    │   └── Intensities
    │       └── BaseCalls
    ├── Images
    ├── InterOp
    ├── Logs
    ├── Recipes
    └── Thumbnail_Images
```

### Requirements
- [BaseMount](https://basemount.basespace.illumina.com/) (Confirmed supported version: Illumina BaseMount v0.15.96.2154)
- Python 3.6

### Installation
`pip install basemountretrieve`

### Usage
```
Usage: basemountretrieve [OPTIONS]

  BaseMountRetrieve will tap into the mounted BaseMount filesystem and
  retrieve an entire Project (or single Run) in the output style of a
  local MiSeq run.

Options:
  -p, --project-dir PATH      Path to the directory on BaseMount for a
                              particular project. e.g.
                              basemount/Projects/[your project].
  -e, --experiment-name TEXT  Searches BaseMount for the name of a
                              run/experiment and attempts to retrieve its
                              contents. Must be used along with the
                              --basemount-dir parameter. Cannot be used
                              alongside the --project-dir flag.
  -b, --basemount-dir PATH    Path to root directory for BaseMount, e.g.
                              ~/basemount/. Must be supplied alongside
                              --experiment-name parameter.
  -o, --out-dir PATH          Directory to dump all runs for project.
                              [required]
  -r, --rename                Use this flag to automatically re-name the R1
                              and R2 files to just include the Sample ID.
  -v, --verbose               Use this flag to enable more verbose output.
  --version                   Use this flag to print the version and exit.
  --help                      Show this message and exit.
```
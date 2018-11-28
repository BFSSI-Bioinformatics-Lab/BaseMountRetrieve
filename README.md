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
  retrieve all of the runs for a given project in the output style of a
  local MiSeq run.

Options:
  -p, --project_dir PATH  Path to the directory on BaseMount for a particular
                          project. e.g. basemount/Projects/[your project]. [required]
  -o, --out_dir PATH      Directory to dump all runs for project.  [required]
  -v, --verbose           Specify this flag to enable more verbose output.
  --version               Specify this flag to print the version and exit.
  --help                  Show this message and exit.
```
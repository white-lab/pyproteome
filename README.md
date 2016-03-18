# pyproteome

[![Build Status](https://img.shields.io/travis/white-lab/pyproteome.svg)](https://travis-ci.org/white-lab/pyproteome)
[![Coverage Status](https://img.shields.io/coveralls/white-lab/pyproteome.svg)](https://coveralls.io/r/white-lab/pyproteome?branch=master)
[![Documentation Status](https://readthedocs.org/projects/pyproteome/badge/?version=latest)](https://readthedocs.org/projects/pyproteome/?badge=latest)
[![Requirements Status](https://requires.io/github/white-lab/pyproteome/requirements.svg?branch=master)](https://requires.io/github/white-lab/pyproteome/requirements/?branch=master)


Python library for analyzing mass spectrometry proteomics data.

## Installation

To install the core pyproteome package, run the following command:

```
pip install pyproteome
```

### CAMV

pyproteome can use CAMV for data validation. If you have the executable
installed on your system, simply add "CAMV.exe" to your system path and
pyproteome will locate it automatically.


## Directory Hierarchy

pyproteome expects a certain directory hierarchy in order to import data files
and interface with CAMV. This pattern is as follows:

```
base_directory/
    BCA Protein Assays/
    CAMV Output/
    CAMV Sessions/
    Mascot XMLs/
    MS RAW/
    MS Searched/
    Scan Lists/
    Scripts/
```

Under this scheme, all of your python code / IPython notebooks should go in the
`Scripts` directory.

See `pyproteome.paths` if you are using a custom directory hierarchy. i.e.:

```
>>> from pyproteome import paths
>>> paths.CAMV_SESS_DIR = "../CAMV Save/"
>>> paths.BCA_ASSAY_DIR = "../BCA/"
```

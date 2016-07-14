# pyproteome

[![Build Status](https://img.shields.io/travis/white-lab/pyproteome.svg)](https://travis-ci.org/white-lab/pyproteome)
[![Coverage Status](https://img.shields.io/coveralls/white-lab/pyproteome.svg)](https://coveralls.io/r/white-lab/pyproteome?branch=master)
[![Documentation Status](https://readthedocs.org/projects/pyproteome/badge/?version=latest)](https://readthedocs.org/projects/pyproteome/?badge=latest)
[![Requirements Status](https://requires.io/github/white-lab/pyproteome/requirements.svg?branch=master)](https://requires.io/github/white-lab/pyproteome/requirements/?branch=master)
[![PyPI](https://img.shields.io/pypi/v/pyproteome.svg)](https://pypi.python.org/pypi/pyproteome)


Python library for analyzing mass spectrometry proteomics data.

## Installation

To install the core pyproteome package, run the following command:

```
pip install pyproteome
```

### Windows

If you are using Windows, it is easiest to use the latest version of
[Anaconda](https://www.continuum.io/downloads) for your Python installation, as
pyproteome requires several hard-to-install packages, such as NumPy and SciPy.
In addition, BioPython should be installed from a [binary or wheel package](http://biopython.org/wiki/Download).

Then, you can simply run the above `pip install pyproteome` command to install
this package and the rest of its dependencies.

### CAMV

pyproteome can use CAMV for data validation. If you have the executable
installed on your system, simply add "CAMV.exe" to your system path and
pyproteome will locate it automatically.

## Examples

The following is an example of code to load a searched run from [Discoverer](https://www.thermofisher.com/order/catalog/product/IQLAAEGABSFAKJMAUH),
normalizing the phosphotyrosine run to the media channel levels in a supernatant
dilution.

```
>>> from pyproteome import analysis, data_sets, levels,
>>> from collections import OrderedDict
>>> ck_channels = OrderedDict([
...     ("126", "3130 CK"),
...     ("127", "3131 CK-p25"),
...     ("128", "3145 CK-p25"),
...     ("129", "3146 CK-p25"),
...     ("130", "3148 CK"),
...     ("131", "3157 CK"),
... ])
>>> ck_groups = OrderedDict([
...     ("CK-p25", ["127", "128", "129"]),
...     ("CK", ["126", "130", "131"]),
... ])
>>> ck_name = "CK-p25 vs. CK, 2 weeks"
>>> ck_h1_py = data_sets.DataSet(
...     mascot_name="2015-09-11-CKH1-pY-imac14-elute-pre35-colAaron250",
...     channels=ck_channels,
...     groups=ck_groups,
...     name="CKH1",
...     enrichments=["pY"],
...     tissues=["Hippocampus"],
... )
... ck_h1_global = data_sets.DataSet(
...     mascot_name="2015-09-18-CKH1-pY-2-sup-10-preRaven-colAaron250",
...     channels=ck_channels,
...     groups=ck_groups,
...     name="CKH1",
...     tissues=["Hippocampus"],
...     merge_duplicates=False,
...     merge_subsets=False,
... )
>>> ck_h1_channel_levels = levels.get_channel_levels(ck_h1_global.filter(ion_score_cutoff=20))
>>> ck_h1_py_norm = ck_h1_py.normalize(ck_h1_channel_levels)
>>> analysis.snr_table(ck_h1_py_norm.filter(p_cutoff=0.05), sort="Fold Change"))
```

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

# pyproteome

[![Build Status](https://img.shields.io/travis/white-lab/pyproteome.svg)](https://travis-ci.org/white-lab/pyproteome)
[![Coverage Status](https://img.shields.io/coveralls/white-lab/pyproteome.svg)](https://coveralls.io/r/white-lab/pyproteome?branch=master)
[![Documentation Status](https://readthedocs.org/projects/pyproteome/badge/?version=latest)](https://pyproteome.readthedocs.io/en/latest/)
[![Requirements Status](https://requires.io/github/white-lab/pyproteome/requirements.svg?branch=master)](https://requires.io/github/white-lab/pyproteome/requirements/?branch=master)
[![PyPI](https://img.shields.io/pypi/v/pyproteome.svg)](https://pypi.python.org/pypi/pyproteome)


Python library for analyzing mass spectrometry proteomics data.

## Installation

To install the core pyproteome package, run the following command:

```
pip install --process-dependency-links pyproteome
```

### Windows

If you are using Windows, it is easiest to use the latest version of
[Anaconda](https://www.continuum.io/downloads) for your Python installation, as
pyproteome requires several hard-to-install packages, such as NumPy and SciPy.

Then, you can simply run the above `pip install pyproteome` command to install
this package and the rest of its dependencies.

### CAMV

pyproteome can use CAMV for data validation. If you have the executable
installed on your system, simply add "CAMV.exe" to your system path and
pyproteome will locate it automatically.

## Examples

The following is an example of code to load a [ProteomeDiscoverer](https://www.thermofisher.com/order/catalog/product/IQLAAEGABSFAKJMAUH),
.msf data set, normalize the phosphotyrosine run to its corresponding
supernatant data set, and merge multiple runs into one final data set:

```
ckh_channels = OrderedDict(
    [
        ("3130 CK Hip",     "126"),
        ("3131 CK-p25 Hip", "127"),
        ("3145 CK-p25 Hip", "128"),
        ("3146 CK-p25 Hip", "129"),
        ("3148 CK Hip",     "130"),
        ("3157 CK Hip",     "131"),
    ]
)
ckx_channels = OrderedDict(
    [
        ("3130 CK Cortex",     "126"),
        ("3131 CK-p25 Cortex", "127"),
        ("3145 CK-p25 Cortex", "128"),
        ("3146 CK-p25 Cortex", "129"),
        ("3148 CK Cortex",     "130"),
        ("3157 CK Cortex",     "131"),
    ]
)
ckp25_groups = OrderedDict(
    [
        (
            "CK",
            [
                "3130 CK Hip",
                "3148 CK Hip",
                "3157 CK Hip",
                "3130 CK Cortex",
                "3148 CK Cortex",
                "3157 CK Cortex",
            ],
        ),
        (
            "CK-p25",
            [
                "3131 CK-p25 Hip",
                "3145 CK-p25 Hip",
                "3146 CK-p25 Hip",
                "3131 CK-p25 Cortex",
                "3145 CK-p25 Cortex",
                "3146 CK-p25 Cortex",
            ],
        ),
    ]
)
# With search data located as follows:
#   Searched/
#       CK-H1-pY.msf
#       CK-H1-pST.msf
#       CK-H1-Global.msf
#       CK-X1-pY.msf
#       CK-X1-pST.msf
#       CK-X1-Global.msf
datas = data_sets.load_all_data(
    chan_mapping={
        "CK-H": ckh_channels,
        "CK-X": ckx_channels,
    },
    # Normalize pY, pST, and Global runs to each sample's global data
    norm_mapping=OrderedDict([
        ("CK-H1", "CK-H1-Global"),
        ("CK-X1", "CK-X1-Global"),
    ]),
    # Merge together normalized hippocampus and cortex runs
    merge_mapping=OrderedDict([
        ("CK Hip", ["CK-H1-pY", "CK-H1-pST", "CK-H1-Global"]),
        ("CK Cortex", ["CK-X1-pY", "CK-X1-pST", "CK-X1-Global"]),
        ("CK All", ["CK Hip", "CK Cortex"]),
    ]),
    groups=ckp25_groups,
)
```

For other functionality, refer to the
[online documentation](https://pyproteome.readthedocs.io/en/latest/). There
are also a set of example analyses located in the [pyproteome-data
repository](https://github.com/white-lab/pyproteome-data/tree/master/examples).

## Directory Hierarchy

pyproteome expects a certain directory hierarchy in order to import data files
and interface with CAMV. This pattern is as follows:

```
base_directory/
    BCA Protein Assays/
    CAMV Output/
    CAMV Sessions/
    Figures/
    MS RAW/
    Scripts/
    Searched/
```

Under this scheme, all of your python code / IPython notebooks should go in the
`Scripts` directory.

See `pyproteome.paths` if you are using a custom directory hierarchy. i.e.:

```
>>> from pyproteome import paths
>>> paths.CAMV_SESS_DIR = "../CAMV Save/"
>>> paths.BCA_ASSAY_DIR = "../BCA/"
```

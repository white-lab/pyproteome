
import os
from setuptools import setup, find_packages

from pyproteome import version


REQUIREMENTS_PATH = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__), "requirements.txt",
    )
)


setup(
    name="pyproteome",
    version=version.__version__,
    description="Package for managing proteomics data",
    url="https://github.com/white-lab/pyproteome",
    author="Nader Morshed",
    author_email="morshed@mit.edu",
    license="BSD",
    packages=find_packages(exclude=["*.tests", "tests"]),
    install_requires=[
        "adjustText>=0.6.1",
        "fastcluster>=1.1.23",
        "genemap",
        "ipython>=5.4.1",
        "matplotlib>=2.2.0",
        "numpy>=1.13.0",
        "numpydoc>=0.7",
        "pandas>=0.22.0",
        "scikit-learn>=0.19.1",
        "scipy>=1.0.0",
        "seaborn>=0.8.1",
        "uniprot==1.3",
        "xlrd>=1.1.0",
        "xlwt>=1.3.0",
        "xlsxwriter>=1.0.2",
    ],
    dependency_links=[
        "git+https://github.com/jrderuiter/genemap.git#egg=genemap",
    ],
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
    ],
    test_suite="tests",
)

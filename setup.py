
from setuptools import setup

from pyproteome import __version__


setup(
    name="pyproteome",
    version=__version__,
    description="Package for managing proteomics data",
    url="https://github.com/naderm/pyproteome",
    author="Nader Morshed",
    author_email="morshed@mit.edu",
    license="BSD",
    packages=["pyproteome"],
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.2",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering",
    ],
    test_suite="tests",
)


import os
from pip.req import parse_requirements
from pip.download import PipSession
from setuptools import setup, find_packages

from pyproteome import __version__


REQUIREMENTS_PATH = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__), "requirements.txt",
    )
)


setup(
    name="pyproteome",
    version=__version__,
    description="Package for managing proteomics data",
    url="https://github.com/white-lab/pyproteome",
    author="Nader Morshed",
    author_email="morshed@mit.edu",
    license="BSD",
    packages=find_packages(exclude=["*.tests", "tests"]),
    install_requires=[
        str(i.req)
        for i in parse_requirements(REQUIREMENTS_PATH, session=PipSession())
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

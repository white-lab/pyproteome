
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
)


import os
from setuptools import setup, find_packages

__dir__ = os.path.dirname(__file__)

with open(
    os.path.join(__dir__, 'pyproteome', 'version.py')
) as f:
    __version__ = '0.0.0'

    for line in f:
        if '#' in line:
            line = line[:line.index('#')]

        if not line.startswith('version ='):
            continue

        __version__ = line.split('=')[1].strip().strip('\'').strip('\"')

REQUIREMENTS = [
    'adjustText>=0.7.3',
    'fastcluster>=1.1.25',
    'goatools>=0.9',
    'ipython>=5.4.1',
    'matplotlib>=3.0',
    'numpy>=1.17',
    'numpydoc>=0.9',
    'pandas>=0.23',
    'pingouin>=0.3',
    'scikit-learn>=0.21',
    'scipy>=1.3.0',
    'seaborn>=0.9.0',
    'uniprot==1.3',
    'xlrd>=1.2',
    'xlwt>=1.3',
    'xlsxwriter>=1.2',
]

if __name__ == '__main__':
    setup(
        name='pyproteome',
        version=__version__,
        description='Package for managing proteomics data',
        url='https://github.com/white-lab/pyproteome',
        author='Nader Morshed',
        author_email='morshed@mit.edu',
        license='BSD-2-Clause',
        packages=find_packages(exclude=['*.tests', 'tests']),
        install_requires=REQUIREMENTS,
        extras_require={
            'photon': [
                'photon_ptm '
                '@ http://github.com/jdrudolph/photon/archive/master.zip',
                'perseuspy>=0.3.8',
            ],
            'genemap': [
                'genemap '
                '@ http://github.com/naderm/genemap/archive/master.zip',
            ],
        },
        classifiers=[
            'License :: OSI Approved :: BSD License',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Topic :: Scientific/Engineering',
        ],
        test_suite='tests',
    )

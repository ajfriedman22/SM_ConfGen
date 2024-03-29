"""
SM_ConfGen
A repository for run temperature replica exchange MD on small molecules
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = "A package for setting up, performing, and analyzing temperature replica exchange molecular dynaimcs simulations using GROMACS for the purpose of small molecule conformer generation".split("\n")[0]

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except FileNotFoundError:
    long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name='SM_ConfGen',
    author='Anika Friedman',
    author_email='Anika.Friedman@colorado.edu',
    description=short_description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',
    project_urls={
        "Source Code": "https://https://github.com/ajfriedman22/SM_ConfGen",
    },
    keywords="molecular mechanics",

    # Describes the project using a list of classifiers (https://pypi.org/classifiers/)
    # This makes the project more searchable.
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows ",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # Add entry points
    entry_points={
        'console_scripts':[
            'run_ConfGen = SM_ConfGen.cli.run_ConfGen:main',
            'analyze_ConfGen = SM_ConfGen.cli.analyze_ConfGen:main',
        ],
    },

    # Required packages, pulls from pip if needed; do not use for Conda deployment
    install_requires=[
        'numpy',
        'argparse',
        'pyyaml',
        'seaborn',
        'matplotlib',
        'mpi4py',
    ],

    platforms=['Linux',
               'Mac OS-X',
               'Unix',
               'Windows'],            # Valid platforms your code works on, adjust to your flavor
    
    python_requires=">=3.8",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
"""This is a setup.py script to install ShakeNBreak."""

import os
import warnings

from setuptools import find_packages, setup
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info
from setuptools.command.install import install

path_to_file = os.path.dirname(os.path.abspath(__file__))

class PostInstallCommand(install):
    """Post-installation for installation mode.

    Subclass of the setup tools install class in order to run custom commands
    after installation. Note that this only works when using 'python setup.py install'
    but not 'pip install .' or 'pip install -e .'.
    """

    def run(self):
        """
        Performs the usual install process and then copies the True Type fonts
        that come with SnB into matplotlib's True Type font directory,
        and deletes the matplotlib fontList.cache.
        """
        # Perform the usual install process
        install.run(self)
        
def package_files(directory):
    """Include package data."""
    paths = []
    for path, _dir, filenames in os.walk(directory):
        paths.extend(os.path.join("..", path, filename) for filename in filenames)
    return paths


input_files = package_files("input_files/")

with open("README.md", encoding="utf-8") as file:
    long_description = file.read()

setup(
    name="dftutils",
    version="1.0",
    description="Package to help with setting up calculations and post processing with VASP",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Alexandre Silva",
    author_email="alexandre.silva@fisica.uminho.pt",
    maintainer="Alexandre Silva",
    maintainer_email="alexandre.silva@fisica.uminho.pt",
    url="https://github.com/Alexandre-M-Silva/dftutils",
    license="MIT",
    license_files=("LICENSE",),
    classifiers=[
        "Development Status :: 1 - In Development",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="chemistry pymatgen dft defects structure-searching distortions symmetry-breaking",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy",  # >=1.21.2" needed for numpy.typing.NDArray?
        "pymatgen>=2022.10.22",
        "pymatgen-analysis-defects>=2022.10.28",
        "matplotlib>=3.6",
        "pandas>=1.1.0",
        "monty",
        "click>8.0",
        "importlib_metadata",
    ],
    # Specify any non-python files to be distributed with the package
    package_data={
        "dftutils": ["dftutils/*", *input_files],
    },
    include_package_data=True,
    # Specify the custom installation class
    zip_safe=False,
    entry_points={
        "console_scripts": [
            "dftutils = dftutils.cli:dftutils",
            "dftutils-polarization = dftutils.cli:polarization",
        ],
    },
    cmdclass={
        "install": PostInstallCommand,
    },
    project_urls={
        "Homepage": "https://github.com/Alexandre-M-Silva/dftutils",
    },
)
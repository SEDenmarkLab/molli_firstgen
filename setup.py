"""
    This module enables `pip install ./`
    Do not run this script standalone.
    Refer to package readme for further details.
"""

from setuptools import setup, find_packages
from glob import glob
import os

setup(
    name="molli",
    packages=[
        "molli",
    ],
    data_files=[],
    version="0.1.0-unstable",
    author="Alexander S. Shved",
    author_email="shvedalx@illinois.edu",
    install_requires=[
        "matplotlib>=3.1.2",
        "numpy>=1.19.0",
        "PyYAML>=5.3",
        "scikit-learn>=0.22.1",
        "scipy>=1.4.1",
        "colorama>=0.4.4",
    ],
    python_requires=">=3.9",
)

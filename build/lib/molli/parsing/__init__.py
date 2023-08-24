"""
    This package natively requires a lot of cross-talking between a lot of modules.
    Hence a generally structured package for parsing purposes is more than desirable.

    The idea is that this folder contains all parsers that Molecule class may potentially need
"""

from .cdxml import split_cdxml
from .xtbout import extract_xtb_atomic_properties

"""
    This module provides hooks to external computational packages
"""

# TODO add checks for package installations maybe?
# TODO add warnings for packages that have not been installed

from .obabel import OpenBabelDriver
from .xtb import XTBDriver, CRESTDriver
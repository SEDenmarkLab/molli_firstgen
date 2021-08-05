"""
    This module provides hooks to external computational packages
"""

# TODO add checks for package installations maybe?
# TODO add warnings for packages that have not been installed


from ._core import AsyncConcurrent as Concurrent

# from .obabel import OpenBabelDriver
from .xtb import AsyncXTBDriver as XTBDriver
from .crest import AsyncCRESTDriver as CRESTDriver
from .obabel import AsyncOpenBabelDriver as OpenBabelDriver
from .orca import AsyncORCADriver as ORCADriver



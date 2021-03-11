"""
    Provides a simple parser for .xyz files (and multi-xyz files)
"""
from typing import List
import numpy as np


def split_xyz(xyzblock: str) -> List[List[str]]:
    """
        Splits an multixyz block into individual xyz blocks
    """
    lines = xyzblock.splitlines(keepends=False)

    xyzblocks = []

    i = 0
    while i in range(len(lines)):
        bs = int(lines[i].strip())  # block size
        xyzblocks.append(lines[i:i + bs + 2])
        i += bs + 2

    return xyzblocks


def parse_xyz_lines(xyzlines: List[str]) -> (np.ndarray, List[str], str):
    """
    xyzlines: as parsed by split_xyz
    Returns a numpy array, list of atoms and comment string
    """
    _lines = [l.strip() for l in xyzlines]
    L = int(_lines[0])

    assert len(_lines) == L + 2

    comment = _lines[1]

    atoms = []
    coord = np.zeros((L, 3))

    for i, l in enumerate(_lines[2:]):
        a, x, y, z = l.split()
        atoms.append(a)
        coord[i] = [float(x), float(y), float(z)]

    return coord, atoms, comment


def parse_xyz(xyzblock: str, single: bool = True, assert_single: bool = False):
    """
        Note that single will return the first block always.
        Even for multixyz files.
        If checks need to be performed, use assert_single
        returns coord, atoms, comment
    """
    blocks = split_xyz(xyzblock)

    if not blocks:
        raise SyntaxError("No valid entries")

    if single and not assert_single:
        return parse_xyz_lines(blocks[0])
    elif single and assert_single:
        if len(blocks) > 1:
            raise AssertionError("Single-xyz check failed")
        else:
            return parse_xyz_lines(blocks[0])
    else:
        return list(map(parse_xyz_lines, blocks))


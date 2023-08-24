"""
    Superclass for all grid based descriptors
"""
from __future__ import annotations
import numpy as np
import re
import json

from ._core import Descriptor


class Grid:
    """
    This is a base class for molecular grid storage
    """

    def __init__(self, gpts: np.ndarray, dtype=np.float32, precision=4):
        """
        Grid is 'just' a list of ordered points
        """
        self.dtype = dtype
        self.precision = precision
        self.gridpoints = np.array(gpts, dtype=dtype)

    def __eq__(self, o: Grid) -> bool:
        if self.gridpoints.shape != o.gridpoints.shape:
            return False
        else:
            diff = self.gridpoints - o.gridpoints

            if np.max(diff) > 10 ** (-self.precision):
                return False
            else:
                return True

    def to_string(self):
        """
        Export grid contents as an importable string
        """
        gs = f"#{len(self.gridpoints)},3:"

        for x, y, z in self.gridpoints:
            p = self.precision
            gs += f"{x:0.{p}f},{y:0.{p}f},{z:0.{p}f};"

        return gs

    @classmethod
    def from_string(cls, s: str, dtype=np.float32, precision=4):
        """
        Import grid contents from a string
        """

        m = re.match(r"#(?P<L>[0-9]+),(?P<D>[0-9]+):(?P<G>.+);", s)

        L = int(m.group("L"))
        D = int(m.group("D"))
        G = m.group("G")

        assert D == 3, "Only 3d coordinates supported for now"

        coord = []

        for xyz in G.split(";"):
            x, y, z = xyz.split(",")
            coord.append([float(x), float(y), float(z)])

        assert len(coord) == L, f"Expected {L} coords, found {len(coord)}"

        return cls(coord, dtype=dtype, precision=precision)


class RectangularGrid(Grid):
    """
    This is an initialization class for a rectangular grid
    """

    def __init__(
        self, x: np.ndarray, y: np.ndarray, z: np.ndarray, dtype=np.float32, precision=4
    ):

        gpts = []

        for _x in x:
            for _y in y:
                for _z in z:
                    gpts.append([_x, _y, _z])

        super(RectangularGrid, self).__init__(gpts, dtype=dtype, precision=precision)


class GridDescriptor(Descriptor):
    def __init__(self, grid: Grid, values: np.ndarray = None):
        self.grid = grid
        if isinstance(values, np.ndarray):
            if values.shape == (grid.gridpoints.shape[0],):
                self.values = values
            else:
                raise ValueError("Inconsistent shape for values array. Double-check!")

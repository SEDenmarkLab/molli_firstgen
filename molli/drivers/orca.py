from ._core import AsyncExternalDriver
from ..dtypes import Atom, Bond, Molecule, CartesianGeometry
from copy import deepcopy
from datetime import datetime
from glob import glob
from warnings import warn
from typing import List, Callable
from math import ceil, pi
from itertools import combinations
import asyncio as aio
import numpy as np
import re

class AsyncORCADriver(AsyncExternalDriver):
    def __init__(self, name="", path="/opt/orca/orca_4_2_1_linux_x86-64_openmpi314/orca", scratch_dir="", nprocs=1, maxcore=3000, encoding="utf8"):
        super().__init__(
            name=name, scratch_dir=scratch_dir, nprocs=nprocs, encoding=encoding
        )
        self.path=path
        self.maxcore = maxcore
    
    async def xyz_energy(self, xyz: str, method="B97-3c sloppyscf"):
        """
            Calculate a single point energy for a given xyz structure
        """
        _inp = f"%pal nprocs {self.nprocs} end\n%maxcore {self.maxcore}\n!rks {method} energy nopop miniprint noprintmos\n*xyzfile 0 1 struct.xyz\n\n"

        code, files, stdout, stderr = await self.aexec(f"{self.path} input", inp_files={f"struct.xyz" : xyz, "input" : _inp})

        for l in stdout.split('\n')[::-1]:
            if m := re.match(r"FINAL SINGLE POINT ENERGY\s+(?P<eh>[0-9.-]+).*", l):
                return float(m['eh'])

    async def conformer_energies(self, mol: Molecule, method="B97-3c sloppyscf"):
        xyzs = mol.confs_to_xyzs()
        energies = []

        for i, xyz in enumerate(xyzs):
            conf_energy = await self.xyz_energy(xyz, method=method)
            energies.append(conf_energy)
        
        ref_energy = energies[0]
          
        return (np.array(energies) - ref_energy)*2625.5 # conversion to kJ/mol







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
import re

class AsyncORCADriver(AsyncExternalDriver):
    def __init__(self, name="", path="/opt/orca/orca_4_2_1_linux_x86-64_openmpi314/orca", scratch_dir="", nprocs=1, maxcore=3000, encoding="utf8"):
        super().__init__(
            name=name, scratch_dir=scratch_dir, nprocs=nprocs, encoding=encoding
        )
        self.path=path
        self.maxcore = maxcore

    
    async def conformer_energies(self, mol: Molecule, method="B97-3c"):
        energies = []
        nn = mol.name
        for i, xyz in enumerate(mol.confs_to_xyzs()):
            _inp = f"""
%pal nprocs {self.nprocs} end
%maxcore {self.maxcore}

!rks {method} energy sloppyscf nopop miniprint noprintmos

*xyzfile 0 1 {nn}.{i}.xyz

            """


        
            code, files, stdout, stderr = await self.aexec(f"{self.path} input", inp_files={"input" : _inp, f"{nn}.{i}.xyz" : xyz + "\n\n"})

            for l in stdout.split('\n')[::-1]:
                if m := re.match(r"FINAL SINGLE POINT ENERGY\s+(?P<eh>[0-9.-]+).*", l):
                    energies.append(float(m['eh']))
        
        return energies







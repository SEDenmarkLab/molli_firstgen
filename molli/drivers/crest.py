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


class AsyncCRESTDriver(AsyncExternalDriver):
    # def __init__(self, name="", scratch_dir="", nprocs=1, encoding="utf8"):
    #     super().__init__(
    #         name=name, scratch_dir=scratch_dir, nprocs=nprocs, encoding=encoding
    #     )

    async def conformer_search(
        self,
        mol: Molecule,
        method: str = "gff",
        ewin: float = 10,
        mdlen: float = 20,
        mddump: float = 0.250,
        vbdump: float = 1,
    ):
        """
        `ewin`: energy window in kcal/mol
        `mdlen`: iMTD-GC molecular dynamics length
        """
        g0_xyz = mol.to_xyz()

        nn = mol.name

        # command that will be used to execute xtb package
        # _cmd = f"""crest {nn}_g0.xyz -{method} -ewin {ewin:0.4f} -mdlen {mdlen:0.4f} -mddump {mddump:0.4f} -vbdump {vbdump:0.4f} -T {self.nprocs}"""

        _cmd = f"""crest {nn}_g0.xyz -{method} -ewin {ewin:0.4f} """

        code, files, stdout, stderr = await self.aexec(
            _cmd,
            inp_files={f"{nn}_g0.xyz": g0_xyz},
            out_files=["crest_conformers.xyz"],
        )

        try:
            ens1 = files["crest_conformers.xyz"]
        except Exception as xc:
            raise xc

        # screen the conformers
        _cmd1 = f"""crest -screen {nn}_ens1.xyz -{method} -ewin {ewin:0.4f} -T {self.nprocs}"""

        code, files, stdout, stderr = await self.aexec(
            _cmd1,
            inp_files={f"{nn}_ens1.xyz": ens1},
            out_files=["crest_ensemble.xyz"],
        )

        try:
            ens2 = files["crest_ensemble.xyz"]
        except Exception as xc:
            raise xc

        parsed_ens2 = CartesianGeometry.from_xyz(ens2)

        for geom, _, _ in parsed_ens2:
            mol.embed_conformers(geom, mode="a")

        return len(parsed_ens2)

        # code, files, stdout, stderr = await self.aexec(
        #     _cmd,
        #     inp_files={f"{nn}_g0.xyz": g0_xyz},
        #     out_files=["crest_conformers.xyz"],
        # )

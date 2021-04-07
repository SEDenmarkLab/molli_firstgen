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


class AsyncXTBDriver(AsyncExternalDriver):
    def __init__(self, name="", scratch_dir="", nprocs=1, encoding="utf8"):
        super().__init__(
            name=name, scratch_dir=scratch_dir, nprocs=nprocs, encoding=encoding
        )

    async def optimize(
        self,
        mol: Molecule,
        method: str = "gff",
        crit: str = "normal",
        xtbinp: str = "",
        maxiter: int = 50,
        in_place: bool = False,
    ):
        """
        Attempt a geometry optimization with parameters specified
        """

        g0_xyz = mol.to_xyz()

        nn = mol.name

        # command that will be used to execute xtb package
        _cmd = f"""xtb {nn}_g0.xyz --{method} --opt {crit} --cycles {maxiter} {"--input param.inp" if xtbinp else ""} -P {self.nprocs}"""

        # pylint: disable=unused-variable
        code, files, stdout, stderr = await self.aexec(
            _cmd,
            inp_files={f"{nn}_g0.xyz": g0_xyz, "param.inp": xtbinp},
            out_files=["xtbopt.xyz"],
        )

        if "xtbopt.xyz" in files:
            nxyz = files["xtbopt.xyz"]

            if not in_place:
                mol1 = deepcopy(mol)
                mol1.update_geom_from_xyz(nxyz, assert_single=True)
                return mol1
            else:
                mol.update_geom_from_xyz(nxyz, assert_single=True)
        else:
            raise FileNotFoundError("Could not locate xtb output file.")

    async def fix_long_bonds(
        self,
        mol: Molecule,
        method: str = "gff",
        rss_length_thresh: float = 4.0,
        rss_steps: float = 15,
        rss_maxcycle: int = 20,
        force_const: float = 0.5,
        target_len: float = 1.5,
        in_place: bool = False,
        constrain_bonds: list = ["C-C", "C-H", "C-F", "C-O", "O-H", "N-Cl"],
        fn_suffix: str = 0,
    ):
        """
        Fix all long bonds in the molecule by doing a relaxed surface scan with coordinates constrained
        """
        tbf = {}  # to be fixed
        for b in mol.bonds:
            if b not in tbf and mol.get_bond_length(b) >= rss_length_thresh:
                tbf[b] = mol.get_bond_length(b)

        if not tbf:
            return mol
        inp = f"$constrain\n  force constant={force_const}\n"

        lb_atoms = set()

        for b in tbf:
            a1, a2 = mol.get_atom_idx(b.a1), mol.get_atom_idx(b.a2)
            inp += f"  distance: {a1+1}, {a2+1}, {tbf[b]:0.4f}\n"
            lb_atoms.add(b.a1)
            lb_atoms.add(b.a2)

        # generate constraints for C-H bonds
        core_bonds = tuple(mol.yield_bonds(*constrain_bonds))
        inp += self.gen_bond_constraints(mol, core_bonds)
        inp += self.gen_angle_constraints(mol, lb_atoms)

        inp += "$scan\n  mode=concerted\n"
        for i, b in enumerate(tbf):
            inp += f"  {i+1}: {tbf[b]:0.4f}, {target_len:0.4f}, {rss_steps}\n"

        inp += f"$opt\n  maxcycle={rss_maxcycle}\n"
        inp += "$end\n"

        m1 = await self.optimize(
            mol,
            method=method,
            crit="crude",
            xtbinp=inp,
            in_place=False,
        )

        return m1

    @staticmethod
    def gen_bond_constraints(mol: Molecule, bonds: List[Bond]):
        """ Generate bond distance constraint list """
        constr = ""
        for b in bonds:
            a1, a2 = mol.get_atom_idx(b.a1), mol.get_atom_idx(b.a2)
            constr += f"  distance: {a1+1}, {a2+1}, {mol.get_bond_length(b):0.4f}\n"
        return constr

    @staticmethod
    def gen_angle_constraints(mol: Molecule, atoms: List[Atom]):
        """ Generate constraints for all angles where atom is the middle atom """
        constr = ""
        for a in atoms:
            neigbors = mol.get_connected_atoms(a)
            for a1, a2 in combinations(neigbors, 2):
                i1 = mol.get_atom_idx(a1) + 1
                i2 = mol.get_atom_idx(a) + 1
                i3 = mol.get_atom_idx(a2) + 1
                angle = mol.get_angle(a1, a, a2) * 180 / pi
                constr += f"  angle: {i1}, {i2}, {i3}, {angle:0.4f}\n"

        return constr

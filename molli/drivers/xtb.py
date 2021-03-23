"""
    This module provides interactions with XTB package
"""
import os
from ._core import ExternalDriver, DriverError
from ..dtypes import Molecule, CartesianGeometry
from copy import deepcopy
from datetime import datetime
from glob import glob
from warnings import warn
from typing import List, Callable
from math import ceil


class XTBDriver(ExternalDriver):
    """
    This driver provides functionality for
    --- XTB package (Grimme group) ---
    https://xtb-docs.readthedocs.io/en/latest/setup.html
    """

    PREFIX = "mx"
    JOB_ID = 0

    def __init__(
        self,
        cwd: str = "/temp_xtb/",
        nprocs: int = 1,
        method: str = "gfn2",
        accuracy: float = 1.0,
        opt_maxiter: int = 100,
        opt_crit: str = "normal",  # crude, sloppy, normal, tight, vtight
    ):
        super().__init__(cwd=cwd, nprocs=nprocs)

        self.method = method
        self.accuracy = accuracy
        self.opt_maxiter = opt_maxiter
        self.opt_crit = opt_crit

    def __call__(self, *args):
        print("xtb", *args, "-P", self.nprocs)
        super().__call__("xtb", *args, "-P", self.nprocs)

    def gen_constraints(
        self,
        mol: Molecule,
        fc: float = 0.5,
        maxcycle: int = 10,
        distances: list = [],
        angles: list = [],
        dihedrals: list = [],
        scans: list = [],
        concerted=True,
    ):
        """
        Generate constraint file
        distances = [(a1, a2, val1), (a3, a4, val2), ...]
        angles = [(a1, a2, a3, val1), ...]
        scan = [(idx, val_start, val_end, steps)]
        """
        result = f"$constrain\n  force constant= {fc}\n"
        for a1, a2, val in distances:
            i1, i2 = mol.get_atom_idx(a1), mol.get_atom_idx(a2)
            result += f"  distance: {i1+1}, {i2+1}, {val}\n"

        for a1, a2, a3, val in angles:
            i1, i2, i3 = (
                mol.get_atom_idx(a1),
                mol.get_atom_idx(a2),
                mol.get_atom_idx(a3),
            )
            result += f"  angle: {i1+1}, {i2+1}, {i3+1}, {val}\n"

        for a1, a2, a3, a4, val in dihedrals:
            i1, i2, i3, i4 = (
                mol.get_atom_idx(a1),
                mol.get_atom_idx(a2),
                mol.get_atom_idx(a3),
                mol.get_atom_idx(a4),
            )
            result += f"  dihedral: {i1+1}, {i2+1}, {i3+1}, {i4+1}, {val}\n"

        if scans:
            result += "$scan\n"
            if concerted and len(scans) > 1:
                result += "  mode = concerted\n"
            for idx, start, end, steps in scans:
                result += f"  {idx}: {start}, {end}, {steps}\n"

        result += f"$opt\n  maxcycle={maxcycle}\n$end\n"

        return result

    def optimize(
        self,
        mol: Molecule,
        crit: str = None,  # overrides criterion in the method, if not none
        in_place: bool = False,
        xtbinp: str = "",
        fn_suffix=0,
    ):
        """
        Attempt a geometry optimization with parameters from the instance.

        """
        cmd = []  # collection of command line arguments for xtb binary
        jobid = self.__class__.JOB_ID
        name = f"{self.PREFIX}.opt.{os.getpid()}.{jobid}.{fn_suffix}"
        with open(f"{self.cwd}/{name}.xyz", "wt") as xyzf:
            xyzf.write(mol.to_xyz())

        if crit != None:
            _crit = crit
        else:
            _crit = self.opt_crit

        cmd.extend(
            (
                f"{name}.xyz",
                self.method,
                "--opt",
                _crit,
                "--cycles",
                self.opt_maxiter,
                "--acc",
                self.accuracy,
                "--namespace",
                name,
            )
        )

        if xtbinp:
            with open(f"{self.cwd}/{name}.inp", "wt") as inpf:
                inpf.write(xtbinp)

            cmd.extend(("--input", f"{name}.inp"))

        self(*cmd)

        with open(f"{self.cwd}/{name}.xtbopt.xyz") as f:
            nxyz = f.read()
            if not in_place:
                mol1 = deepcopy(mol)
                mol1.update_geom_from_xyz(nxyz, assert_single=True)
                return mol1
            else:
                mol.update_geom_from_xyz(nxyz, assert_single=True)

    def fix_long_bonds(
        self,
        mol: Molecule,
        rss_length_thresh: float = 3.0,
        rss_steps: float = 20,
        force_const: float = 0.05,
        target_len: float = 1.75,
        in_place: bool = False,
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

        for b in tbf:

            a1, a2 = mol.get_atom_idx(b.a1), mol.get_atom_idx(b.a2)
            inp += f"  distance: {a1+1}, {a2+1}, {tbf[b]:0.6f}\n"

        inp += "$scan\n  mode=concerted\n"
        for i, b in enumerate(tbf):
            inp += f"  {i+1}: {tbf[b]}, {target_len}, {rss_steps}\n"

        inp += "$end\n"

        m1 = self.optimize(mol, crit="sloppy", in_place=False, xtbinp=inp, fn_suffix=0)

        return m1


class CRESTDriver(ExternalDriver):
    """
    This driver provides functionality for
    --- CREST package (Grimme group) ---
    https://xtb-docs.readthedocs.io/en/latest/setup.html
    """

    WORKFLOWS = {
        "quick:v3:gfn2//gfnff": [],
    }
    PREFIX = "molli.crest"
    JOB_ID = 0

    def find_conformers(
        self,
        mol: Molecule,
        workflow: str = "gfn2-xtb",
        ewin: float = 15.0,  # in kJ/mol
        keep_rotamers: bool = True,
    ):
        """
        Take the molecule's geometry and find conformers
        ewin: energy window of 15 kJ/mol
        """
        raise NotImplementedError
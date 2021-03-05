"""
    This module provides interactions with XTB package
"""
import os
from ._core import ExternalDriver, DriverError
from ..dtypes import Molecule
from copy import deepcopy
from datetime import datetime
from glob import glob


class XTBDriver(ExternalDriver):
    """
        This driver provides functionality for
        --- XTB package (Grimme group) ---
        https://xtb-docs.readthedocs.io/en/latest/setup.html
    """
    METHODS = {"gfn2-xtb": "--gfn2", "gfn1-xtb": "--gfn1", "gfn-ff": "--gfnff"}
    PREFIX = "molli.xtb"
    JOB_ID = 0

    def minimize(self,
                 mol: Molecule,
                 method: str = 'gfn2-xtb',
                 crit: str = "tight",
                 maxiter: int = 200,
                 in_place: bool = False,
                 constraints: str = ""):
        """
            Perform energy minimization
            WARNING: BONDING TABLE IS NOT GOING TO BE UPDATED!
            Warning: make sure that the bond lengths are physically sound.
                        alternatively, use the fix_long_bonds routine
            GFN family of methods does not fix unphysical geometries (without constraints).

            if in_place == True: modify the input molecule 
        """
        assert method in self.METHODS

        jobid = self.__class__.JOB_ID

        xyz = mol.to_xyz()
        name = f'{self.PREFIX}.{os.getpid()}.{jobid:0>3}.minimize'  # pid for multiprocessing safety
        with open(f"{self.cwd}/{name}.g0.xyz", "wt") as f:
            f.write(xyz)

        if constraints:
            with open(f"{self.cwd}/{name}.c.inp", "wt") as f:
                f.write(constraints)

            self(
                "xtb",
                f"{name}.g0.xyz",
                self.METHODS[method],
                "--opt",
                crit,
                "--input",
                f"{name}.c.inp",
                "--namespace",
                name,
            )
        else:
            self("xtb", f"{name}.g0.xyz", self.METHODS[method], "--opt", crit,
                 "--cycles", maxiter, "--namespace", name)

        if in_place:
            mol1 = mol
        else:
            mol1 = deepcopy(mol)

        with open(f"{self.cwd}/{name}.xtbopt.xyz") as f:
            nxyz = f.read()
            if not in_place:
                mol1.update_geom_from_xyz(nxyz, assert_single=True)
            else:
                mol.update_geom_from_xyz(nxyz, assert_single=True)

        self.cleanup(regex=f"{name}.*")

        if not in_place:
            return mol1

    def gen_constraints(self,
                        mol: Molecule,
                        fc: float = 0.5,
                        maxcycle: int = 10,
                        distances: list = [],
                        angles: list = [],
                        dihedrals: list = [],
                        scans: list = [],
                        concerted=True):
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
            i1, i2, i3 = mol.get_atom_idx(a1), mol.get_atom_idx(
                a2), mol.get_atom_idx(a3)
            result += f"  angle: {i1+1}, {i2+1}, {i3+1}, {val}\n"

        for a1, a2, a3, a4, val in dihedrals:
            i1, i2, i3, i4 = mol.get_atom_idx(a1), mol.get_atom_idx(
                a2), mol.get_atom_idx(a3), mol.get_atom_idx(a4)
            result += f"  dihedral: {i1+1}, {i2+1}, {i3+1}, {i4+1}, {val}\n"

        if scans:
            result += "$scan\n"
            if concerted and len(scans) > 1:
                result += "  mode = concerted\n"
            for idx, start, end, steps in scans:
                result += f"  {idx}: {start}, {end}, {steps}\n"

        result += f"$opt\n  maxcycle={maxcycle}\n$end\n"

        return result

    def fix_long_bonds(self,
                       mol: Molecule,
                       method: str = 'gfn-ff',
                       lthresh: float = 5.0,
                       target: float = 1.5,
                       fc: float = 0.5,
                       nsteps: int = 30,
                       in_place: bool = True):
        """
            Detects long bonds by lower length threshold (artifacts of molecule joining) and
            performs a concerted relaxed surface scan to bring those fragments together.
            Geometry
        """

        dists = []
        scans = []

        for b in mol.bonds:
            l = mol.get_bond_length(b)
            if l >= lthresh:
                dists.append((b.a1, b.a2, l))
                scans.append((len(dists), l, target, nsteps))

        constr = self.gen_constraints(mol, distances=dists, scans=scans)

        m1 = self.minimize(mol,
                           method=method,
                           crit="sloppy",
                           in_place=in_place,
                           constraints=constr)

        if m1:
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
            workflow: str = 'gfn2-xtb',
            ewin: float = 15.0,  # in kJ/mol
            keep_rotamers: bool = True):
        """
            Take the molecule's geometry and find conformers
            ewin: energy window of 15 kJ/mol
        """
        raise NotImplementedError

from ._core import AsyncExternalDriver
from ..dtypes import Atom, Bond, Molecule, CartesianGeometry
from copy import deepcopy
from datetime import datetime
from glob import glob
from warnings import warn
from typing import List, Callable
from math import ceil, pi
import numpy as np
from itertools import combinations
import asyncio as aio
import re

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
        optimized = await self.xyz_optimize(g0_xyz, method = method, crit = crit, xtbinp = xtbinp,
                                                maxiter = 50,xyz_name = nn + '_g0.xyz')
        
        if not in_place:
            mol1 = deepcopy(mol)
            mol1.update_geom_from_xyz(optimized, assert_single=True)
            return mol1
        else:
            mol.update_geom_from_xyz(optimized, assert_single=True)

    async def xyz_optimize(
        self,
        xyz: str,
        method: str = 'gff',
        crit: str = "normal",
        xtbinp: str='',
        maxiter: int = 50,
        xyz_name: str = 'mol'
    ):
        # command that will be used to execute xtb package
        _cmd = f"""xtb {xyz_name}.xyz --{method} --opt {crit} --cycles {maxiter} {"--input param.inp" if xtbinp else ""} -P {self.nprocs}"""

        # pylint: disable=unused-variable
        code, files, stdout, stderr = await self.aexec(
            _cmd,
            inp_files={f"{xyz_name}.xyz": xyz, "param.inp": xtbinp},
            out_files=["xtbopt.xyz"],
        )

        if "xtbopt.xyz" in files:
            nxyz = files["xtbopt.xyz"]
            return nxyz
        else:
            raise FileNotFoundError("Could not locate xtb output file.")

    async def optimize_conformers(
        self,
        mol: Molecule,
        method: str = 'gff',
        crit: str = "normal",
        xtbinp: str='',
        maxiter: int = 50,
        in_place: bool = False,
    ):
        xyzs = mol.confs_to_xyzs()
        nn = mol.name

        optimized_confs = []
        for i, xyz in enumerate(xyzs):
            mol_name = nn + f'_{i}'
            optimized = await self.xyz_optimize(xyz, method=method, crit=crit, xtbinp=xtbinp,
                                                    maxiter=maxiter, xyz_name=mol_name)
            optimized_confs.append(optimized)

        geoms = [CartesianGeometry.from_xyz(conf)[0][0] for conf in optimized_confs]
        if in_place:
            mol.embed_conformers(*geoms, mode='w')
        else:
            mol1 = deepcopy(mol)
            mol1.embed_conformers(*geoms, mode='w')
            return mol1



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

    async def xyz_energy(self, xyz: str, method: str='gfn2', accuracy: float = 1.0):
        _cmd = f"""xtb struct.xyz --{method} --acc {accuracy:0.2f}"""

        code, files, stdout, stderr = await self.aexec(_cmd, inp_files={f"struct.xyz": xyz})
        
        # This is what we are trying to find in the output file
        # | TOTAL ENERGY             -172.541095318001 Eh   |

        for l in stdout.split('\n')[::-1]:
            if m := re.match(r"\s+\|\s+TOTAL ENERGY\s+(?P<eh>[0-9.-]+)\s+Eh\s+\|.*", l):
                return float(m['eh'])

    async def conformer_energies(self, mol: Molecule, method: str='gfn2', accuracy:float=1.0):
        """
        Returns relative conformer energies in kJ/mol
        The relative energies are referenced to the first conformer
        """
        xyzs = mol.confs_to_xyzs()
        nn = mol.name
        energies = []

        for i, xyz in enumerate(xyzs):
            conf_energy = await self.xyz_energy(xyz, method=method, accuracy=accuracy)
            energies.append(conf_energy)
        
        ref_energy = energies[0]
          
        return (np.array(energies) - ref_energy)*2625.5 # conversion to kJ/mol
    
    async def charges(self, mol: Molecule, method: str="gfn2", accuracy: float=0.5, net_charge: int = 0):
        """
            Compute atomic charges using XTB methodology

            FIXED 2022-02-23: added at least a partial support for total charge of the molecule.
            DO NOT USE unless you know what you are doing. -SAS
        """
        xyz = mol.to_xyz(n=0)
        _cmd = f"""xtb struct.xyz --sp --{method} --acc {accuracy:0.2f} --chrg {net_charge}"""

        code, files, stdout, stderr = await self.aexec(_cmd, inp_files={f"struct.xyz": xyz}, out_files=["charges"])
        
        charges = np.array(list(map(float, files["charges"].split())), dtype=np.float32)

        return charges 


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

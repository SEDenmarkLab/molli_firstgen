from ._core import AsyncExternalDriver
from ..dtypes import Atom, Bond, Molecule, CartesianGeometry
from ..dtypes.molecule import Orca_Out_Recognize
from copy import deepcopy
from datetime import datetime
from glob import glob
from warnings import warn
from typing import List, Callable
from itertools import combinations
import asyncio as aio
import numpy as np
import re

class AsyncORCADriver(AsyncExternalDriver):
    def __init__(self, name="", scratch_dir="", nprocs=1, encoding="utf8"):
        super().__init__(
            name=name, scratch_dir=scratch_dir, nprocs=nprocs, encoding=encoding)
    async def orca_basic_calc(
        self, 
        mol: Molecule,
        orca_path: str = '/opt/share/orca/5.0.2/orca', 
        ram_setting: str = '900',
        kohn_sham_type = 'rks',
        method: str = 'b3lyp',
        basis_set: str = 'def2_tzvp',
        calc_type = 'sp',
        addtl_settings = 'rijcosx def2/j tightscf nopop miniprint',
        charge: int = 0,
        spin_multiplicity: int = 1,
        already_completed: bool = False
        ):
        """
            General Orca Driver to Create a File and run calculations. Currently usable for the following calculations: "sp","opt","freq", "opt freq".
            This currently cannot recognize different calculation types in a backup directory since the files are built from the "Molecule Object" name.
            Consider doing different calculations in different folders to prevent loading incorrect files.
        """
        
        #Corrects xyz file to be usable in Orca
        full_xyz = mol.to_xyz()
        split_xyz_list = full_xyz.split('\n')
        split_xyz_list_fixed = [f'{x}\n' if i <= len(split_xyz_list)-3 else f'{x}' for i,x in enumerate(split_xyz_list)][2:]
        dft_xyz = ''.join(split_xyz_list_fixed)


        nn = mol.name
        
        _inp = f'''#{str.upper(calc_type)} {mol.name}

%maxcore {ram_setting}

%pal nprocs {self.nprocs} end

!{kohn_sham_type} {method} {basis_set} {calc_type} {addtl_settings}

*xyz {charge} {spin_multiplicity}
{dft_xyz}
*


'''
        #Removes spaces in opt_freq
        if calc_type == 'opt freq':
            calc_type = 'opt_freq'

        if already_completed:
            _cmd = f':'
        else:
            _cmd = f"""{orca_path} {nn}_{calc_type}.inp"""
        
        # print(_cmd)
        code, files, stdout, stderr = await self.aexec(
            _cmd,
            inp_files={f"{nn}_{calc_type}.inp": _inp},
            out_files = [f'{nn}_{calc_type}.hess',f'{nn}_{calc_type}.gbw'],
        )

        try:
            _hess = files[f'{nn}_{calc_type}.hess']
        except:
            _hess = None

        try:
            _gbw = files[f'{nn}_{calc_type}.gbw']
        except:
            _gbw = None

        try:
            _out = stdout
        except:
            # print(FileNotFoundError(f'{nn}_{calc_type}.out'))
            _out = None

        orca_obj = Orca_Out_Recognize(name = f'{nn}', output_file = _out, calc_type = calc_type, hess_file = _hess, gbw_file = _gbw)
        
        if orca_obj.orca_failed:
            print(f'{orca_obj.name}')
            print("Orca failed to converge. Here is the error given")
            print(orca_obj.end_lines)

        return orca_obj
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







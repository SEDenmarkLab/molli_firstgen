from ._core import AsyncExternalDriver
from ..dtypes import Atom, Bond, Molecule, CartesianGeometry
from glob import glob
import numpy as np
from itertools import combinations
# Ian checked/added
import pandas as pd
import io
import scipy.spatial.distance as dist
from asyncio.subprocess import PIPE
from typing import Callable, Dict, List, Awaitable
from tempfile import TemporaryDirectory as TempDir, NamedTemporaryFile as TempFile
import asyncio as aio
import os

def default_nw_inp():

    return """
title '{name}'
scratch_dir /scratch/nir2/nwchem/
memory 3000 mb
echo
start

geometry units angstroms noautoz
symmetry C1
{trimmed_xyz_block}
end

charge {charge}

driver 
    loose
    maxiter {geomiter}
end

basis
* library def2-svp 
end

dft
    xc b3lyp
    maxiter 100
end

{opt_or_not}

esp
    recalculate
    range 0.2
    probe 0.1
    spacing 0.025
end

task esp
"""

class AsyncNWCHEMDriver(AsyncExternalDriver):
    def __init__(self,
                name="",
                scratch_dir="",
                nprocs=1,
                # memory_total="800 mb",
                # memory_global="400 mb",
                encoding="utf8"):

        ''' 
        This is a prototype which relies on nwcrun alias for submission of nwchem jobs
        The
        '''
        super().__init__(
                name=name,
                scratch_dir=scratch_dir,
                nprocs=nprocs,
                encoding=encoding
        )
        # self.memory_total = memory_total
        # self.memory_global = memory_global
        with open(self.name+'_out.csv','a') as f:
            f.write('mol,espmin,espmax,min_a_chg,max_a_chg\n')


    async def optimize_esp(
        self,
        mol: Molecule,
        basis: str = "6-31g*",
        functional: str = "B3LYP",
        maxiter: int = 100,
        update_geom: bool = False,
        in_place: bool = False):

        ''' 
        This will update geometry for the structures only if update_geom is True.
        This is not implemented yet [7/21/2021]

        Default output is a tuple of espmin, espmax, minimum atomic charge, max atomic charge
        These are written to a file as they are generated so that names are correct even when files fail
        The .esp output (xyz + atomic charges format) is output at the end
        returns mol.name,espmin,espmax,min_a_chg,max_a_chg,files[f"{nn}.esp"]
        ''' 
        ## Need to get input file from mol with the necessary nwchem information in it
        self.maxiter = maxiter
        trim_xyz = mol.to_xyz()

        xyz_lines = '\n'.join(trim_xyz.splitlines(keepends=False)[2:])
        nw_lines = self._nwchem_composition(xyz_=xyz_lines,opt_=True)
        
        nn = mol.name #name of mol for input line fstring

        # Command will use mpirun to run nwchem parallelized
        # OR NOT if no parallelization requested
        if self.nprocs > 1:        
            _cmd = f"""mpiexec --prefix /opt/openmpi/openmpi314 -n {self.nprocs} nwchem {nn}.nw"""
        else:
            _cmd = f"""nwchem {nn}.nw """                

        # Later: stdout is log file, parse for energies, frequencies, or whatever
        # files dict has esp and grid files. Use to parse them for specific esp values of interest and 
        # return those as the final output 

        code,files,stdout,stderr = await self.aexec(
            _cmd,
            inp_files={f"{nn}.nw": nw_lines},
            out_files=[f"{nn}.esp",f"{nn}.grid"]
        )

        espmin,espmax,min_a_chg,max_a_chg = self.get_min_max_esp(files[f"{nn}.grid"],files[f"{nn}.esp"])
        with open(self.name+'_out.csv','a') as f:
            f.write(f"{mol.name},{espmin},{espmax},{min_a_chg},{max_a_chg}\n")
        return mol.name,espmin,espmax,min_a_chg,max_a_chg,files[f"{nn}.esp"]


        # see molli.drivers.xtb AsyncXTBDriver.fix_long_bonds code for inspiration

        # This next part uses xyz. The NWChem output is pretty easily converted to that. 
        # I could use this existing machinery to handle geometry updating. Might as well
        # make an xyz AND an nwchem input file, and use the xyz to update and save the mol
        # geometry at the end. 

        # This is for later: make it possible to return updated xyz geometries for a list of mol objects

    async def aexec(
        self, cmd: str, inp_files: Dict[str, str] = dict(), out_files: List[str] = []
    ):
        """
        Coroutine that asynchronously schedules a shell command to be executed
        Before the command is executed it writes temporary files (`inp_files`)
        After the command is executed, output files are harvested (`out_files`)
        returns: code, files, stdout, stderr
        """

        with TempDir(prefix=self.prefix, dir=self.scratch_dir) as td:
            self.writefiles(td, inp_files)
            # Asynchronous process spawning code goes here

            _cmd = f"cd {td}; {cmd}"

            proc = await aio.create_subprocess_shell(_cmd, stdout=PIPE, stderr=PIPE)
            # code = await proc.wait()
            while True:
                await aio.sleep(10)
                _out = await proc.stdout.read()
                stdout = _out.decode(self.encoding)

                if "Total times" in stdout[-1]:
                    await aio.sleep(60)
                    proc.terminate()
                    code = 0
                    break

                else:
                    if not (proc.returncode is None):
                        code = proc.returncode
                        break

            _out = await proc.stdout.read()
            _err = await proc.stderr.read()

            stdout = _out.decode(self.encoding)
            stderr = _err.decode(self.encoding)

            files = self.getfiles(td, out_files)

            if code:
                td_base = os.path.basename(td)
                with open(f"{self.scratch_dir}/{td_base}_dump_stdout.log", "wt") as f:
                    f.write(stdout)
                with open(f"{self.scratch_dir}/{td_base}_dump_stderr.log", "wt") as f:
                    f.write(stderr)
                    for fout in inp_files:
                        f.write(f"\n\n === {fout} ===\n\n")
                        f.write(inp_files[fout])

        return code, files, stdout, stderr

    def _nwchem_composition(self,xyz_,opt_=True,):
        '''
        This helper function takes xyz file string to nwchem file string with optional
        optimization or direct esp calculation. This is specifically designed for function
        optimize_esp for getting ESP min and max values from ammonium fragments.
        '''
        if opt_ == True:
            nw_string = default_nw_inp().format(name=self.name,
                                charge='1',
                                geomiter=str(self.maxiter),
                                opt_or_not = r'task dft optimize',
                                trimmed_xyz_block=xyz_,
                                # opt_criteria='loose'
                                )
            return nw_string
        elif opt_ == False:
            nw_string = default_nw_inp().format(name=self.name,
                                charge='1',
                                geomiter=str(self.maxiter),
                                opt_or_not = r'task dft energy',
                                trimmed_xyz_block=xyz_,
                                # opt_criteria='loose'
                                )
            return nw_string

    def _energy_parser(*args, **kwargs):
        '''
        This will parse the stdout (log file) to obtain useful structural information
        '''
        raise NotImplementedError("In development")
        

    def get_min_max_esp(self,grid_file,esp_file):
        '''
        This returns a tuple (ESPmin,ESPmax) in kcal/mol
        '''
        h_t_kcal = 627.509469
        mxyz_df,mc_df,mxyz_arr,mc_arr = self._process_nwespfile(esp_file,typ_='esp')
        gxyz_df,gc_df,gxyz_arr,gc_arr = self._process_nwespfile(grid_file,typ_='grid')
        esp_array = gc_arr*h_t_kcal
        espmin = np.min(esp_array)
        espmax = np.max(esp_array)
        max_a_chg = np.max(mc_arr)
        min_a_chg = np.min(mc_arr)
        return espmin,espmax,min_a_chg,max_a_chg


    def _measure_distance(self,grid_coord,atom_coord):
        '''
        This takes grid [x,y,z] floats and atom [x,y,z] floats and returns a scalar distance.
        Make sure to pass a numpy array row for each.
        '''
        p1 = grid_coord
        p2 = atom_coord
        stack = np.stack((p1,p2))
        distance_ = dist.pdist(stack)
        return distance_


    def _process_nwespfile(self,nw_xyz_file,typ_:str):
        '''
        This accepts a file string for .esp or .grid, and returns tuple xyz_df,c_df,xyz_arr,c_arr for various calculation tasks
        When parsing esp file, the _df outputs retain the atom label as the pandas index list
        '''
        if typ_ == '': 
            raise ValueError("Must pass typ_ argument as string esp or grid")

        elif typ_ == 'esp':
            trimmed_esp_block = '\n'.join(nw_xyz_file.splitlines(keepends=False)[2:]) #Two lines thrown out
            # Read in as dataframe; this will pull x, y, z, chg column data in
            esp_df = pd.read_csv(io.StringIO(trimmed_esp_block),delim_whitespace=True,header=None,names=['x','y','z','chg'])
            xyz = esp_df[['x','y','z']] #just coordinates for computing distances
            c_df = esp_df['chg'] #charge for corresponding xyz data point
            #esp coords
            xyz_arr = xyz.to_numpy()
            #charge coords
            c_arr = c_df.to_numpy()
            return xyz,c_df,xyz_arr,c_arr

        elif typ_ == 'grid':
            trimmed_grid_block = '\n'.join(nw_xyz_file.splitlines(keepends=False)[1:]) #One line thrown out
            # Read in as dataframe; this will pull x, y, z, chg column data in
            esp_df = pd.read_csv(io.StringIO(trimmed_grid_block),delim_whitespace=True,header=None,names=['x','y','z','chg'])
            xyz = esp_df[['x','y','z']] #just coordinates for computing distances
            c_df = esp_df['chg'] #charge for corresponding xyz data point
            #esp coords
            xyz_arr = xyz.to_numpy()
            #charge coords
            c_arr = c_df.to_numpy()
            return xyz,c_df,xyz_arr,c_arr

        else: 
            raise ValueError("You can only pass esp or grid as typ_ argument for this function")







# Old idea for get_min_max_esp built-in parse. Depreciated. 
# with open(grid_file,'r') as f:
#     grid_chgs = []
#     for line in f.readlines():
#         x, y, z, chrg = line.strip().split()
#         grid_chgs.append(float(chrg))
# Remove first two lines of standard nwchem output block
# with open(esp_file,'r') as f:
#     esplines = [g.strip() for g in f.readlines()[2:]]
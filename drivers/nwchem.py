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
# Ian checked/added
import pandas as pd
import io
import scipy.spatial.distance as dist

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
    maxiter {geomiter}
end

basis
* library 6-31g*
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
    def __init__(self,name="",
                scratch_dir="",
                nprocs=1,
                memory_total="800 mb",
                memory_global="400 mb",
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
        self.memory_total = memory_total
        self.memory_global = memory_global

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

        Default output is a tuple of esp min, esp max, and energy
        ''' 
        ## Need to get input file from mol with the necessary nwchem information in it
        
        trim_xyz = mol.to_xyz()

        xyz_lines = '\n'.join(trim_xyz.splitlines(keepends=False)[2:])
        nw_lines = self._nwchem_composition(xyz_=trim_xyz,opt_=True)
        
        nn = mol.name #name of mol for input line fstring

        # Command will use mpirun to run nwchem parallelized
        # OR NOT if no parallelization requested
        if self.nprocs > 1:        
            _cmd = f"""mpirun -np {self.nprocs} nwchem {nn}.nw"""
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


        # see molli.drivers.xtb AsyncXTBDriver.fix_long_bonds code for inspiration

        # This next part uses xyz. The NWChem output is pretty easily converted to that. 
        # I could use this existing machinery to handle geometry updating. Might as well
        # make an xyz AND an nwchem input file, and use the xyz to update and save the mol
        # geometry at the end. 

        # This is for later: make it possible to return updated xyz geometries for a list of mol objects

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
                                trimmed_xyz_block=xyz_
                                )
            return nw_string
        elif opt_ == False:
            nw_string = default_nw_inp().format(name=self.name,
                                charge='1',
                                geomiter=str(self.maxiter),
                                opt_or_not = r'',
                                trimmed_xyz_block=xyz_
                                )
            return nw_string

    def _energy_parser(*args, **kwargs):
        '''
        This will parse the stdout (log file) to obtain useful structural information
        '''
        raise NotImplementedError("In development")
        
    def _esp_parser(self,esp_file,grid_file):
        '''
        This will parse out esp descriptors.
        ESp first, then grid file
        '''
        raise NotImplementedError("FUCK OPENBABEL")

    def _convert_grid_file(self,grid_file_lines):
        ...
    def get_min_max_esp(self,grid_file,esp_file):
        '''
        This returns a tuple (ESPmin,ESPmax)
        '''


        mxyz_df,mc_df,mxyz_arr,mc_arr = self._process_nwespfile(esp_file,typ_='esp')
        gxyz_df,gc_df,gxyz_arr,gc_arr = self._process_nwespfile(grid_file,typ_='grid')



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
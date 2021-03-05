"""
    Contains all definitions for an open babel driver
"""
import os
from ..dtypes import Molecule, CartesianGeometry

from ._core import ExternalDriver


class OpenBabelDriver(ExternalDriver):
    """
        Perform some molecular operations with open babel
    """
    def __init__(self, cwd: str = './'):
        super().__init__(cwd=cwd)

    def convert(self,
                mol_text: str,
                *args,
                src: str = "pdb",
                dest: str = "mol2"):
        """ 
        Convert string molecule representation into something else.
        Useful argument: '-h' adds hydrogens
        """
        # self._prep_str_input(mol_text)
        return self("obabel", f"-i{src}", f"-o{dest}", *args, inp=mol_text)

    def add_hydrogens(self, mol_text, fmt: str = 'mol2'):
        """
            Add hydrogen atoms
        """
        return self.convert(mol_text, "-h", src=fmt, dest=fmt)

    def minimize(self,
                 mol_text: str,
                 src: str = "mol2",
                 dest: str = "mol2",
                 ff: str = "UFF",
                 n: int = 500,
                 c: float = 1.0e-4):
        """
            Perform energy minimization.
            Since the authors of openbabel didn't bother make the input strategy consistent
            We need to prepare an input file first
        """
        with open(os.path.join(self.cwd, f"temp01min.{src}"), "wt") as f:
            f.write(mol_text)

        pdb = self("obminimize", "-ff", ff, "-n", n, "-c", c,
                   f"temp01min.{src}")

        return self.convert(pdb, src="pdb", dest=dest)

    def hadd_opt(self,
                 mol: Molecule,
                 ff='UFF',
                 n=500,
                 c=1e-4,
                 update_source=True) -> Molecule:
        """ 
            Take a Molecule instance, add hydrogens and optimize.
            If conformers are present, they are ignored in this version.
            NOTE This will change in the future.

            WARNING: if mol contains hydrogen atoms, they will not be added per OBABEL behavior

            if update_source == True: modifies the original molecule [not implemented]
        """

        mol2 = mol.to_mol2()
        mol2_h = self.add_hydrogens(mol2, fmt="mol2")

        res = Molecule.from_mol2(mol2_h)

        xyz_h_opt = self.minimize(mol2_h,
                                  src='mol2',
                                  dest='xyz',
                                  ff=ff,
                                  n=n,
                                  c=c)
        res.update_geom_from_xyz(xyz_h_opt, assert_single=True)
        return res

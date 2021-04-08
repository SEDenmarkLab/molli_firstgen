from ._core import AsyncExternalDriver
from ..dtypes import Molecule


class AsyncOpenBabelDriver(AsyncExternalDriver):
    """
    Asynchronous version of the OpenBabel driver
    """

    async def convert(self, mol_text: str, *args, src: str = "pdb", dest: str = "mol2"):
        """
        Convert string molecule representation into something else.
        Useful argument: '-h' adds hydrogens
        """
        _cmd = f"obabel -i{src} source -o{dest} -O converted " + " ".join(args)

        # pylint: disable=unused-variable
        code, files, stdin, stdout = await self.aexec(
            _cmd, inp_files={"source": mol_text}, out_files=["converted"]
        )
        res = files["converted"]
        return res

    async def add_hydrogens(self, mol_text, fmt: str = "mol2"):
        """
        Add hydrogen atoms
        """
        return await self.convert(mol_text, "-h", src=fmt, dest=fmt)

    async def minimize(
        self,
        mol_text: str,
        src: str = "mol2",
        dest: str = "mol2",
        ff: str = "UFF",
        n: int = 500,
        c: float = 1.0e-4,
    ):
        """
        Perform energy minimization of the input structure
        """
        _cmd = f"obminimize -ff {ff} -n {n} -c {c} input.{src}"

        # pylint: disable=unused-variable
        code, files, stdout, stderr = await self.aexec(
            _cmd, inp_files={f"input.{src}": mol_text}
        )

        return await self.convert(stdout, src="pdb", dest=dest)

    async def hadd_opt(
        self, mol: Molecule, ff="UFF", n=500, c=1e-4, update_source=True
    ) -> Molecule:
        """
        Take a Molecule instance, add hydrogens and do preliminary optimization.
        If conformers are present, they are ignored in this version.
        NOTE This will change in the future.

        WARNING: if mol contains hydrogen atoms, they will not be added per OBABEL behavior

        if update_source == True: modifies the original molecule [not implemented]
        """

        mol2 = mol.to_mol2()
        mol2_h = await self.add_hydrogens(mol2, fmt="mol2")

        res = Molecule.from_mol2(mol2_h)

        xyz_h_opt = await self.minimize(mol2_h, src="mol2", dest="xyz", ff=ff, n=n, c=c)
        res.update_geom_from_xyz(xyz_h_opt, assert_single=True)

        return res
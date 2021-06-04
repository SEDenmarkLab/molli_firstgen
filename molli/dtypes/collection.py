from __future__ import annotations
from .molecule import Molecule
from datetime import datetime
from typing import List, Callable, Any, Awaitable
import json
from zipfile import ZipFile
from itertools import product, combinations_with_replacement
from multiprocessing import Pool
import asyncio as aio
import numpy as np


class Collection:
    """
    This class provides convenience when handling molecule collections (zip files)
    Performance of this class is limited by the time to load and parse massive objects
    Therefore it is refactored a little bit
    """

    def __init__(self, name: str = "", molecules: List[Molecule] = []):
        self.name = name
        self.molecules = molecules
        self.mol_index = [x.name for x in molecules]

    def add(self, m: Molecule):
        self.molecules.append(m)
        self.mol_index.append(m.name)

    def extend(self, c: Collection):
        for m in c:
            self.add(m)

    def __len__(self):
        return len(self.molecules)

    def __iter__(self):
        return iter(self.molecules)

    def __getitem__(self, item) -> Molecule:
        if isinstance(item, int):
            return self.molecules[item]
        elif isinstance(item, str):
            return self.molecules[self.mol_index.index(item)]
        else:
            return Collection(self.name, self.molecules[item])

    def bounding_box(self):
        """
        Get the rectangular space that encompasses all atoms in all conformers
        """
        mins = []
        maxs = []

        for m in self.molecules:
            rmin, rmax = m.bounding_box()
            mins.append(rmin)
            maxs.append(rmax)

        rmin = np.min(mins, axis=(0,))
        rmax = np.max(maxs, axis=(0,))

        return rmin, rmax

    def to_multixyz(self, fn: str = None):
        """
        Return a multixyz representation of molecules
        """
        result = ""
        for m in self:
            result += m.to_xyz()

        if fn:
            with open(fn, "wt") as f:
                f.write(result)

        return result

    def applyawt(self, aw: Awaitable, timeout=None, show_progress=True, update=1000):
        """
        This function is designed to mimic applyfx, but be useful with awaitables.
        Right now does not support failure handling, this is something
        """

        def inner(*args, **kwargs):
            async def f():
                results = []
                L = len(self.mol_index)
                start = datetime.now()
                if show_progress:
                    print(f"\nApplying [{aw.__name__}] to {L} molecules:")

                for i, m in enumerate(self.molecules):
                    results.append(
                        await aio.wait_for(aw(m, *args, **kwargs), timeout=timeout)
                    )
                    if show_progress and not (i + 1) % update:
                        print(
                            f"{i+1:>10} molecules processed ({(i+1)/L:>6.2%}) Total WCT: {datetime.now() - start}",
                            flush=True,
                        )

                if show_progress:
                    print(f"Complete! Total WCT: {datetime.now() - start}\n")

                return results

            return aio.run(f())

        return inner

    def applyfx(self, fx: Callable[[Molecule], Any]):
        """
        An incredibly useful *decorator* that applies a given function to all molecules in the library.
        May take advantage of multiprocessing by specifying nprocs
        Provides visual confirmation

        if fx returns Molecule objects, they are then assembled in a collection.
        Otherwise, it returns a list.
        """

        def inner(workers=1, show_progress=True, update=1000):
            result = []

            L = len(self.mol_index)

            start = datetime.now()

            if show_progress:
                print(f"\nApplying [{fx.__name__}] to {L} molecules:")

            if workers == 1:
                for i, m in enumerate(self.molecules):
                    result.append(fx(m))
                    if show_progress and not (i + 1) % update:
                        print(
                            f"{i+1:>10} molecules processed ({(i+1)/L:>6.2%}) Total WCT: {datetime.now() - start}",
                            flush=True,
                        )

            if workers > 1:

                total = 0

                batch_size = max(update, workers)
                n_batches = L // batch_size + 1

                result = []

                for i in range(n_batches):
                    chunk = self.molecules[
                        i * batch_size : min((i + 1) * batch_size, L)
                    ]

                    with Pool(workers) as pool:
                        result.extend(pool.map(fx, chunk))

                    total += len(chunk)

                    if show_progress:
                        print(
                            f"{total:>10} molecules processed ({total/L:>6.2%}) Total WCT: {datetime.now() - start}",
                            flush=True,
                        )

            if show_progress:
                print(f"Complete! Total WCT: {datetime.now() - start}\n")

            if all(isinstance(x, Molecule) for x in result):
                return Collection(name=fx.__name__ + self.name, molecules=result)
            else:
                return result

        return inner

    def to_zip(self, fpath=""):
        """
        Create a molecule archive
        """
        meta = {"name": self.name, "idx": self.mol_index, "files": []}
        with ZipFile(fpath, mode="w", compression=0, allowZip64=True) as zf:
            for i, m in enumerate(self.molecules):
                zf.writestr(f"{i+1}.xml", m.to_xml())

                meta["files"].append(f"{i+1}.xml")

            zf.writestr("__molli__", json.dumps(meta))

    @classmethod
    def merge(cls, *collections: Collection, name="merged"):
        """
        Self-explanatory
        """
        res = cls(name=name)

        for c in collections:
            res.extend(c)

        return res

    @classmethod
    def from_zip(cls, fpath: str):

        # restricted = ['__molli__']
        molecules = []

        with ZipFile(fpath, mode="r", compression=0, allowZip64=True) as zf:
            with zf.open("__molli__") as f:
                meta = json.load(f)

            if "files" in meta:
                for fn in meta["files"]:
                    with zf.open(fn, "r") as f:
                        molecules.append(Molecule.from_file(f))

            else:
                for fn in zf.namelist():
                    if fn in ["__molli__"]:
                        pass
                    else:
                        with zf.open(fn, "r") as f:
                            molecules.append(Molecule.from_file(f))

        return cls(name=meta["name"], molecules=molecules)

    @classmethod
    def join(
        cls,
        mc1: Collection,
        mc2: Collection,
        ap1: str,
        ap2: str,
        dist: float = 10.0,
        nprocs=1,
    ):
        """
        A great function that joins fragments in collections (!)
        VERY USEFUL FOR IN SILICO LIBRARY GENERATION
        """
        molecules = []
        for m1, m2 in product(mc1, mc2):
            molecules.append(Molecule.join_ap(m1, m2, ap1=ap1, ap2=ap2, dist=dist))

        return cls(name=f"{mc1.name}_{mc2.name}", molecules=molecules)


class CollectionFile:
    """
    This context manager provides access to Molecule items from a collection file.
    It provides lower level interactions with the file, such as modifications to individual files.

    It is also more convenient to use when only a few items need to be processed (big libraries tend to be cumbersome to load)

    Use this when you need to faster access to individual files rather than the entire collection
    if save: upon exiting the molecule objects in the zip file will be updated

    """

    def __init__(self, fpath: str, save_on_exit: bool = False):
        self.fpath = fpath
        self._save_on_exit = save_on_exit
        self._to_be_updated = []

    def __enter__(self):
        self._fstream = ZipFile(self.fpath, "r")
        with self._fstream.open("__molli__") as _mf:
            self._meta = json.load(_mf)
        self._collection = Collection(self._meta["name"])
        self.name = self._collection.name
        return self

    def __getitem__(self, item: str):
        if hasattr(self, "_collection") and item in self._collection.mol_index:
            # if collection (a buffer for molecule objects) exists
            return self._collection[item]
        elif hasattr(self, "_fstream") and hasattr(self, "_meta"):
            # ie if the file is open
            # and the item was not located in the existing collection
            idx = self._meta["idx"].index(item)
            with self._fstream.open(f"{idx+1}.xml") as mf:
                m = Molecule.from_file(mf)
                self._collection.add(m)

            return m

        else:
            raise IOError(
                f"Unable to import a molecule {item}. Not found in buffer, and the file stream is closed."
            )

    def save(self):
        with ZipFile(self.fpath, "a") as zf:
            for m in self._collection.mol_index:
                if m not in self._meta["idx"]:
                    raise IndexError("Malformed zip")
                else:
                    fn = self._meta["files"][self._meta["idx"].index(m)]
                    with zf.open(fn, "w") as mf:
                        mf.write(bytes(self._collection[m].to_xml(), encoding="utf8"))

    def __exit__(self, *args):
        self._fstream.close()

        if self._save_on_exit:
            self.save()

        del self._fstream

from __future__ import annotations
from .molecule import Molecule
from datetime import datetime
from typing import List, Callable, Any
import json
from zipfile import ZipFile
from itertools import product, combinations_with_replacement
from multiprocessing import Pool


class MolCollection:
    """
        This class provides convenience when handling molecule collections (zip files)
        Performance of this class is limited by the time to load and parse massive objects
        Therefore it is refactored a little bit
    """
    def __init__(self, name: str = "", molecules: List[Molecule] = []):
        self.name = name
        self.molecules = molecules
        self.mol_index = [x.name for x in molecules]

    def __len__(self):
        return len(self.molecules)

    def __iter__(self):
        return iter(self.molecules)

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.molecules[item]
        elif isinstance(item, str):
            return self.molecules[self.mol_index.index(item)]

    def applyfx(self, fx: Callable[[Molecule], Any]):
        """
        An incredibly useful *decorator* that applies a given function to all molecules in the library.
        May take advantage of multiprocessing by specifying nprocs
        Provides visual confirmation
        """
        def inner(nprocs=1, show_progress=True, update=1000):
            result = []

            L = len(self.mol_index)

            start = datetime.now()

            if show_progress:
                print(f"\nApplying [{fx.__name__}] to {L} molecules:")

            if nprocs == 1:
                for i, m in enumerate(self.molecules):
                    result.append(fx(m))
                    if show_progress and not (i + 1) % update:
                        print(
                            f"{i+1:>10} molecules processed ({(i+1)/L:>6.2%}) Total WCT: {datetime.now() - start}"
                        )

            if nprocs > 1:
                # This is where the multiprocessing fun starts. Buckle up!
                if show_progress:
                    print(f"(in {nprocs} parallel processes)")
                workers = Pool(nprocs)

                total = 0

                batch_size = max(update, nprocs)
                n_batches = L // batch_size + 1

                result = []

                for i in range(n_batches):
                    chunk = self.molecules[i * batch_size :\
                         min((i + 1) * batch_size, L)]
                    result.extend(workers.map(fx, chunk))
                    total += len(chunk)

                    if show_progress:
                        print(
                            f"{total:>10} molecules processed ({total/L:>6.2%}) Total WCT: {datetime.now() - start}"
                        )

            if show_progress:
                print(f"Complete! Total WCT: {datetime.now() - start}\n")
            return result

        return inner

    def to_zip(self, fpath=''):
        """
            Create a molecule archive
        """
        meta = {"name": self.name, "idx": self.mol_index, "files": []}
        with ZipFile(fpath, mode="w", compression=0, allowZip64=True) as zf:
            for i, m in enumerate(self.molecules):
                zf.writestr(f"{i+1}.xml", m.to_xml())

                meta['files'].append(f"{i+1}.xml")

            zf.writestr("__molli__", json.dumps(meta))

    @classmethod
    def from_zip(cls, fpath: str):

        # restricted = ['__molli__']
        molecules = []

        with ZipFile(fpath, mode="r", compression=0, allowZip64=True) as zf:
            with zf.open("__molli__") as f:
                meta = json.load(f)

            if "files" in meta:
                for fn in meta["files"]:
                    with zf.open(fn, 'r') as f:
                        molecules.append(Molecule.from_file(f))

            else:
                for fn in zf.namelist():
                    if fn in ['__molli__']:
                        pass
                    else:
                        with zf.open(fn, 'r') as f:
                            molecules.append(Molecule.from_file(f))

        return cls(name=meta['name'], molecules=molecules)

    @classmethod
    def join(cls,
             mc1: MolCollection,
             mc2: MolCollection,
             ap1: str,
             ap2: str,
             dist: float = 10.0):
        """
            A great function that joins fragments in collections (!)
            VERY USEFUL FOR IN SILICO LIBRARY GENERATION 
        """

        molecules = []

        for m1, m2 in product(mc1, mc2):
            molecules.append(
                Molecule.join_ap(m1, m2, ap1=ap1, ap2=ap2, dist=dist))

        return cls(name=f"{mc1.name}_{mc2.name}", molecules=molecules)


class MolCollectionFile:
    """
        This context manager provides access to Molecule items from a collection 
        when loading all of them in the memory is not an optimal strategy

        Use this when you need to faster access to individual files rather than the entire collection
        if save: upon exiting the molecule objects in the zip file will be updated

    """
    def __init__(self, fpath: str, save_on_exit: bool = False):
        self.fpath = fpath
        self._save_on_exit = save_on_exit
        self._to_be_updated = []

    def __enter__(self):
        self.__collection: MolCollection = None
        self.__fstream = ZipFile(self.fpath, "r")
        self._meta = self.__fstream.read('__molli__')
        return self

    def __getitem__(self, item: str):
        if hasattr(self,
                   '__collection') and item in self.__collection.mol_index:
            # if collection (a buffer for molecule objects) exists
            return self.__collection[item]
        if hasattr(self, '__fstream') and hasattr(self, '_meta'):
            # ie if the file is open
            # and the item was not located in the existing collection
            pass
        else:
            raise IOError(
                f"Unable to import a molecule {item}. Not found in buffer, and the file stream is closed."
            )

    def __exit__(self, *args):
        self.__fstream.close()

        if self._save_on_exit:
            ...
            # TODO: code that updates the molecules if

        del self.__fstream

"""

"""

import molli as ml
from glob import glob



def main(argv: list):
    print("Compiling molecules")
    print(*argv)
    mfs = glob(argv[0])
    print(f"Found: {len(mfs)} files")

    _mols = []
    for f in mfs:
        _mols.append(ml.Molecule.from_xml(f)) 

    _mols.sort(key = lambda x: x.name) 

    mols = ml.Collection("compilation", _mols)

    mols.to_zip(argv[1])      


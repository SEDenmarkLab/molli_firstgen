"""

"""

import molli as ml
from glob import glob

def main(argv: list):

    print("=== Compiling molecules ===")

    print("Locating input files ...")

    if len(argv) > 2:
        print(f"\t... found {len(argv[1:])} input files. Using those.")
        mfs = argv
    else:
        mfs = glob(argv[1])
        print(f"\t... found one expression matching {len(mfs)} files.")

    print(f"Searching for unique molecule containers")

    molecules = dict()

    for f in mfs:
        try:
            m = ml.Molecule.from_xml(f)
        except:
            print(f"[!f] {f}", flush=True)
        else:
            if m.name not in molecules:
                molecules[m.name] = m
                print(f"[+m] {len(m.conformers):>5} {m.name}", flush=True)
            elif m.name in molecules and len(molecules[m.name].conformers) < len(m.conformers):
                molecules[m.name] = m
                print(f"[*m] {len(m.conformers):>5} {m.name}", flush=True)
    
    print(f"Found {len(molecules)} unique molecules.")
    
    _mols = [molecules[x] for x in molecules]

    _mols.sort(key = lambda x: x.name) 

    mols = ml.Collection("compilation", _mols)

    mols.to_zip(argv[0])      


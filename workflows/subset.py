import molli as ml
from sys import argv
from subprocess import run, DEVNULL
from tempfile import mktemp
import argparse as ap

# c = ml.Collection.from_zip(argv[1])

parser = ap.ArgumentParser("show")

parser.add_argument(
    "zipfile",
)

meg_show = parser.add_mutually_exclusive_group()

meg_show.add_argument(
    "-itx", "--input_txt",
    action="store",
    metavar="<list.txt>",
    help="compile a subset of library to show (using default mol representation). Reads from a text file where every line is a molecule name."    
)

parser.add_argument(
    "-o", "--output",
    action="store",
    metavar="<output.zip>",
    help="use the following command to visualize a temporary file",
)

def main(argv):

    parsed = parser.parse_args(argv)

    
    if parsed.input_txt:
        c = ml.CollectionFile(parsed.zipfile)

        _mols = []

        with c, open(parsed.input_txt) as inpf:
            allnames = [x.strip() for x in inpf.read().split("\n")]
            for n in allnames:
                _mols.append(c[n])

        out_collection = ml.Collection(name="subset", molecules=_mols)

        tf = mktemp(suffix=".zip", prefix="ml_subset", dir='.')
        out_collection.to_zip(tf)
        
        print(f"A subset of {len(out_collection)} molecules was written to <{tf}>")
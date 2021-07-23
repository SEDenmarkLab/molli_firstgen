import molli as ml
from sys import argv
from subprocess import run, DEVNULL
from tempfile import NamedTemporaryFile
import argparse as ap

# c = ml.Collection.from_zip(argv[1])

parser = ap.ArgumentParser("show")

parser.add_argument(
    "zipfile",
)

meg_show = parser.add_mutually_exclusive_group()

meg_show.add_argument(
    "-c", "--show_conformers",
    action="store",
    metavar="<mol_name>",
    help="Show conformers for a given molecule"
)

meg_show.add_argument(
    "-l", "--show_library",
    action="store_true",
    help="Shows the entire library of molecules using the default conformer"
)

meg_show.add_argument(
    "-s", "--show_subset",
    action="store",
    help="compile a subset of library to show (using default mol representation)"    
)

parser.add_argument(
    "-cmd", "--command",
    action="store",
    metavar="<vis_program>",
    default="avogadro",
    help="use the following command to visualize a temporary file",
)

def main(argv):

    parsed = parser.parse_args(argv)


    if parsed.show_library:    
        c = ml.Collection.from_zip(parsed.zipfile)
        c: ml.Collection
        with NamedTemporaryFile("w+t", dir=".", suffix=".xyz") as f:
            f.write(c.to_multixyz())
            run([parsed.command, f.name], stderr=DEVNULL, stdout=DEVNULL)
    
    elif parsed.show_conformers:
        c = ml.CollectionFile(parsed.zipfile)
        with c, NamedTemporaryFile("w+t", dir=".", suffix=".mol2") as tf:
            m = c[parsed.show_conformers]

            allconfs = m.confs_to_molecules()
            for x in allconfs:
                tf.write(x.to_mol2())
                tf.write("\n")
            
            run([parsed.command, tf.name], stderr=DEVNULL, stdout=DEVNULL)
    
    elif parsed.show_subset:
        c = ml.CollectionFile(parsed.zipfile)

        with c, NamedTemporaryFile("w+t", dir=".", suffix=".mol2") as tf:

            with open(parsed.show_subset) as inpf:
                for l in inpf:
                    nn = l.strip()            
                    tf.write(c[nn].to_mol2())
                    tf.write("\n")
            
            run([parsed.command, tf.name], stderr=DEVNULL, stdout=DEVNULL)
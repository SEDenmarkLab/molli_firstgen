import molli as ml
from sys import argv
from subprocess import run, DEVNULL
from tempfile import NamedTemporaryFile
import argparse as ap

parser = ap.ArgumentParser("show")

parser.add_argument("zipfile")
parser.add_argument("molecule")

parser.add_argument(
    "-cmd", "--command",
    action="store",
    metavar="<vis_program>",
    default="avogadro",
    help="use the following command to visualize a temporary file",
)

def main(argv):

    parsed = parser.parse_args(argv)

    c = ml.CollectionFile(parsed.zipfile)
    
    with NamedTemporaryFile("w+t", suffix=".mol2") as f, c:
        f.write(c[parsed.molecule].to_mol2())
        run([parsed.command, f.name], stderr=DEVNULL, stdout=DEVNULL)

        if input(f"Do you want to update the molecule? [y/n]") == "y":
            print("Will eventually update, stub for now")

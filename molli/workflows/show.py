import molli as ml
from sys import argv
from subprocess import run, DEVNULL
from tempfile import NamedTemporaryFile
# c = ml.Collection.from_zip(argv[1])

def main(argv):

    c = ml.Collection.from_zip(argv[0])

    c: ml.Collection

    with NamedTemporaryFile("w+t", dir=".", suffix=".xyz") as f:
        f.write(c.to_multixyz())

        run(['jmol', f.name], stderr=DEVNULL, stdout=DEVNULL)
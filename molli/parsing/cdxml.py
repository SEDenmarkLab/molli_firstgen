from ..dtypes import Molecule, Atom, Bond, CartesianGeometry, Collection
from typing import List
from xml.etree import cElementTree as et
import os
import numpy as np
import os
import re


def parse_pos(p: str):
    return list(map(float, p.split()))


def dist(a, b):
    return np.sqrt(np.sum((np.array(a) - np.array(b)) ** 2))


def split_cdxml(file_path: str) -> Collection:
    """
    Split a single cdxml file into a collection of molecules
    """

    name = os.path.basename(file_path).rsplit(".")[0]

    with open(file_path) as f:
        xml = et.parse(f)

    rt = xml.getroot()

    # =======================================
    zu = float(rt.attrib["BondLength"])

    # Find all molecular fragments and text boxes
    fragments = rt.findall(".//fragment")

    # So apparently chemdraw occasionally places fragments inside fragments...
    # This makes sure that those don't get counted
    for x in rt.findall(".//fragment//fragment"):
        fragments.remove(x)

    textboxes = rt.findall("./*/t") + rt.findall("./*/group/t")

    assert len(fragments) == len(
        textboxes
    ), f"Received {len(fragments)} fragments and {len(textboxes)} tboxes"

    # Extract labels and coordinates of said labels
    # Coordinates are required to assign labels to correct coordinates

    labels = []
    label_coord = []

    for tb in textboxes:
        l, t, r, b = parse_pos(tb.attrib["BoundingBox"])
        label_coord.append([(r + l) / 2, (b + t) / 2])
        labels.append(tb.find("./s").text)

    # Iterate over fragments
    # Convert 2d geometry to Molecule files

    molecules = []

    for frag in fragments:

        l, t, r, b = parse_pos(frag.attrib["BoundingBox"])
        frag_centroid = [(r + l) / 2, (b + t) / 2]

        atoms = []
        atom_ids = []
        bonds = []
        # bond_ids = []
        geom = []

        # Iterate over nodes
        for i, n in enumerate(frag.findall("./n")):  # pylint: disable=unused-variable
            x, y = parse_pos(n.attrib["p"])
            atom_id = n.attrib["id"]
            a = n.find("./t/s")

            if a == None:
                atom = Atom("C", "C", "C")
            elif a != None and a.text[0] == "#":
                atom = Atom("Cl", a.text[1:], "Cl", ap=True)
            else:
                atom = Atom(a.text, a.text, a.text)

            atoms.append(atom)
            atom_ids.append(atom_id)

            geom.append([x, y, 0.0])

        # Iterate over bonds
        for b in frag.findall("./b"):
            id1, id2 = b.attrib["B"], b.attrib["E"]
            bt = "1" if not "Order" in b.attrib else b.attrib["Order"]

            # ==========================================================
            #      If the bond contains stereochemical indication:
            #           Add z-coordinate hints
            # ==========================================================

            if "Display" in b.attrib:
                f1, f2 = atom_ids.index(id1), atom_ids.index(id2)
                ZSC = -1

                if b.attrib["Display"] == "WedgeBegin":
                    geom[f2][2] = ZSC * zu + geom[f1][2]

                if b.attrib["Display"] == "WedgedHashBegin":
                    geom[f2][2] = -ZSC * zu + geom[f1][2]

                if b.attrib["Display"] == "WedgeEnd":
                    geom[f1][2] = ZSC * zu + geom[f2][2]

                if b.attrib["Display"] == "WedgedHashEnd":
                    geom[f1][2] = -ZSC * zu + geom[f2][2]

                if b.attrib["Display"] == "Bold":
                    geom[f1][2] = ZSC * zu
                    geom[f2][2] = ZSC * zu

                if b.attrib["Display"] == "Hashed":
                    geom[f1][2] = -ZSC * zu
                    geom[f2][2] = -ZSC * zu

            bond = Bond(atoms[atom_ids.index(id1)], atoms[atom_ids.index(id2)], bt)
            bonds.append(bond)

        closest = np.argmin([dist(frag_centroid, x) for x in label_coord])
        mol_name = labels[closest]

        mol = Molecule(mol_name, atoms=atoms, bonds=bonds, geom=CartesianGeometry(geom))
        mol.fix_geom()

        molecules.append(mol)

    return Collection(name=name, molecules=molecules)
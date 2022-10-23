from __future__ import annotations
from .geometry import CartesianGeometry, rotation_matrix
from typing import List, Dict, Union, Set, Any, Callable, overload
import os
import json

# from ..parsing.mol2 import get_mol2_block_lines
from copy import deepcopy
import numpy as np
from ..ftypes import parse_xyz
from xml.dom import minidom as xmd
from xml.etree.cElementTree import parse as xparse
from io import IOBase


def yield_mol2_block_lines(title, text):
    """
    Yield lines delimited by @TRIPOS<some_name> [...]
    """

    lines = [x.strip() for x in str(text).splitlines()]

    s = lines.index("@<TRIPOS>{}".format(str(title).strip().upper()))

    for l in lines[s + 1 :]:
        if l and l[0] == "@":
            break
        else:
            yield l


def get_mol2_block_lines(title: str, text: str) -> list:
    return list(yield_mol2_block_lines(title, text))


class Atom:
    """
    Symbol, label, atom type.

    ! Additional treatment of stereogenicity.

    """

    def __init__(
        self,
        symbol: str,
        label: str,
        atom_type: str = "C",
        stereo: str = "U",
        ap: bool = False,
    ):

        self.label = str(label)
        self.symbol = str(symbol)
        self.atom_type = str(atom_type)
        self.stereo = stereo  # Maybe use in the future?
        self.ap = ap  # Maybe use in the future?

    def __str__(self):
        return f"{self.label}"

    def __repr__(self):
        return f"{self.label}"

    def set_attachment_point(self, v: bool = True):
        self.ap = v


class Bond:
    """
    Chemical bond
    """

    def __init__(self, a1: Atom, a2: Atom, bond_type: str = None):
        self.a1 = a1
        self.a2 = a2
        self.bond_type = bond_type

    def __contains__(self, a: Atom):
        if self.a1 == a or self.a2 == a:
            return True
        else:
            return False

    # def __eq__(self, other: Bond):
    #     if self.a1 in other or self.a2 in other:
    #         return True
    #     else:
    #         return False

    def __repr__(self):
        return f"Bond ({self.bond_type}) {self.a1}-{self.a2}"


class Fragment:
    """
    This class is not defined yet
    """

    ...


class Property:
    """
    Any property that can be attributable to a whole molecule / atom / bond / conformer
    """

    def __init__(self, ref, ptype: str, value: str, source: str = ""):
        self.ref = ref
        # check reference object for type
        if not isinstance(ref, (Atom, Bond, Fragment, CartesianGeometry, Molecule)):
            raise NotImplementedError(
                "Reference must be an object in (Atom, Bond, Fragment, CartesianGeometry, Molecule)"
            )

        self.ptype = ptype
        self.value = value
        self.source = source


def structure_clone(atoms: List[Atom], bonds: List[Bond]) -> (List[Atom], List[Bond]):
    """
    This functions allows deep copying of atoms and bonds, keeping the connectivity table intact
    Bonds are not kept if their atoms are not in the list!
    """

    old_new_map = {}
    new_atoms = []
    new_bonds = []

    for a in atoms:
        _a = deepcopy(a)
        old_new_map[a] = _a
        new_atoms.append(_a)

    for b in bonds:
        if b.a1 in old_new_map and b.a2 in old_new_map:
            _a1 = old_new_map[b.a1]
            _a2 = old_new_map[b.a2]

            new_bonds.append(Bond(_a1, _a2, b.bond_type))

    return new_atoms, new_bonds, old_new_map

class Orca_Out_Recognize:
    """
    This builds a quick Orca object that is used with the Orca driver
    """
    def __init__(
        self,
        name: str,
        output_file: str,
        calc_type: str,
        hess_file: str,
    ):
        self.name = name
        self.output_file = output_file
        self.calc_type = calc_type
        self.hess_file = hess_file
        self.end_line_list = output_file.split('\n')[-11:]
        self.fixed_err = [f'{x}\n' for x in self.end_line_list]
        self.end_lines = ''.join(self.fixed_err)

        if any("ORCA TERMINATED NORMALLY" in x for x in self.end_line_list):
            self.orca_failed = False
        else:
            self.orca_failed = True

class Molecule:
    """
    Molecule is a class that is supposed to be a cornerstone in all cross-talk between different pieces of code.
    XML-serializable.
    """

    def __init__(
        self,
        name: str,
        atoms: List[Atom],
        bonds: List[Bond],
        geom: CartesianGeometry,
        conformers: List[CartesianGeometry] = [],
        clone: bool = True,
    ):
        """"""
        self.name = name
        self.atoms, self.bonds = atoms, bonds
        self.geom = geom
        self.conformers = conformers

    def __contains__(self, x):
        if isinstance(x, Atom):
            return x in self.atoms
        if isinstance(x, Bond):
            return x in self.bonds

    @classmethod
    def from_mol2(cls, mol2s: str | IOBase, name: str = None):
        """
        Parse a mol2 string and generate a Molecule object
        """

        if isinstance(mol2s, str) and os.path.isfile(mol2s):
            with open(mol2s) as f:
                mol2block = f.read()
        elif hasattr(mol2s, "read"):
            mol2block = mol2s.read()
        elif isinstance(mol2s, str):
            mol2block = mol2s
        else:
            raise NotImplementedError

        if isinstance(mol2block, bytes):
            mol2block = mol2block.decode()

        ## Retrieving molecule metadata
        mol2_header = get_mol2_block_lines("molecule", mol2block)
        _name = mol2_header[0] if name == None else name

        ## Generating the list of atoms and molecular geometry
        mol2_atoms = get_mol2_block_lines("ATOM", mol2block)

        _atoms = []
        _geom = []
        for line in mol2_atoms:
            ls = line.split()
            sym = ls[5].rsplit(".")[0]
            _atoms.append(Atom(sym, ls[1], ls[5]))
            _geom.append(list(map(float, ls[2:5])))

        ## Generating the list of bonds
        mol2_bonds = get_mol2_block_lines("BOND", mol2block)
        _bonds = []

        for line in mol2_bonds:
            ls = line.split()
            a1, a2, bt = int(ls[1]) - 1, int(ls[2]) - 1, ls[3]
            _bonds.append(Bond(_atoms[a1], _atoms[a2], bond_type=bt))

        return cls(
            name=_name, atoms=_atoms, bonds=_bonds, geom=CartesianGeometry(_geom)
        )

    def has_confomers(self):
        """"""
        return True if self.conformers else False

    def get_atom(self, label: Atom | str):
        """
        Returns the first atom with matching label
        """
        # Attempt to do the list index
        return self.atoms[self.get_atom_idx(label)]

    def get_atoms_by_symbol(self, symbol: str):
        """
        Returns all atoms with matching symbol
        """
        atoms = []
        for a in self.atoms:
            if a.symbol == symbol:
                atoms.append(a)

        return atoms

    def get_atom_idx(self, label: Atom | str):
        if isinstance(label, Atom):
            return self.atoms.index(label)

        for a in self.atoms:
            if a.label == label:
                return self.atoms.index(a)

        raise IndexError(f"Cannot find {label} in {self.name}")

    def add_bond(self, a1: Atom, a2: Atom, bond_type: str = "1"):
        """
        Create an additional bond
        """
        self.bonds.append(Bond(a1, a2, bond_type=bond_type))

    def get_bonds_with_atom(self, a: Atom):
        res = []
        for b in self.bonds:
            if a in b:
                res.append(b)
        return res

    def get_atom_valence(self, a: Atom):
        """
        Return sum of bond orders that exist with a given atom
        """
        assert a in self.atoms
        a_bonds = self.get_bonds_with_atom(a)

        orders = []

        for b in a_bonds:
            try:
                o = int(b.bond_type)
            except:
                if b.bond_type == "ar":
                    orders.append(1)
            else:
                orders.append(o)

        valence = sum(orders)
        return valence

    def get_connected_atoms(self, a: Atom | str) -> Set[Atom]:
        """
        Return a `set` of atoms connected to a given atom
        """
        atoms = set()
        for b in self.get_bonds_with_atom(a):
            atoms.add(b.a1 if b.a2 == a else b.a2)

        return atoms

    def get_subgeom(self, atoms: List[Atom], conformer=-1):
        """
        Return a subset of atomic position
        if conformer = -1: return default geometry
        if integer >= 0: return subgeometry of conformer with that index
        """
        subgeom = []
        if conformer == -1:
            for a in atoms:
                c = self.geom.get_coord(self.get_atom_idx(a))
                subgeom.append(c)
        else:
            for a in atoms:
                c = self.conformers[conformer].get_coord(self.get_atom_idx(a))
                subgeom.append(c)

        return CartesianGeometry(subgeom)

    def update_geom_from_xyz(self, xyzblock: str, assert_single=False):
        """
        Update geometry from an xyz block.
        Assert single ensures that the xyz block is a single-molecule xyz file
        """

        coord, atoms, _ = parse_xyz(xyzblock, single=True, assert_single=assert_single)

        if atoms == [x.symbol for x in self.atoms]:
            self.geom.coord = coord
        else:
            raise ValueError(
                "The xyz file signature does not match the current atom list."
            )

    def to_xyz(self, n=-1):
        """
        Return a .xyz block of n-th conformer (if n >= 0), else return default geometry
        """

        if n == -1:
            g = self.geom
        else:
            g = self.conformers[n]

        N = len(self.atoms)
        res = f"{N}\n{self.name}\n"
        for i, a in enumerate(self.atoms):
            x, y, z = g.coord[i]
            res += f"{a.symbol} {x:>10.4f} {y:>10.4f} {z:>10.4f}\n"
        return res

    def get_bond_length(self, b: Bond):
        i1, i2 = self.atoms.index(b.a1), self.atoms.index(b.a2)
        return self.geom.get_distance(i1, i2)

    def get_angle(self, a1: Atom, a2: Atom, a3: Atom):
        "Calculate and return an angle between three atoms. a2 is the middle one"
        i1 = self.get_atom_idx(a1)
        i2 = self.get_atom_idx(a2)
        i3 = self.get_atom_idx(a3)
        return self.geom.get_angle(i1, i2, i3)

    def yield_bonds(self, *b_types):
        """
        Get bonds that match the atom symbols
        b_type is a string of type "Cs-Br"
        """
        for bt in b_types:
            as1, as2 = bt.split("-")
            for b in self.bonds:
                if {as1, as2} == {b.a1.symbol, b.a2.symbol}:
                    yield b

    def fix_geom(
        self, s1: str = "C", s2: str = "C", dist: float = 1.5, center: bool = True
    ):
        """
        Rescale geometry so that the average length of bonds with selected labels is equal to dist.
        Also, translate the geometry to the geometric center if center == True
        """
        dists = []
        for b in self.bonds:
            if set((s1, s2)) == set((b.a1.symbol, b.a2.symbol)):
                dists.append(self.get_bond_length(b))

        if len(dists) == 0:
            # Emergency scenario: just average the bond lengths what we already have
            for b in self.bonds:
                dists.append(self.get_bond_length(b))

        factor = dist / np.average(dists)
        self.geom.scale(factor)

        if center:
            self.geom.center_geom()

    def to_mol2(self):
        """
        Return a .mol2 block
        """
        mol2 = f"@<TRIPOS>MOLECULE\n{self.name}\n{len(self.atoms)} {len(self.bonds)} 0 0 0\nSMALL\nGASTEIGER\n\n"

        mol2 += "@<TRIPOS>ATOM"
        for i, a in enumerate(self.atoms):
            x, y, z = self.geom.coord[i]
            mol2 += f"\n{i+1:>6} {a.label:<3} {x:>10.4f} {y:>10.4f} {z:>10.4f} {a.atom_type:<10} 1 {a.label if a.ap else 'UNL1'} 0.0"

        mol2 += "\n@<TRIPOS>BOND"
        for i, b in enumerate(self.bonds):
            a1, a2 = self.atoms.index(b.a1), self.atoms.index(b.a2)
            mol2 += f"\n{i+1:>6} {a1+1:>6} {a2+1:>6} {b.bond_type:>10}"

        mol2 += "\n\n"

        return mol2

    def embed_conformers(self, *confs: CartesianGeometry, mode="a"):
        """
        This function embeds alternative geometries (conformers)
        if mode == 'a': append conformers to existing list
        if mode == 'w': overwrite the list of conformers
        """
        if mode == "a":
            self.conformers.extend(deepcopy(confs))
        elif mode == "w":
            self.conformers = deepcopy(confs)
        else:
            raise ValueError("Mode can only be 'w' or 'a'")

    def confs_to_multixyz(self):
        labels = [x.symbol for x in self.atoms]
        allxyz = ""
        for i, conf in enumerate(self.conformers):
            xyz = conf.to_xyz(labels, f"{self.name}:{i+1}")
            allxyz += xyz
        return allxyz

    def confs_to_xyzs(self):
        labels = [x.symbol for x in self.atoms]
        allxyz = []
        for i, conf in enumerate(self.conformers):
            xyz = conf.to_xyz(labels, f"{self.name}:{i+1}")
            allxyz.append(xyz)
        return allxyz

    def confs_to_molecules(self, name_fmt="{name}_cf{n}"):
        mols = []
        for cn, cg in enumerate(self.conformers):
            name = name_fmt.format(name=self.name, n=cn)
            m = Molecule(name, self.atoms, self.bonds, cg, clone=True)
            mols.append(m)

        return mols

    def confs_to_mol2_files(self, path="", name_fmt="{name}_cf{n}"):
        """
        This function exports all conformers from current molecule file into mol2 files.

        `path`: directory into which they should be exported.
            Will be created if not existent!

        `name_fmt`: name formatter. Accepts the following variables:
            - `{name}`: molecule name
            - `{n}`: conformer number (conformers will be numcered starting with 0)

        """
        # check if the folder exists and create if needed
        if not os.path.isdir(path):
            os.makedirs(path)

        if self.conformers:
            for m in self.confs_to_molecules(name_fmt=name_fmt):

                fn = os.path.normpath(os.path.join(path, f"{m.name}.mol2"))

                with open(fn, "wt") as f:
                    f.write(m.to_mol2())
        else:
            print(f"{self.name}: no conformers to export")

    def remove_atoms(self, *atoms: Atom):
        """
        Remove selected atoms from the molecule.
        Also deletes all bonds to and between selected atoms, and the respective coordinates.
        """
        for a in atoms:
            for b in self.get_bonds_with_atom(a):
                self.bonds.remove(b)

            aidx = self.get_atom_idx(a)
            self.geom.delete(aidx)
            for conf in self.conformers:
                conf.delete(aidx)

            self.atoms.remove(a)

    @classmethod
    def join(
        cls,
        m1: Molecule,
        m2: Molecule,
        a11: Atom,
        a12: Atom,
        a21: Atom,
        a22: Atom,
        dist: float = 10.0,
    ):
        """
        Join two molecular fragments with bond a11--a21, delete atoms a21 and a22
        """

        if m1.has_confomers() or m2.has_confomers():
            # TODO: implement conformer joining, because that could be pretty powerful
            raise NotImplementedError("Currently conformer joining is not supported")

        ## Pre-flight checks
        if not len(_bm1 := m1.get_bonds_with_atom(a12)) == 1:
            raise SyntaxError(
                f"Attachment poing should only have one bond, found {len(_bm1)}"
            )
        if not len(_bm2 := m2.get_bonds_with_atom(a22)) == 1:
            raise SyntaxError(
                f"Attachment poing should only have one bond, found {len(_bm1)}"
            )

        ## STEP 1. Determine the rotation matrix.
        #   In this algorithm, we rotate `m2`

        i11 = m1.atoms.index(a11)
        i12 = m1.atoms.index(a12)
        i21 = m2.atoms.index(a21)
        i22 = m2.atoms.index(a22)

        v1 = m1.geom.coord[i12] - m1.geom.coord[i11]
        v2 = m2.geom.coord[i22] - m2.geom.coord[i21]

        R = rotation_matrix(v2, -v1)

        ## STEP 2. Transform and translate a copy of geom-2
        g1 = deepcopy(m1.geom)
        g2 = deepcopy(m2.geom)

        g1.set_origin(i11)
        g2.set_origin(i21)

        g2.transform(R)

        # vT = (v1 * (dist - np.linalg.norm(v1) - np.linalg.norm(v2)) /
        #   np.linalg.norm(v1))

        vT = v1 * dist / np.linalg.norm(v1)
        g2.translate(vT)

        ## STEP 3. Start assembling the new class
        g1.delete(i12)
        g2.delete(i22)

        name = f"{m1.name}_{m2.name}"
        atoms = []
        bonds = []

        m1_atoms, m1_bonds, m1_map = structure_clone(m1.atoms, m1.bonds)
        m2_atoms, m2_bonds, m2_map = structure_clone(m2.atoms, m2.bonds)

        for a in m1_atoms + m2_atoms:
            if not a in (m1_map[a12], m2_map[a22]):
                atoms.append(a)

        for b in m1_bonds + m2_bonds:
            if (m2_map[a22] not in b) and (m1_map[a12] not in b):
                bonds.append(b)

        geom = np.concatenate((g1.coord, g2.coord), axis=0)

        result = cls(name=name, atoms=atoms, bonds=bonds, geom=CartesianGeometry(geom))
        result.add_bond(m1_map[a11], m2_map[a21])
        # result.relabel_atoms()

        return result

    @classmethod
    def join_ap(
        cls,
        m1: Molecule,
        m2: Molecule,
        ap1: str = "A0",
        ap2: str = "A1",
        dist: float = 10.0,
    ):
        """
        Join two molecules at the attachment points. Attachment points are defined as atoms with distinct labels, such as A0.
        """

        a12 = m1.get_atom(ap1)
        a22 = m2.get_atom(ap2)

        bonds1 = m1.get_bonds_with_atom(a12)
        bonds2 = m2.get_bonds_with_atom(a22)

        assert (
            len(bonds1) == 1
        ), f"Doesn't look like ap1 is an attachment point: {len(bonds1)} bonds"
        assert (
            len(bonds2) == 1
        ), f"Doesn't look like ap2 is an attachment point: {len(bonds2)} bonds"

        b1 = bonds1[0]
        b2 = bonds2[0]

        a11 = b1.a1 if b1.a2 == a12 else b1.a2
        a21 = b2.a1 if b2.a2 == a22 else b2.a2

        return cls.join(m1, m2, a11, a12, a21, a22, dist=dist)

    @classmethod
    def _get_join_ap_fx(cls, ap1, ap2, dist):
        def fx(mol_tuple):
            return cls.join_ap(mol_tuple[0], mol_tuple[1], ap1=ap1, ap2=ap2, dist=dist)

        return fx

    def bounding_box(self):
        """
        Get the rectangular space that encompasses all atoms in all conformers
        """
        mins = []
        maxs = []

        for g in self.conformers:
            rmin, rmax = g.bounding_box()
            mins.append(rmin)
            maxs.append(rmax)

        rmin = np.min(mins, axis=(0,))
        rmax = np.max(maxs, axis=(0,))

        return rmin, rmax

    def to_xml(self, pretty=True):
        """
        Save the molecule object in an xml format
        """

        xdoc = xmd.Document()

        xdoc.appendChild(xdoc.createComment("MOLLI PACKAGE EXPERIMENTAL XML FORMAT"))

        xmol = xdoc.createElement("molecule")
        xmol.setAttribute("name", self.name)
        xdoc.appendChild(xmol)
        xatoms = xdoc.createElement("atoms")
        xbonds = xdoc.createElement("bonds")
        xgeom = xdoc.createElement("geometry")
        xconfs = xdoc.createElement("conformers")
        xprops = xdoc.createElement("properties")

        for x in xatoms, xbonds, xgeom, xconfs, xprops:
            xmol.appendChild(x)

        ids = {}
        for i, a in enumerate(self.atoms):
            xa = xdoc.createElement("a")
            xa.setAttribute("id", f"{i+1}")
            xa.setAttribute("s", a.symbol)
            xa.setAttribute("t", a.atom_type)
            xa.setAttribute("l", a.label)
            ids[a] = f"{i+1}"
            xatoms.appendChild(xa)

        for i, b in enumerate(self.bonds):
            xb = xdoc.createElement("b")
            xb.setAttribute("id", f"{i+1}")

            id1, id2 = ids[b.a1], ids[b.a2]
            xb.setAttribute("c", f"{id1} {id2}")
            xb.setAttribute("t", b.bond_type)

            xbonds.appendChild(xb)

        xg0 = xdoc.createElement("g")
        xg0.setAttribute("u", "A")
        xg0.setAttribute("t", "cart/3d")
        xg0.appendChild(xdoc.createTextNode(self.geom.dumps()))
        xgeom.appendChild(xg0)

        for i, conf in enumerate(self.conformers):
            cg = xdoc.createElement("g")
            cg.setAttribute("id", f"{i+1}")
            cg.setAttribute("u", "A")
            cg.setAttribute("t", "cart/3d")
            cg.appendChild(xdoc.createTextNode(conf.dumps()))
            xconfs.appendChild(cg)

        if pretty:
            return xdoc.toprettyxml()
        else:
            return xdoc.toxml()

    @classmethod
    def from_xml(cls, fp: str) -> Molecule:
        """
        Parse a molli xml file and create a Molecule instance
        """
        et = xparse(fp)
        # rt = et.getroot()

        mol = et.getroot()
        name = mol.attrib["name"]

        xatoms = mol.findall("./atoms/a")
        xbonds = mol.findall("./bonds/b")
        xgeom = mol.find("./geometry/g")
        xconfs = mol.findall("./conformers/g")
        # xprops = mol.findall("./properties/p")  # not really implemented yet

        atoms = []
        ids = []
        bonds = []
        conformers = []

        for a in xatoms:
            aid, s, l, at = a.attrib["id"], a.attrib["s"], a.attrib["l"], a.attrib["t"]
            ids.append(aid)
            atoms.append(Atom(s, l, at))

        for b in xbonds:
            ia1, ia2 = map(ids.index, b.attrib["c"].split())
            bt = b.attrib["t"]
            bonds.append(Bond(atoms[ia1], atoms[ia2], bt))

        geom = CartesianGeometry.from_str(xgeom.text)
        conformers = [CartesianGeometry.from_str(g.text) for g in xconfs]

        return cls(name, atoms=atoms, bonds=bonds, geom=geom, conformers=conformers)

    @classmethod
    def from_file(cls, fref: str | IOBase):

        # Determine file extension
        if isinstance(fref, str):
            ext = fref.rsplit(".", 1)[1]
        elif hasattr(fref, "name"):
            ext = fref.name.rsplit(".", 1)[1]

        if ext not in ["mol2", "xml"]:
            raise ValueError("Unknown file extension")

        if ext == "mol2":
            return cls.from_mol2(fref)
        if ext == "xml":
            return cls.from_xml(fref)
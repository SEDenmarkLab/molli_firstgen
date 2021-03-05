import molli as ml

mymol = ml.dtypes.Molecule.from_file("tests/bp1.mol2")
xml = mymol.to_xml(pretty=True)

with open("b.xml", "wt") as f:
    f.write(xml)

newmol = ml.dtypes.Molecule.from_xml_file("b.xml")

print(newmol.to_xyz())
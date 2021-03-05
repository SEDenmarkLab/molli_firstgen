import molli as ml

a = ml.dtypes.MolCollection.from_zip("tests/collection.zip")

print(a.name)
print(a.mol_index)
print(len(a.molecules))

a.to_zip("collect-out.zip")


def wm2(m: ml.dtypes.Molecule):
    mol2 = m.to_mol2()
    with open(m.name + ".mol2", "w") as f:
        f.write(mol2)


a.applyfx(wm2)(update=1)

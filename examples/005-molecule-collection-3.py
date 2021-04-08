# Sample code snippet for splitting a cdxml file into a zip file
# And also produce a folder of mol2 files (because why not)

import molli as ml

collection = ml.parsing.split_cdxml("ac1.cdxml", enum=True)

collection.to_multixyz("a.xyz")
# add hydrogens and perform a crude optimization
# returns a list, not collection!

obabel = ml.OpenBabelDriver("ob", scratch_dir="obabel_tmp")
concurrent = ml.Concurrent(
    collection, backup_dir="backup", update=1.0, concurrent=1, timeout=1.0
)

_preoptimized = concurrent(obabel.hadd_opt)()
# recreate them into molecular collections
# preoptimized = ml.Collection("your-files-preoptimized", _preoptimized)

# preoptimized.to_zip("preopt.zip")

# for mol in preoptimized:
#     with open(f"{mol.name}.mol2", "wt") as f:
#         f.write(mol.to_mol2())

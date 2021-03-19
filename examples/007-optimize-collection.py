"""
In this example we will optimize a collection of phosphoramidite-CuCl things
"""
import molli as ml
import os

# loading molecular collection
FILES = os.path.join(os.path.dirname(__file__), "files")
TEMP = os.path.join(os.path.dirname(__file__), "temp")

binols = ml.dtypes.MolCollection.from_zip(f"{FILES}/binol-p-cucl.zip")

#setting up the optimization driver
xtb = ml.drivers.XTBDriver("temp")
opts = binols.applyfx(xtb.minimize)(nprocs=1)

new = ml.dtypes.MolCollection(name="binols-cucl-optimized", molecules=opts)
new.to_zip(f"{FILES}/binol-p-cucl-opt.zip")
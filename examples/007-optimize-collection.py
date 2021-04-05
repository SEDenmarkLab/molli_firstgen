"""
In this example we will optimize a molecule collection
Using an asynchronous aggregator
"""
import molli as ml  # pylint: disable=import-error

collection = ml.Collection.from_zip("tests/ac.zip")
xtb = ml.XTBDriver("myxtb", scratch_dir="scratch/", nprocs=1)
opt = ml.Concurrent(collection, update=1)(xtb.optimize)(method="gfn2", crit="tight")
col_opt = ml.Collection(name="gfnopt-test", molecules=opt)

col_opt.to_zip("tests/ac-gfo.zip")
col_opt.to_multixyz("tests/ac-gfo.xyz")

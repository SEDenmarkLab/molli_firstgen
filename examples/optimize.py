import molli as ml

bp_aa = ml.Collection.from_zip("bp_aa_opt1.zip")

xtb = ml.XTBDriver("sxtb", scratch_dir="/scratch/shvedalx/xtb", nprocs=4)

_opt = ml.Concurrent(bp_aa, scratch_dir="/scratch/shvedalx/c", concurrent=24, update=30)(xtb.optimize)(method='gff', maxiter=100)

bp_aa_opt = ml.Collection("bp_aa_opt", _opt)

bp_aa_opt.to_zip("bp_aa_opt2.zip")


import molli as ml

with ml.utils.WCTimer("Assembling acyclic amines with binol-P cores"):
    # amines-cyclic
    aa = ml.Collection.from_zip('preopt/aa.zip')
    bp = ml.Collection.from_zip('preopt/bp-cu.zip')

    bp_aa = ml.Collection.join(bp, aa, 'A0', 'A1', dist=15, nprocs=128)

bp_aa.to_zip("bp_aa_preopt.zip")

xtb = ml.XTBDriver("sxtb", scratch_dir="/scratch/shvedalx/xtb", nprocs=4)

_opt = ml.Concurrent(bp_aa, scratch_dir="/scratch/shvedalx/c", concurrent=16, update=30)(xtb.fix_long_bonds)(force_const=0.1, in_place=False)

bp_aa_opt = ml.Collection("bp_aa_opt", _opt)

bp_aa_opt.to_zip("bp_aa_opt.zip")


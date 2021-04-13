import molli as ml
from pprint import pprint

bp_aa = ml.Collection.from_zip("bp_aa_opt2.zip")
crest = ml.CRESTDriver("screst", scratch_dir="/scratch/shvedalx/xtb", nprocs=8)

while True:
    try:
        concur = ml.Concurrent(bp_aa, backup_dir="/scratch/shvedalx/conf-backup/", concurrent=16, update=300, timeout=3600)
        _opt = concur(crest.conformer_search)(mdlen=20, mddump=250, method='gfnff', chk_topo=False, constr_val_angles=['P', 'Cu'])
        # _opt = concur(crest.confomer_screen)(method='gfn2')

        bp_aa_opt = ml.Collection("bp_aa_conf", _opt)
        bp_aa_opt.to_zip("bp_aa_conf-pre2.zip")
        
    except:
        pass
    else:
        print("Appears to be done... Congrats boi!")


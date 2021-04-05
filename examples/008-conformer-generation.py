import molli as ml  # pylint: disable=import-error

collection = ml.Collection.from_zip("tests/ac.zip")
crest = ml.CRESTDriver("mycrest", scratch_dir="scratch/", nprocs=4)
ml.Concurrent(collection[8:], concurrent=1, update=30, timeout=3600)(
    crest.conformer_search
)(method="gff")

collection.to_zip("tests/ac-confs.zip")

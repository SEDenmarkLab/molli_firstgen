import molli as ml
from importlib import import_module
from sys import argv, stderr


print("MOLLI", ml.__version__)

if len(argv) > 1:
    wf = argv[1]

    if wf in ml.workflows.__all__:
        print(f"Requested workflow <{wf}>")
        try:
            mod = import_module(f"molli.workflows.{wf}")
            mod.main(argv[2:])
        except:
            print(f"Failed to import {wf}")
            raise
    else:
        print("Error: workflow not recognized.")
        print("List of available workflows:")
        for wfa in ml.workflows.__all__:
            print(f"\t{wfa}")
            
            if wf not in ["help", "--help", "-h"]:
                exit(1)
            else:
                exit(0)

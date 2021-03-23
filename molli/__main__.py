import molli as ml
from sys import argv, stderr


print("=" * 80)
print("\tMOLLI", ml.__version__)
print("=" * 80)


if len(argv) > 1:
    wf = argv[1]

    if wf in ml.workflows.__all__:
        pass
    else:
        print(f"Unable to find workflow <{wf}>")
        print("List of available workflows:")
        for wfa in ml.workflows.__all__:
            print(f"\t{wfa}")
            exit(0)
        if wf not in ["help", "--help", "-h"]:
            exit(1)
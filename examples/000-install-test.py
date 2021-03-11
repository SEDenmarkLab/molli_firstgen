print("Trying to import molli:")
try:
    import molli
except Exception as xc:
    print("- Failed to import molli.")
    print("  Exception:", xc)
    exit(1)
else:
    print("+ Successfully imported molli", molli.__version__)
    exit(0)

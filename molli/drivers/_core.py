import os
from glob import glob
from subprocess import run


class DriverError(Exception):
    "Base class for driver exceptions"
    ...


class ExternalDriver:
    """
        This class defines basic functionality for using external packages as executables
    """
    JOB_ID = 0

    def __init__(self, cwd: str = "/", nprocs=1, encoding='utf8'):
        """
            cwd: working directory from which the external program is called
        """
        # self.cmd = cmd
        self.cwd = cwd
        self.encoding = encoding
        self.nprocs = nprocs
        if not os.path.isdir(cwd):
            os.makedirs(cwd)

    def __call__(self, *args, inp: str = None):
        """
            Low-lever caller
        """
        self.__class__.JOB_ID += 1

        try:
            if inp != None:
                proc = run([str(x) for x in args],
                           capture_output=True,
                           cwd=self.cwd,
                           input=inp,
                           encoding=self.encoding)
            else:
                proc = run([str(x) for x in args],
                           capture_output=True,
                           cwd=self.cwd,
                           encoding=self.encoding)
        except:
            raise DriverError(f"Invalid command: {' '.join(args)}")

        if proc.returncode or len(proc.stdout) == 0:
            raise DriverError(
                f"Process ended in error state.\nCall: {' '.join(map(str, args))}\n Stderr:\n{proc.stderr}"
            )
        else:
            return str(proc.stdout)

    def cleanup(self, regex: str = ''):
        """
            Clean up the working directory
        """
        file_candidates = \
            glob(f"{self.cwd}/{regex}") + \
            glob(f"{self.cwd}/.{regex}")

        for fc in file_candidates:
            if os.path.isfile(fc):
                os.remove(fc)

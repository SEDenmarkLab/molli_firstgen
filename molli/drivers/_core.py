import os
from glob import glob
from subprocess import run
from shutil import rmtree
from warnings import warn
import asyncio as aio

class DriverError(Exception):
    "Base class for driver exceptions"
    ...

class ExternalDriver:
    """
        This class defines basic functionality for using external packages as executables
    """
    JOB_ID = 0

    def __init__(self, scratch_dir: str = "", nprocs=1, encoding='utf8'):
        """
            cwd: scratch directory from which the external program is called
        """
        self.scratch_dir = os.path.normpath(scratch_dir)
        self.encoding = encoding
        self.nprocs = nprocs
    
    def mktemp(self, fmt="{pid}_{jobid}"):
        """ Make a temporary folder """
        
        base = fmt.format(pid = os.getpid(), jobid = self.__class__.JOB_ID)
        dr = os.path.normpath(os.path.join(self.scratch_dir, base))
        os.makedirs(dr, exist_ok=False)
        self.__class__.JOB_ID += 1
        self.cwd = dr

        return dr
    
    def __call__(self, *args, inp: str = None):
        """
            Low-lever caller
        """

        if not hasattr(self, 'cwd'):
            raise DriverError("CWD must be set first. Mktemp")
        
        try:
            if inp != None:
                proc = run([str(x) for x in args],
                           capture_output=True,
                           cwd=self.cwd,
                           input=inp,
                           encoding=self.encoding,
                           env=os.environ
                           )
            else:
                proc = run([str(x) for x in args],
                           capture_output=True,
                           cwd=self.cwd,
                           encoding=self.encoding,
                           env=os.environ
                           )
        except:
            raise DriverError(f"Invalid command: {' '.join(args)}")

        if proc.returncode or len(proc.stdout) == 0:
            raise DriverError(
                f"Process ended in error state.\nCall: {' '.join(map(str, args))}\n Stderr:\n{proc.stderr}"
            )
        else:
            return str(proc.stdout)
    
    async def _async_cmd(self, cmd, fx_default=None, on_exception=None):
        """ 
            This is a coroutine for asyncio. 
            returns (code, stdout, stderr, )
        """
        try:
            proc = await aio.create_subprocess_shell(f'cd {self.cwd};', cmd, stderr=aio.subprocess.PIPE, stdout=aio.subprocess.PIPE)
            code = await proc.wait()
            stdout = await proc.stdout.read()
            stderr = await proc.stderr.read()
        except Exception as xc:
            if callable(on_exception):
                on_exception(xc)
            else:
                raise xc
        else:
            if callable(on_success):
                on_success()
            

    def cleanup(self, sdr):
        """
            Clean up the current working directory
        """
        if os.path.commonpath([sdr, self.scratch_dir]) == self.scratch_dir:
            rmtree(sdr, ignore_errors=True)
        else:
            warn("Could not remove directory", sdr)



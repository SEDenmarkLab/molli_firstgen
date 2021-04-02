import asyncio as aio
from asyncio.subprocess import PIPE
from typing import Callable, Dict, List, Awaitable
from tempfile import TemporaryDirectory as TempDir, NamedTemporaryFile as TempFile
import os

from ..dtypes import CollectionFile, Collection


class AsyncExternalDriver:
    """
    This driver supports asynchronous programming.
    My attempt to bypass multiprocessing.

    The idea is that on machines with a lot of cores we can offload the processing to those cores using OS functionality
    And subsequently treat the problem as if it were I/O bound within the main code (actually it is CPU bound, but that is now OS's problem)

    nprocs is only respected where packages can make use of it
    """

    PREFIX = "molli-ad"

    def __init__(
        self, name="", scratch_dir: str = "", nprocs: int = 1, encoding: str = "utf8"
    ):
        self.scratch_dir = scratch_dir
        self.nprocs = nprocs
        self.encoding = encoding
        self.prefix = f"{self.__class__.PREFIX}.{name}."

    async def aexec(
        self, cmd: str, inp_files: Dict[str, str] = dict(), out_files: List[str] = []
    ):
        """
        Coroutine that asynchronously schedules a shell command to be executed
        Before the command is executed it writes temporary files (`inp_files`)
        After the command is executed, output files are harvested (`out_files`)
        """

        with TempDir(prefix=self.prefix, dir=self.scratch_dir) as td:
            self.writefiles(td, inp_files)
            # Asynchronous process spawning code goes here

            _cmd = f"cd {td}; {cmd}"

            proc = await aio.create_subprocess_shell(_cmd, stdout=PIPE, stderr=PIPE)
            code = await proc.wait()

            _out = await proc.stdout.read()
            _err = await proc.stdout.read()

            stdout = _out.decode(self.encoding)
            stderr = _err.decode(self.encoding)

            files = self.getfiles(td, out_files)

        return code, files, stdout, stderr

    def _exec(self, *args, **kwargs):
        return aio.run(self.aexec(*args, **kwargs))

    @staticmethod
    def writefiles(dr: str, files: Dict[str, str], overwrite=False):
        """
        Write files in directory `d` from text. files dict example:
        {
            "myfile.txt": "Hello, World!\n",
            ...
        }
        """
        for f in files:
            path = os.path.join(dr, f)

            if os.path.isfile(path) and not overwrite:
                raise FileExistsError(path)
            else:
                with open(path, "wt") as fs:
                    fs.write(files[f])

    @staticmethod
    def getfiles(dr: str, files: List[str], strict: bool = False):
        """
        Retrieves files listed in the files list, returns dictionary
        by default (Strict=False), does not raise an exception if file is not found
        """
        result = dict()
        for f in files:
            path = os.path.join(dr, f)
            if not os.path.isfile(path) and strict:
                raise FileNotFoundError(f)
            elif os.path.isfile(path):
                with open(path) as fs:
                    result[f] = fs.read()
        return result


class AsyncAggregator:
    """
    This class provides functionality for processing of molecule collections
    using asynchronous driver awaitables that
    """

    PREFIX = "molli-agg"

    def __init__(
        self,
        collection: Collection | CollectionFile,
        scratch_dir: str = "",
        logfile: str = "",
        dumpfile: str = "",
    ):
        """
        `scratch_dir` is the directory that will be used for temporary file storage

        `logfile` is the file where brief text results of processing will be stored (use tail -f on Linux to monitor if needed)

        `dumpfile` is the new collection that will be dumped as a result of the processing
        """
        pass

    def __call__(self, fx: Callable):
        pass

    async def _worker(self):
        """
        *Non-public coroutine*

        This coroutine processes the queue.
        """
        pass

    async def _add_task(self, coro: Awaitable, **kwargs):
        """
        *Non-public coroutine*

        This coroutine defines a task that needs to be performed
        kwargs define all parameters that are passed alongside a molecule
        """
        pass

    async def _construct_queue(self):
        pass

    async def _process_monitor(self):
        pass

    def __enter__(self):
        pass

    def __exit__(self, *args):
        pass

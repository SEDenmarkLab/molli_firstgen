from __future__ import annotations
import asyncio as aio
from asyncio.subprocess import PIPE
from typing import Callable, Dict, List, Awaitable
from tempfile import TemporaryDirectory as TempDir, NamedTemporaryFile as TempFile
import functools
import os
from datetime import datetime
from tempfile import mkstemp
import pickle
from ..dtypes import CollectionFile, Collection, Molecule


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


class AsyncConcurrent:
    """
    This class should be used as a decorator that allows to define a function on a molecule object,
    but apply to a collection
    """

    PREFIX = "molli-acc"

    def __init__(
        self,
        collection: Collection,
        /,
        backup_dir: str = "",
        logfile: str = "out.log",
        update: float = 1.0,
        timeout: float = 600.0,
        concurrent=4,
    ):
        """
        `scratch_dir` is the directory that will be used for temporary file storage     
        `logfile` is the file where brief text results of processing will be stored (use tail -f on Linux to monitor if needed)     
        `dumpfile` is the new collection that will be dumped as a result of the processing      
        """
        self.collection = collection
        self.backup_dir = backup_dir
        self.logfile = logfile
        self.concurrent = concurrent
        # self._queue: aio.Queue
        # self._construct_queue(self.collection)
        self.update = update
        self.timeout = timeout
        self._result = [None] * len(collection)

    def _construct_queue(self, collection: Collection):
        self._queue = aio.Queue()
        for i, m in enumerate(collection):
            self._queue.put_nowait((i, m))

    async def _worker(self, fx: Awaitable):
        while True:
            i, mol = await self._queue.get()
            try:
                res = await aio.wait_for(fx(mol), timeout=self.timeout)
            except Exception as xc:
                res = xc
                with open(self.logfile, "at") as f:
                    f.write(f"Error at {mol.name}", str(xc))
            finally:
                self._queue.task_done()
                self._result[i] = res

                if isinstance(res, Molecule):
                    _, fn = mkstemp(prefix=mol.name, suffix='.xml', dir=self.backup_dir)
                    with open(fn, "wt") as f:
                        f.write(res.to_xml())

                    
                    

    def _spawn_workers(self, fx: Awaitable):
        if hasattr(self, "_worker_pool"):
            raise Exception("")
        else:
            self._worker_pool = []
            for _ in range(self.concurrent):
                task = aio.create_task(self._worker(fx))
                self._worker_pool.append(task)

    def _cancel_workers(self):
        for w in self._worker_pool:
            w.cancel()

    def get_status(self):

        total = len(self.collection)
        success = 0
        timed_out = 0
        other_err = 0
        not_started = 0

        for x in self._result:
            if isinstance(x, Exception) and isinstance(x, aio.TimeoutError):
                timed_out += 1
            elif isinstance(x, Exception):
                other_err += 1
            elif x == None:
                not_started += 1
            else:
                success += 1

        return total, success, timed_out, other_err, not_started
        # return f"{datetime.now()} : {success:>7}(++) {timed_out:>7}(??) {other_err:>7}(**) {not_started:>7}(..)"

    async def aexec(self, fx: Awaitable):
        self._queue: aio.Queue
        self._construct_queue(self.collection)
        self._spawn_workers(fx)
        
        start = datetime.now()
        print(f"Starting to process {self._queue.qsize()} molecules:")

        while True:
            try:
                await aio.wait_for(self._queue.join(), timeout=self.update)
            except aio.TimeoutError:
                pass
            except aio.CancelledError:
                print("One or more of the workers was cancelled. Exiting...")
                self._cancel_workers()
                break
            except KeyboardInterrupt:
                print("User interrupted the script. Exiting...")
                self._cancel_workers()
                break
            else:
                break
            finally:
                total, success, timed_out, other_err, not_started = self.get_status()
                s = success
                sr = s / total
                e = timed_out + other_err
                er = e / total
                if (sr + er) > 0:
                    eta = start + (datetime.now() - start) / (sr + er)
                else:
                    eta = "..."

                print(
                    f"{datetime.now() - start} --- successful {sr:>6.1%} ({s:>7}) --- failed {er:>6.1%} ({e:>7}) --- ETA {eta}"
                )
                  
        self._cancel_workers()
        del self._queue

        return self._result
     

    def _exec(self, fx: Awaitable):
        return aio.run(self.aexec(fx))

    def __call__(self, fx: Awaitable):
        """
        This is a decorator
        """

        def inner(**kwargs):
            """
            `kwargs` are applied to `fx` !
            """

            async def f(m):
                return await fx(m, **kwargs)

            return self._exec(f)

        return inner

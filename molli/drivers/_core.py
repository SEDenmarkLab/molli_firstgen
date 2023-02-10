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
from dataclasses import dataclass
from ..dtypes.molecule import Orca_Out_Recognize
from glob import glob


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

        if not os.path.isdir(scratch_dir):
            os.makedirs(scratch_dir)

    async def aexec(
        self, cmd: str, inp_files: Dict[str, str] = dict(), out_files: List[str] = []
    ):
        """
        Coroutine that asynchronously schedules a shell command to be executed
        Before the command is executed it writes temporary files (`inp_files`)
        After the command is executed, output files are harvested (`out_files`)
        returns: code, files, stdout, stderr
        """

        with TempDir(prefix=self.prefix, dir=self.scratch_dir) as td:
            self.writefiles(td, inp_files)
            # Asynchronous process spawning code goes here

            _cmd = f"cd {td}; {cmd}"

            with open(f"{td}/_output", "wt") as f: 
                proc = await aio.create_subprocess_shell(_cmd, stdout=f, stderr=PIPE)
                code = await proc.wait()
                _err = await proc.stderr.read()

            with open(f"{td}/_output") as f:
                stdout = f.read()
            
            stderr = _err.decode(self.encoding)

            files = self.getfiles(td, out_files)

            if code:
                td_base = os.path.basename(td)
                with open(f"{self.scratch_dir}/{td_base}_dump_stdout.log", "wt") as f:
                    f.write(stdout)
                with open(f"{self.scratch_dir}/{td_base}_dump_stderr.log", "wt") as f:
                    f.write(stderr)
                    for fout in inp_files:
                        f.write(f"\n\n === {fout} ===\n\n")
                        f.write(inp_files[fout])

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
                if '.gbw' in path:
                    with open(path, "rb") as fs:
                        result[f] = fs.read()
                else:
                    with open(path, "rt") as fs:
                        result[f] = fs.read()
        return result

@dataclass(repr=True, init=True, frozen=True)
class OrcaJobDescriptor:
    '''
    Used to avoid opening entire out file
    '''
    mol_name: str
    out_name: str
    failed: bool
    calc_type: str 
    end_lines: str
    hess_file_name: str
    gbw_file_name: str

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
        qsize=None,
        concurrent=4,
    ):
        """
        `scratch_dir` is the directory that will be used for temporary file storage
        `logfile` is the file where brief text results of processing will be stored (use tail -f on Linux to monitor if needed)
        `dumpfile` is the new collection that will be dumped as a result of the processing
        `reset`: reset the workers every so often to make sure that there are no unexpected crashes
        """
        self.collection = collection
        self.backup_dir = backup_dir
        self.logfile = logfile
        self.concurrent = concurrent
        # self._queue: aio.Queue
        # self._construct_queue(self.collection)
        self.update = update
        self.timeout = timeout
        self.qsize = qsize
        self._result = [None] * len(collection)
        self._bypassed = 0

        if not os.path.isdir(self.backup_dir):
            os.makedirs(self.backup_dir)

    def _construct_queue(self, collection: Collection):
        self._queue = aio.Queue()
        for i, m in enumerate(collection):
            # THIS IS A STUB CODE FOR RESTARTING CALCULATIONS
            if backed_up := glob(f"{self.backup_dir}/{m.name}.*.xml"):
                # print(f"Molecule {m.name} has already been done. trying to read the file ", end='')
                try:
                    res = Molecule.from_xml(backed_up[0])
                except:
                    # print("... not good. Queuing up.")
                    self._queue.put_nowait((i, m))
                else:
                    # print("... success! Bypassing.")
                    if len(res.conformers):
                        self._result[i] = res
                        self._bypassed += 1
                    else:
                        self._queue.put_nowait((i, m))
            elif backed_up := glob(f"{self.backup_dir}/{m.name}*.out"):
                try:
                    if len(backed_up) >= 1:
                        orca_failed = True
                        for out_file in backed_up:
                            general_file_path_list = out_file.split('.')
                            #Returns the most recent calculation done on this
                            calc_type = general_file_path_list[0].split('_')[-1]
                            print(general_file_path_list)
                            with open(out_file, 'r') as f:
                                end_line_list = f.readlines()[-11:]
                                #This is a check to see if the ORCA output terminated normally
                                #In the future, this should be further elaborated upon based on different erroring out messages, especially those regarding SCF non-convergence
                                if any('ORCA TERMINATED NORMALLY' in x for x in end_line_list):
                                    orca_failed = False
                                    fixed_err = [f'{x}\n' for x in end_line_list]
                                    end_lines = ''.join(fixed_err)
                                    print(out_file)
                                    #Creates general file path not including hash
                                    general_file_path = f"{general_file_path_list[0]}."
                                    backed_up = [out_file]
                                    break
                        if orca_failed:
                            raise FileNotFoundError('ORCA DID NOT TERMINATE NORMALLY')
                    else:
                        raise FileNotFoundError('Output File not found! Restarting Calculation')
                    
                    assert len(backed_up) == 1, 'The length of backed_up must be equal to 1.' 
                    gbw_file_name = None

                    if len(gbw_back := glob(f"{general_file_path}*.gbw")) == 1:
                        gbw_file_name = gbw_back[0]
                    else:
                        raise FileNotFoundError("GBW file not found in calculation, Calculation must be restared. Restarting...")


                    hess_file_name = None
                    #Checks to see if frequency was the most recent calculation run
                    if '_freq.' in general_file_path:
                        #Glob is used to scrape any files and use a special character "*" to allow for collection of all files
                        if len(hess_back :=glob(f"{general_file_path}*.hess")) == 1:
                            hess_file_name = hess_back[0]
                        elif len(hess_back) < 1:
                            raise FileNotFoundError('Hessian not found for FREQ calculation. FREQ calculation must be restarted. Restarting...')
                        else:
                            raise NameError('Multiple files found with the general file name, unclear which one is associated with which file, restarting calculation')
                except Exception as e:
                    # print("... not good. Queuing up.")
                    print(e)
                    self._queue.put_nowait((i, m))
                else:
                    # print("... success! Bypassing.")  
                    #This returns the result
                    self._result[i] = OrcaJobDescriptor(mol_name =m.name, out_name = backed_up[0], failed=orca_failed, calc_type = calc_type, end_lines = end_lines, hess_file_name = hess_file_name, gbw_file_name=gbw_file_name)
                    self._bypassed += 1

            else:
                self._queue.put_nowait((i, m))


        print(
            f"=== REQUESTED {len(collection)} :: IN QUEUE {self._queue.qsize()} :: BYPASSED {self._bypassed} ===",
            flush=True,
        )

    async def _worker(self, fx: Awaitable):
        while True:
            try:
                i, mol = self._queue.get_nowait()
            except:
                break

            try:
                res = await aio.wait_for(fx(mol), timeout=self.timeout)
            except KeyboardInterrupt:
                raise KeyboardInterrupt
            except aio.TimeoutError as xc:
                res = xc
            except Exception as xc:
                res = xc
            finally:
                self._queue.task_done()
                self._result[i] = res

                if isinstance(res, Molecule):
                    fd, fn = mkstemp(
                        prefix=f"{res.name}.", suffix=".xml", dir=self.backup_dir
                    )
                    with open(fn, "wt") as f:
                        f.write(res.to_xml())

                    os.close(fd)
                elif isinstance(res, Orca_Out_Recognize):
                    if res.output_file is not None:
                        out_fd, out_fn = mkstemp(
                            prefix=f"{res.name}_{res.calc_type}.", suffix=".out", dir=self.backup_dir
                        )
                        with open(out_fn, "wt") as f:
                            f.write(res.output_file)
                            # out_file_name = out_fn
                        os.close(out_fd)
                    else:
                        out_fn = None
                    if res.hess_file is not None:
                        hess_fd, hess_fn = mkstemp(
                        prefix=f"{res.name}_{res.calc_type}.", suffix=".hess", dir=self.backup_dir
                        )
                        with open(hess_fn, "wt") as f:
                            f.write(res.hess_file)
                        os.close(hess_fd)
                    else:
                        hess_fn = None
                    if res.gbw_file is not None:
                        gbw_fd, gbw_fn = mkstemp(
                        prefix=f"{res.name}_{res.calc_type}.", suffix=".gbw", dir=self.backup_dir
                        )
                        with open(gbw_fn, "wb") as f:
                            f.write(res.gbw_file)
                        os.close(gbw_fd)
                    else:
                        gbw_fn = None

                    #Returning tuple as result including name, if orca failed, calc type, and the last 11 lines
                    self._result[i] = OrcaJobDescriptor(mol_name=res.name,out_name = out_fn ,failed = res.orca_failed,  calc_type = res.calc_type, end_lines = res.end_lines, hess_file_name = hess_fn, gbw_file_name=gbw_fn)
                ####Unblake Dev Shit Again
                    
    def _spawn_workers(self, fx: Awaitable, n=1):
        if hasattr(self, "_worker_pool"):
            raise Exception("")
        else:
            self._worker_pool = []

        for _ in range(n):
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
            elif x is None:
                not_started += 1
            else:
                success += 1

        return total, success, timed_out, other_err, not_started
        # return f"{datetime.now()} : {success:>7}(++) {timed_out:>7}(??) {other_err:>7}(**) {not_started:>7}(..)"

    async def aexec(self, fx: Awaitable):
        self._queue: aio.Queue
        self._construct_queue(self.collection)
        self._spawn_workers(fx, self.concurrent)

        start = datetime.now()
        print(f"Starting to process {self._queue.qsize()} molecules:")

        while True:
            try:
                await aio.wait_for(self._queue.join(), timeout=self.update)

            except aio.TimeoutError:
                # check if tasks are doing okay
                for w in self._worker_pool:
                    if w.cancelled() and self._queue.qsize():
                        print("Ghost worker was replaced with a newly spawned one.")
                        await aio.sleep(10)
                        self._spawn_workers(fx, n=1)
            else:
                break
            finally:
                
                total, success, timed_out, other_err, not_started = self.get_status()
                b = self._bypassed
                br = self._bypassed / total

                if total - b == 0:
                    break

                s = success - b
                try:
                    sr = s / (total - b)
                except:
                    break
                e = timed_out + other_err
                er = e / (total - b)

                if (sr + er) > 0:
                    eta = start + (datetime.now() - start) / (sr + er)
                else:
                    eta = "..."

                print(
                    f"{datetime.now() - start} --- successful {sr:>6.1%} ({s:>7}) --- failed {er:>6.1%} ({e:>7}) --- ETA {eta} WORKERS {self.concurrent}",
                    flush=True,
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

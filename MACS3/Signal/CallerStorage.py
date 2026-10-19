# cython: language_level=3
"""Internal, caller-owned pileup storage; no scoring or caller imports.

The parent owns the directory until destroy(). Forked workers inherit its path
and write private copies of the chromosome-file mapping. Returned paths are
registered by the parent; directory-wide cleanup also finds unreturned files.
Destroy must run only after execution scopes have joined every worker. It is
idempotent and the same owner can allocate a fresh directory on the next call.
"""

import os
from stat import S_ISREG
import struct
from tempfile import mkdtemp, mkstemp

import cython
from cython.cimports.posix.unistd import close
from cython.cimports.libc.stdlib import calloc, free
import numpy as np

from MACS3.Utilities.Execution import resolve_execution_policy
from MACS3.Utilities.Performance import DEFAULTS
from MACS3.Utilities.Threads import map_chromosomes


# header of a pileup file: the dtype and length of each of the three
# arrays [pos, treat, ctrl], which follow it as raw bytes
PILEUP_FILE_HEADER = struct.Struct("<8sq8sq8sq")


def write_pileup_file(fd: int, arrays: list):
    """Write a chromosome's [pos, treat, ctrl] arrays to the open file
    descriptor ``fd`` and close it."""
    header: list = []

    with open(fd, "wb") as f:
        for a in arrays:
            if a.ndim != 1 or a.dtype.hasobject:
                raise TypeError("pileup arrays must be 1-D and numeric")
            header += [a.dtype.str.encode(), a.shape[0]]
        f.write(PILEUP_FILE_HEADER.pack(*header))
        for a in arrays:
            f.write(a)
    return


def read_pileup_file(filename) -> list:
    """Read back the [pos, treat, ctrl] arrays written by
    write_pileup_file."""
    arrays: list = []

    with open(filename, "rb") as f:
        header = PILEUP_FILE_HEADER.unpack(f.read(PILEUP_FILE_HEADER.size))
        for i in range(0, 6, 2):
            a = np.empty(header[i + 1], dtype=header[i].rstrip(b"\0").decode())
            if f.readinto(a) != a.nbytes:
                raise EOFError("pileup file %s is truncated" % filename)
            arrays.append(a)
    return arrays


# PileupStorage.destroy keeps a descriptor on at most this many
# unlinked pileup files, those of at least this many bytes, for a
# thread to close
CLOSE_LATER_MAX_FILES = cython.declare(
    cython.Py_ssize_t, DEFAULTS.deferred_close_max_files)
CLOSE_LATER_MIN_BYTES = cython.declare(
    cython.longlong, DEFAULTS.deferred_close_min_bytes)


def close_unlinked_files(fds: list):
    """Close ``fds``, the last references to pileup files that are
    already unlinked; the kernel frees each file's pages in this
    close. The closes run without the GIL, so that a thread running
    this one waits for the GIL once, not once per file."""
    n: cython.Py_ssize_t = len(fds)
    i: cython.Py_ssize_t
    c_fds: cython.pointer(cython.int)

    c_fds = cython.cast(cython.pointer(cython.int),
                        calloc(n + 1, cython.sizeof(cython.int)))
    if c_fds == cython.NULL:
        raise MemoryError()
    for i in range(n):
        c_fds[i] = fds[i]
    with cython.nogil:
        for i in range(n):
            close(c_fds[i])
    free(c_fds)
    return


def close_unlinked_files_later(fds: list, sizes: list, execution_policy=None):
    """Close unlinked files within the budget, joining before returning.

    The compatibility name is retained. At budget one no thread starts.
    Unclaimed descriptors are closed here even if worker startup fails.
    """
    policy = resolve_execution_policy(execution_policy=execution_policy)
    if not fds:
        return
    n = min(policy.threads, len(fds))
    groups = [[] for _ in range(n)]
    loads = [0] * n
    for j in sorted(range(len(fds)), key=sizes.__getitem__, reverse=True):
        k = loads.index(min(loads))
        groups[k].append(fds[j])
        loads[k] += sizes[j]

    def finish(k):
        group = groups[k]
        # Ownership transfers to this invocation before native close releases
        # the GIL; the finalizer handles only groups never dispatched.
        groups[k] = None
        try:
            close_unlinked_files(group)
        except MemoryError:
            for fd in group:
                os.close(fd)

    try:
        map_chromosomes(finish, list(range(n)), loads,
                        execution_policy=policy)
    finally:
        for k in range(n):
            if groups[k] is not None:
                finish(k)


class PileupStorage:
    """One caller's files, reused until explicit or failure-triggered destroy.

    Only pileup caching is best effort: ordinary write failures discard the
    incomplete file and leave scoring arrays usable. BaseException failures
    such as KeyboardInterrupt propagate after the partial file is unlinked.
    The writer owns the mkstemp descriptor even when it raises.
    """

    def __init__(self, execution_policy):
        self.policy = execution_policy
        self.files = {}
        self.directory = None

    def ensure(self):
        if self.directory is None:
            self.directory = mkdtemp(prefix="macs3-pileup-")

    def read(self, chrom):
        if chrom in self.files:
            try:
                return read_pileup_file(self.files[chrom])
            except Exception:
                pass
        return None

    def write(self, chrom, arrays):
        # Preserve descriptor-based writes: reopening/truncating a path can
        # force ext4 to flush the file to disk when it is closed.
        fd, filename = mkstemp(dir=self.directory)
        self.files[chrom] = os.fsencode(filename)
        try:
            write_pileup_file(fd, arrays)
        except BaseException as exc:
            os.unlink(filename)
            self.files.pop(chrom, None)
            if not isinstance(exc, Exception):
                raise

    def destroy(self):
        """Remove the caller's owned pileup directory, including worker files
        whose names were never returned. Cached paths are cleared so calling
        again recreates storage using the standard tempfile location.

        Each file is unlinked here, so the temporary directory is left
        as before. Freeing a file's pages (2.75 GB for a large control)
        is the slow part, and the kernel does it when the last
        reference to the file goes: a descriptor opened before the
        unlink holds that reference. Budgeted threads close the descriptors
        and are joined here before the next phase can start. At budget one
        the calling thread closes every descriptor.
        """
        fds: list = []
        sizes: list = []

        try:
            for name in (os.listdir(self.directory)
                         if self.directory is not None else []):
                f = os.fsencode(os.path.join(self.directory, name))
                # os.path.isfile(f), keeping the size
                try:
                    st = os.stat(f)
                except (OSError, ValueError):
                    continue
                if not S_ISREG(st.st_mode):
                    continue
                # a bounded number of descriptors, on files large
                # enough to be worth it
                if (st.st_size >= CLOSE_LATER_MIN_BYTES and
                        len(fds) < CLOSE_LATER_MAX_FILES):
                    try:
                        fds.append(os.open(f, os.O_RDONLY | os.O_CLOEXEC))
                        sizes.append(st.st_size)
                    except OSError:
                        pass
                os.unlink(f)
        except BaseException:
            close_unlinked_files(fds)
            raise
        if fds:
            close_unlinked_files_later(fds, sizes, self.policy)
        if self.directory is not None:
            os.rmdir(self.directory)
            self.directory = None
        self.files.clear()
        return

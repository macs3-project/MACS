"""Per-invocation CPU budgets, independent of third-party thread pools.

Resolve at the API/command boundary and pass the immutable policy to consumers.
Resource discovery is a snapshot; it never changes affinity or environment.
"""

from dataclasses import dataclass
import multiprocessing
import os
from pathlib import Path, PurePosixPath
import re
import sys
import threading

from MACS3.Utilities.Performance import DEFAULTS


UNSET = object()
_MAX_COUNT = (1 << 31) - 1


def _positive_count(value, name):
    if type(value) is not int or not 0 < value <= _MAX_COUNT:
        raise ValueError(f"{name} must be a positive integer <= {_MAX_COUNT}")
    return value


def _read(path):
    try:
        return Path(path).read_text()
    except (OSError, UnicodeError):
        return None


def _unescape_mount(value):
    return re.sub(r"\\([0-7]{3})", lambda m: chr(int(m[1], 8)), value)


def _quota_limits():
    """Visible cgroup CPU bandwidth limits, including mount ancestors."""
    memberships = _read('/proc/self/cgroup')
    mounts = _read('/proc/self/mountinfo')
    if memberships is None or mounts is None:
        return [], ['cgroup metadata unavailable']
    groups = {}
    for line in memberships.splitlines():
        fields = line.split(':', 2)
        if len(fields) == 3:
            for controller in fields[1].split(','):
                groups[controller] = fields[2]
    limits, notes = [], []
    seen = set()
    for line in mounts.splitlines():
        fields = line.split()
        try:
            sep = fields.index('-')
            fs = fields[sep + 1]
            controllers = fields[sep + 3].split(',')
            if fs == 'cgroup2':
                member = groups.get('')
            elif fs == 'cgroup' and 'cpu' in controllers:
                member = groups.get('cpu')
            else:
                continue
            if member is None:
                continue
            root = PurePosixPath(_unescape_mount(fields[3]))
            member = PurePosixPath(member)
            mount = Path(_unescape_mount(fields[4]))
            if (not member.is_absolute() or '..' in member.parts
                    or not root.is_absolute() or not mount.is_absolute()):
                continue
            # In a cgroup namespace membership can be relative to a hidden
            # ancestor, whereas mountinfo's root still names that ancestor.
            relative = (member.relative_to(root) if member.is_relative_to(root)
                        else member.relative_to('/'))
            current = mount.joinpath(*relative.parts)
        except (ValueError, IndexError):
            continue
        while True:
            key = (fs, str(current))
            if key not in seen:
                seen.add(key)
                if fs == 'cgroup2':
                    raw = _read(current / 'cpu.max')
                    values = raw.split() if raw is not None else []
                else:
                    quota = _read(current / 'cpu.cfs_quota_us')
                    period = _read(current / 'cpu.cfs_period_us')
                    values = [quota, period] if quota and period else []
                try:
                    if len(values) != 2:
                        raise ValueError()
                    period = int(values[1])
                    if period <= 0:
                        raise ValueError()
                    unlimited = (values[0] == 'max' if fs == 'cgroup2'
                                 else int(values[0]) < 0)
                    if not unlimited:
                        quota = int(values[0])
                        if quota <= 0:
                            raise ValueError()
                        limits.append((f'{fs}:{current}', max(1, quota // period)))
                except ValueError:
                    notes.append(f'quota unavailable: {current}')
            if current == mount:
                break
            current = current.parent
    return limits, notes


@dataclass(frozen=True, slots=True)
class CPUAvailability:
    """Detected whole-CPU capacity and the signals used to obtain it."""

    count: int
    limits: tuple
    notes: tuple

    def __post_init__(self):
        _positive_count(self.count, 'available CPUs')
        if (type(self.limits) is not tuple or
                any(type(item) is not tuple or len(item) != 2 or
                    not isinstance(item[0], str) or type(item[1]) is not int or
                    item[1] <= 0 for item in self.limits)):
            raise ValueError('CPU limits must be immutable (name, positive count) pairs')
        if type(self.notes) is not tuple or any(not isinstance(x, str) for x in self.notes):
            raise ValueError('CPU detection notes must be an immutable tuple of strings')


def detect_cpu_resources():
    """Combine host, process, affinity and visible Linux cgroup limits."""
    limits, notes = [], []
    for name in ('process_cpu_count', 'cpu_count'):
        function = getattr(os, name, None)
        if function is not None:
            try:
                value = function()
                if type(value) is int and value > 0:
                    limits.append((name, value))
            except (OSError, NotImplementedError):
                notes.append(f'{name} unavailable')
    affinity = getattr(os, 'sched_getaffinity', None)
    if affinity is not None:
        try:
            count = len(affinity(0))
            if count:
                limits.append(('affinity', count))
        except (OSError, NotImplementedError):
            notes.append('affinity unavailable')
    if sys.platform == 'linux':
        quotas, quota_notes = _quota_limits()
        limits.extend(quotas)
        notes.extend(quota_notes)
    count = min([value for _, value in limits] + [_MAX_COUNT]) if limits else 1
    return CPUAvailability(count, tuple(limits), tuple(notes))


def available_cpus():
    """Return at least one CPU, tolerating unsupported OS resource queries."""
    return detect_cpu_resources().count


@dataclass(frozen=True, slots=True)
class ExecutionPolicy:
    """Validated budget including the caller; child processes use B-1 slots."""

    requested_threads: int
    threads: int
    availability: CPUAvailability
    source: str = 'api'
    requested_process_backend: str = 'auto'

    def __post_init__(self):
        _positive_count(self.requested_threads, 'threads')
        _positive_count(self.threads, 'effective threads')
        if not isinstance(self.availability, CPUAvailability):
            raise TypeError('availability must be CPUAvailability')
        _positive_count(self.availability.count, 'available CPUs')
        if self.threads != min(self.requested_threads, self.availability.count):
            raise ValueError('effective threads must equal the clamped request')
        if self.requested_process_backend not in ('auto', 'serial'):
            raise ValueError("process_backend must be 'auto' or 'serial'")
        if not isinstance(self.source, str):
            raise TypeError('policy source must be a string')

    @property
    def process_backend(self):
        return 'serial' if self.threads == 1 else self.requested_process_backend

    @property
    def process_workers(self):
        return self.threads - 1 if self.process_backend == 'auto' else 0

    def describe(self):
        return (f'Execution: requested={self.requested_threads}, source={self.source}, '
                f'available={self.availability.count}, threads={self.threads} '
                f'(including caller), process_backend={self.process_backend}, '
                f'max_process_workers={self.process_workers}')


def resolve_execution_policy(threads=None, process_backend=UNSET, *,
                             execution_policy=None, default_threads=None):
    """Explicit policy/count wins over MACS3_NUM_THREADS, then default eight.

    Supplied policies are never re-resolved, including after environment or
    resource changes. Explicit conflicting overrides are rejected.
    """
    if process_backend is not UNSET and process_backend not in ('auto', 'serial'):
        raise ValueError("process_backend must be 'auto' or 'serial'")
    if execution_policy is not None:
        if not isinstance(execution_policy, ExecutionPolicy):
            raise TypeError('execution_policy must be ExecutionPolicy')
        if threads is not None:
            _positive_count(threads, 'threads')
            if threads != execution_policy.requested_threads:
                raise ValueError('threads conflicts with execution_policy')
        if (process_backend is not UNSET and
                process_backend != execution_policy.requested_process_backend):
            raise ValueError('process_backend conflicts with execution_policy')
        return execution_policy
    source = 'explicit'
    if threads is None:
        value = os.environ.get('MACS3_NUM_THREADS')
        if value is not None:
            if not value.isascii() or not value.isdecimal():
                raise ValueError('MACS3_NUM_THREADS must be a positive integer')
            try:
                threads = int(value)
            except ValueError:
                raise ValueError('MACS3_NUM_THREADS must be a positive integer') from None
            source = 'MACS3_NUM_THREADS'
        else:
            threads = DEFAULTS.cpu_limit if default_threads is None else default_threads
            source = 'default'
    _positive_count(threads, 'threads')
    availability = detect_cpu_resources()
    return ExecutionPolicy(threads, min(threads, availability.count), availability,
                           source, 'auto' if process_backend is UNSET else process_backend)


def process_backend_reason(process_backend):
    """Return None when fork is eligible now, otherwise the serial reason.

    Recheck immediately before creating a pool. Native threads/concurrent
    thread creation cannot be detected reliably; embedded hosts may require
    an explicit serial policy.
    """
    if process_backend != 'auto':
        return 'serial policy'
    if sys.platform != 'linux':
        return 'automatic fork requires Linux'
    if threading.current_thread() is not threading.main_thread():
        return 'not the main Python thread'
    if threading.active_count() != 1:
        return 'other Python threads are active'
    if multiprocessing.current_process().daemon:
        return 'daemon process'
    if 'fork' not in multiprocessing.get_all_start_methods():
        return 'fork is unavailable'
    return None

"""Target-interpreter probes for benchmark_performance.py.

Executed from a temporary working directory. Does not add the checkout to
sys.path: the selected interpreter/PYTHONPATH determines the installation.
"""

import gc
import hashlib
import importlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import resource
import runpy
import sys
import sysconfig
import time
import weakref


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def peak_rss():
    scale = 1 if sys.platform == 'darwin' else 1024
    return {'self_peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale,
            'largest_child_peak_rss_bytes': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * scale}


def constrain(spec):
    if spec.get('affinity'):
        os.sched_setaffinity(0, spec['affinity'])
    if spec.get('address_space_mib'):
        amount = spec['address_space_mib'] * 1024 * 1024
        _, hard = resource.getrlimit(resource.RLIMIT_AS)
        resource.setrlimit(resource.RLIMIT_AS, (amount, hard))


def metadata(spec):
    import MACS3
    import numpy as np
    from MACS3.Utilities.Constants import MACS_VERSION
    modules = {}
    for name in ('IO.Parser', 'Signal.Prob', 'Signal.PairedEndTrack', 'Signal.CallPeakUnit'):
        path = importlib.import_module('MACS3.' + name).__file__
        if not path.endswith(('.so', '.pyd')):
            raise RuntimeError(f'Benchmark requires compiled extensions: {path}')
        modules[name] = path
    try:
        from numpy._core.multiarray import get_handler_name
    except ImportError:
        from numpy.core.multiarray import get_handler_name
    from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII
    handlers, cleanup_notes = {}, {}
    for cls in (PETrackI, PETrackII):
        track = cls(buffer_size=17)
        starts = np.arange(300000, dtype='i4')
        if cls is PETrackI:
            track.add_loc_arrays(b'chr1', starts, starts + 150)
        else:
            track.barcode_id(b'BC')
            track.add_loc_arrays(b'chr1', starts, starts + 150,
                                 np.ones(len(starts), dtype='u2'),
                                 np.zeros(len(starts), dtype='i4'))
        array = track.get_locations_by_chr(b'chr1')
        handlers[cls.__name__] = get_handler_name(array)
        del array
        try:
            track.destroy()
        except AttributeError as exc:
            # Both fixed historical refs have this unrelated public cleanup
            # bug. Preserve the diagnostic; this isolated process exits after
            # metadata. Actual benchmark CLI failures still fail their run.
            if (cls is not PETrackII or
                    importlib.util.find_spec('MACS3.Utilities.Execution') is not None or
                    str(exc) != "'dict' object has no attribute 'resize'"):
                raise
            cleanup_notes[cls.__name__] = str(exc)
    try:
        from MACS3.Utilities.Execution import resolve_execution_policy, process_backend_reason
        policy = resolve_execution_policy()
        execution = {'default_policy': policy.describe(),
                     'auto_fallback_reason': process_backend_reason('auto'),
                     'limits': policy.availability.limits}
    except ImportError:
        execution = {'default_policy': 'legacy fixed policy; no per-invocation CPU API'}
    resource_files = {}
    for name in ('/sys/kernel/mm/transparent_hugepage/enabled',
                 '/sys/kernel/mm/transparent_hugepage/hpage_pmd_size',
                 '/sys/fs/cgroup/cpu.max', '/sys/fs/cgroup/memory.max',
                 '/proc/self/cgroup'):
        try:
            resource_files[name] = Path(name).read_text().strip()
        except OSError:
            pass
    return {'python': sys.version, 'executable': sys.executable,
            'macs_version': MACS_VERSION,
            'platform': platform.platform(), 'machine': platform.machine(),
            'cpu_count': os.cpu_count(), 'package': MACS3.__file__,
            'affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
            'address_space_limit': resource.getrlimit(resource.RLIMIT_AS),
            'resource_files': resource_files,
            'compiled_modules': modules, 'allocator_probe': handlers,
            'allocator_cleanup_notes': cleanup_notes,
            'compiled_module_sha256': {name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                                       for name, path in modules.items()},
            'execution': execution,
            'cli': str(Path(sysconfig.get_path('scripts')) / 'macs3'),
            'build': {key: sysconfig.get_config_var(key)
                      for key in ('CC', 'CFLAGS', 'LDSHARED')},
            'environment': {key: value for key, value in os.environ.items()
                            if key in ('CC', 'CFLAGS', 'CPPFLAGS', 'LDFLAGS', 'PYTHONPATH',
                                       'MACS3_NUM_THREADS', 'OMP_NUM_THREADS',
                                       'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'LD_PRELOAD')},
            'dependencies': sorted(f'{d.metadata["Name"]}=={d.version}'
                                   for d in importlib.metadata.distributions())}


def cache(spec):
    import numpy as np
    from unittest.mock import patch
    from MACS3.Signal.Prob import (poisson_cdf_Q_log10_prepare as prepare,
                                   poisson_cdf_Q_log10_many as many, poisson_cdf)
    original = np.empty
    allocations = []

    def allocate(*args, **kwargs):
        result = original(*args, **kwargs)
        allocations.append((weakref.ref(result), result.nbytes))
        return result

    # Stay inside the smallest prepared prefix so scoring cannot perform
    # uninstrumented growth between the allocation/lifetime observations.
    n = np.array([0, 7, 100, 500], dtype='u4')
    lam = np.array([0.5, 0.25, 100., 500.], dtype='f8')
    out = np.empty(len(n), dtype='f8')
    phases = []
    for maximum in (500, 65536, 2**22 - 2, 2**22 - 2, 500):
        start = time.perf_counter()
        with patch.object(np, 'empty', allocate):
            ready = prepare(maximum)
        elapsed = time.perf_counter() - start
        if not ready:
            raise RuntimeError('log-table preparation failed')
        # Repeated scoring jobs in one interpreter retain the current cache.
        for _ in range(spec.get('jobs', 3)):
            many(n, lam, out)
            expected = np.array([poisson_cdf(int(k), v, False, True) for k, v in zip(n, lam)])
            if not np.array_equal(out.view('u8'), expected.view('u8')):
                raise RuntimeError('scalar/vector probability mismatch')
        gc.collect()
        phases.append({'max_n': maximum, 'prepare_seconds': elapsed,
                       'live_table_bytes': sum(size for ref, size in allocations if ref() is not None),
                       'live_table_arrays': sum(ref() is not None for ref, _ in allocations),
                       'allocations_to_date': len(allocations),
                       'scores_hex': [float(x).hex() for x in out], **peak_rss()})
    instrumentation = 'weak references to NumPy empty allocations during public prepare only'
    if not allocations:
        # 170556d owns native realloc buffers rather than NumPy arrays. RSS
        # still measures growth, but this observer cannot count that payload.
        instrumentation = 'native cache payload unobserved; process high-water RSS only'
        for phase in phases:
            for field in ('live_table_bytes', 'live_table_arrays', 'allocations_to_date'):
                phase[field] = None
    caller_jobs = []
    try:
        from MACS3.Utilities import Execution
    except ImportError:
        caller_note = 'legacy lacks explicit serial policy; repeated scoring jobs only'
    else:
        from MACS3.Signal.CallPeakUnit import CallerFromAlignments
        from MACS3.Signal.FixWidthTrack import FWTrack
        caller_note = 'repeated independent serial callers in the warm-cache interpreter'
        for _ in range(spec.get('jobs', 3)):
            track = FWTrack(fw=50, threads=1)
            for chrom in (b'chr1', b'chr2'):
                for pos in [1000] * 40 + [9000]:
                    track.add_loc(chrom, pos, 0)
            track.finalize()
            caller = CallerFromAlignments(track, None, d=200, ctrl_d_s=[],
                                          ctrl_scaling_factor_s=[], lambda_bg=0.5,
                                          threads=1, process_backend='serial')
            try:
                peaks = caller.call_peaks(['q'], [3.0], min_length=100)
                rows = [(chrom.hex(), [tuple(row[k] for k in
                         ('start', 'end', 'summit', 'score', 'pileup', 'pscore', 'qscore', 'fc'))
                         for row in peaks.peaks[chrom]]) for chrom in sorted(peaks.peaks)]
                if not rows:
                    raise RuntimeError('Repeated caller produced no peaks')
            finally:
                caller.destroy()
                track.destroy()
            gc.collect()
            caller_jobs.append({'output_sha256': hashlib.sha256(json.dumps(rows).encode()).hexdigest(),
                                'live_table_bytes': sum(size for ref, size in allocations if ref() is not None),
                                **peak_rss()})
        if len({job['output_sha256'] for job in caller_jobs}) != 1:
            raise RuntimeError('Repeated caller output changed')
    return {'phases': phases, 'jobs_per_phase': spec.get('jobs', 3),
            'caller_jobs': caller_jobs, 'caller_note': caller_note,
            'instrumentation': instrumentation}


def main():
    operation, filename = sys.argv[1:]
    spec = json.loads(Path(filename).read_text())
    constrain(spec)
    if operation == 'metadata':
        result = metadata(spec)
    elif operation == 'cache':
        result = cache(spec)
    elif operation == 'cli':
        sys.argv = [spec['cli'], *spec['args']]
        try:
            runpy.run_path(spec['cli'], run_name='__main__')
        except SystemExit as exc:
            if exc.code not in (None, 0):
                raise
        result = peak_rss()
    else:
        raise ValueError(operation)
    dump(spec['result'], result)


if __name__ == '__main__':
    main()

"""Bound native table lifetime and allocation-failure probes in subprocesses."""
import subprocess
import sys

import pytest

from conftest import macs3_env


SCRIPT = r'''
import math
import sys
import threading
import weakref
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch

import numpy as np
from MACS3.Signal.Prob import (
    poisson_cdf, poisson_cdf_Q_log10_many as many,
    poisson_cdf_Q_log10_prepare as prepare,
    poisson_cdf_Q_log10_range as score_range,
)


def reference(n, lam):
    # The ordered, uncached scalar algorithm, including uint32 wrap and round.
    m = (int(n) + 1) % (2**32)
    total = 0.0
    for i in range(1, m + 1):
        total += math.log(i)
    ln_lam = math.log(lam)
    logx = m * ln_lam - total
    residue = logx
    while True:
        m += 1
        logy = logx + ln_lam - math.log(m)
        hi, lo = max(residue, logy), min(residue, logy)
        new = hi + math.log1p(math.exp(lo - hi))
        if abs(residue - new) < 1e-5:
            return round((new - lam) / math.log(10), 5)
        residue = new
        logx = logy


mode = sys.argv[1]
if mode == 'failure':
    allocation, warm = map(int, sys.argv[2:])
    if warm:
        assert prepare(500)
    n = np.array([0, 7, 500, 4096, 2**32 - 1], dtype='u4')
    lam = np.array([0.5, 0.25, 500., 100., 0.5], dtype='f8')
    expected = np.array([reference(k, v) for k, v in zip(n, lam)])
    out = np.full(len(n), 123., dtype='f8')
    real_empty = np.empty
    calls = 0
    allocated = []

    def failing_empty(*args, **kwargs):
        global calls
        calls += 1
        if calls == allocation:
            raise MemoryError('injected table allocation failure')
        result = real_empty(*args, **kwargs)
        allocated.append(weakref.ref(result))
        return result

    with patch.object(np, 'empty', failing_empty):
        assert not prepare(4096)
    assert calls == allocation
    assert all(ref() is None for ref in allocated), 'partial generation leaked'

    # Keep every later growth failing: neither scalar nor vector fallback may
    # depend on a successful retry, and the old cached prefix must stay valid.
    with patch.object(np, 'empty', side_effect=MemoryError('still failing')):
        for k, v, value in zip(n, lam, expected):
            assert poisson_cdf(int(k), v, False, True).hex() == float(value).hex()
        many(n, lam, out)
        assert list(out.view('u8')) == list(expected.view('u8'))
        out[:] = 123.
        score_range(n, lam, out, 1, len(n) - 1)
        assert list(out.view('u8')[1:-1]) == list(expected.view('u8')[1:-1])
        assert out[0] == out[-1] == 123.
        score_range(n, lam, out, 2, 2)
        assert out[0] == 123.
    assert prepare(4096), 'failed growth must allow retry'
    many(n, lam, out)
    np.testing.assert_array_equal(out.view('u8'), expected.view('u8'))

elif mode == 'lifetime':
    real_empty = np.empty
    allocated = []
    reenter = False

    def tracking_empty(*args, **kwargs):
        global reenter
        if reenter:
            reenter = False
            # Simulate another caller publishing while allocation releases
            # the GIL. The outer call must not replace this larger generation.
            assert prepare(16384)
        result = real_empty(*args, **kwargs)
        allocated.append(weakref.ref(result))
        return result

    with patch.object(np, 'empty', tracking_empty):
        assert prepare(500)
        assert prepare(2048)
        assert sum(ref() is not None for ref in allocated) == 2
        reenter = True
        assert prepare(8192)
    assert sum(ref() is not None for ref in allocated) == 2
    with patch.object(np, 'empty', side_effect=AssertionError('cache regressed')):
        assert prepare(16384)

elif mode == 'concurrent':
    # Fresh process: serial comparisons happen AFTER all growth. Otherwise
    # preparing the reference could silently remove the race from the fixture.
    assert prepare(500)
    barrier = threading.Barrier(5, timeout=15)
    counts = [1024, 4096, 16384, 65536, 262144, 1048576, 4194302]

    def worker(role):
        records = []
        for count in counts:
            n = np.full(16384 if role == 0 else 8,
                        500 if role == 0 else count, dtype='u4')
            lam = np.full(len(n), 500. if role == 0 else 100., dtype='f8')
            out = np.full(len(n), 123., dtype='f8')
            barrier.wait()
            if role == 1:
                records.append((count, 100., poisson_cdf(count, 100., False, True)))
            elif role == 2:
                many(n, lam, out)
                assert np.all(out == out[0])
                records.append((count, 100., out[0]))
            elif role == 3:
                assert prepare(count)
            else:
                score_range(n, lam, out, 1, len(n) - 1)
                assert out[0] == out[-1] == 123.
                assert np.all(out[1:-1] == out[1])
                records.append((int(n[1]), lam[1], out[1]))
        return records

    with ThreadPoolExecutor(max_workers=5) as pool:
        results = list(pool.map(worker, range(5)))
    expected = {}
    for records in results:
        for n, lam, value in records:
            key = n, lam
            if key not in expected:
                expected[key] = reference(n, lam)
            assert float(value).hex() == expected[key].hex(), (key, value, expected[key])
            assert float(value).hex() == poisson_cdf(n, lam, False, True).hex()

elif mode == 'caller':
    from MACS3.Signal.CallPeakUnit import CallerFromAlignments
    from MACS3.Signal.FixWidthTrack import FWTrack

    broad = bool(int(sys.argv[2]))

    def call():
        track = FWTrack(fw=50)
        for index in range(3):
            for pos in [1000 + index * 200] * (20 + index * 7) + [9000]:
                track.add_loc(('chr%d' % index).encode(), pos, 0)
        track.finalize()
        calc = CallerFromAlignments(
            track, None, d=200, ctrl_d_s=[], ctrl_scaling_factor_s=[],
            lambda_bg=0.5, process_backend='serial')
        try:
            if broad:
                peaks = calc.call_broadpeaks(['q'], [3.0], [1.0], min_length=100)
            else:
                peaks = calc.call_peaks(['q'], [3.0], min_length=100)
            keys = ['start', 'end', 'score', 'pileup', 'pscore', 'qscore', 'fc']
            keys += (['thickStart', 'thickEnd', 'blockNum', 'blockSizes', 'blockStarts']
                     if broad else ['summit'])
            rows = [(chrom, [tuple(row[key] for key in keys)
                             for row in peaks.peaks[chrom]])
                    for chrom in sorted(peaks.peaks)]
            assert rows and all(items for _, items in rows)
            return rows
        finally:
            calc.destroy()

    baseline = call()
    barrier = threading.Barrier(2, timeout=15)

    def grow():
        for count in [4096, 65536, 1048576, 4194302]:
            barrier.wait()
            assert prepare(count)

    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(grow)
        for _ in range(4):
            barrier.wait()
            assert call() == baseline
        future.result()

elif mode == 'boundary':
    # Small arrays exercise the table cap and uncached sum beyond it.
    n = np.array([0, 1022, 1023, 2047, 2**22 - 2, 2**22 - 1,
                  2**22, 2**32 - 1], dtype='u4')
    lam = np.full(len(n), 100., dtype='f8')
    expected = np.array([reference(k, v) for k, v in zip(n, lam)])
    out = np.empty(len(n), dtype='f8')
    # Unprepared range first, then prepared and batch paths.
    score_range(n, lam, out, 0, len(n))
    np.testing.assert_array_equal(out.view('u8'), expected.view('u8'))
    assert prepare(2**32 - 1)
    many(n, lam, out)
    np.testing.assert_array_equal(out.view('u8'), expected.view('u8'))
    for k, v, value in zip(n, lam, expected):
        assert poisson_cdf(int(k), v, False, True).hex() == float(value).hex()
    empty = np.empty(0, dtype='u4')
    empty_float = np.empty(0, dtype='f8')
    many(empty, empty_float, empty_float)
    score_range(empty, empty_float, empty_float, 0, 0)
else:
    raise AssertionError(mode)
print('ok')
'''


@pytest.mark.parametrize('mode,args', [
    ('failure', ('1', '0')),
    ('failure', ('2', '0')),
    ('failure', ('1', '1')),
    ('failure', ('2', '1')),
    ('concurrent', ()),
    ('lifetime', ()),
    ('boundary', ()),
    ('caller', ('0',)),
    ('caller', ('1',)),
])
def test_log_table_public_calls(tmp_path, mode, args):
    script = tmp_path / 'probe.py'
    script.write_text(SCRIPT)
    result = subprocess.run(
        [sys.executable, str(script), mode, *args], cwd=tmp_path,
        env=macs3_env(tmp_path), capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.strip() == 'ok'

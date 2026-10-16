"""Defaults contracts and public buffer/layout boundary regressions for #767."""

from dataclasses import FrozenInstanceError, fields, replace
import io
import subprocess
import sys

import pytest

from conftest import macs3_env
from MACS3.Utilities.Performance import DEFAULTS
from test_Parser_frag_eof import write_frag


def test_defaults_are_immutable_and_capacities_are_valid():
    with pytest.raises(FrozenInstanceError):
        DEFAULTS.cpu_limit = 1
    with pytest.raises(FrozenInstanceError):
        del DEFAULTS.cpu_limit
    assert not hasattr(DEFAULTS, "__dict__")
    assert all(type(getattr(DEFAULTS, field.name)) is int
               and getattr(DEFAULTS, field.name) > 0 for field in fields(DEFAULTS))
    for capacity in (DEFAULTS.frag_interner_initial_slots,
                     DEFAULTS.pscore_initial_slots, DEFAULTS.qscore_min_slots):
        assert capacity & (capacity - 1) == 0
    assert DEFAULTS.poisson_log_cache_min_entries <= DEFAULTS.poisson_log_cache_max_entries


@pytest.mark.parametrize("changes,diagnostic", [
    ({"cpu_limit": 0}, "cpu_limit"),
    ({"frag_tasks_per_thread": True}, "frag_tasks_per_thread"),
    ({"bgzf_batch_input_bytes": 1.5}, "bgzf_batch_input_bytes"),
    ({"pscore_initial_slots": 1000}, "power of two"),
    ({"poisson_log_cache_max_entries": 1 << 31}, "native-size|32-bit"),
    ({"poisson_log_cache_min_entries": DEFAULTS.poisson_log_cache_max_entries + 1},
     "minimum exceeds"),
])
def test_defaults_reject_invalid_policy(changes, diagnostic):
    # replace constructs another immutable record; it cannot change DEFAULTS.
    with pytest.raises(ValueError, match=diagnostic):
        replace(DEFAULTS, **changes)


FRAG_SCRIPT = r'''
from pathlib import Path
import sys
from MACS3.IO.Parser import FragParser
from MACS3.Signal.PairedEndTrack import PETrackII

path, reference, operation = sys.argv[1:]
rows = [line.split(b'\t') for line in Path(reference).read_bytes().splitlines()]
barcodes = {}
expected = {}
if operation == 'append':
    barcodes[b'SEED'] = 0
    expected[b'seed'] = [(1, 41, 1, 0)]
for chrom, left, right, barcode, count in rows:
    code = barcodes.setdefault(barcode, len(barcodes))
    expected.setdefault(chrom, []).append((int(left), int(right), int(count), code))
parser = FragParser(path)
if operation == 'build':
    track = parser.build_petrack()
else:
    track = PETrackII(buffer_size=2)
    track.add_loc(b'seed', 1, 41, b'SEED', 1)
    assert parser.append_petrack(track) is track
assert parser.n == len(rows) and parser.d == 40.0
track.finalize()
assert dict(track.barcode_dict) == barcodes
assert track.get_chr_names() == set(expected)
assert track.total == sum(row[2] for values in expected.values() for row in values)
for chrom, values in expected.items():
    values.sort()
    locations = track.get_locations_by_chr(chrom).tolist()
    assert locations == [row[:3] for row in values], chrom
    assert track.barcodes[chrom].tolist() == [row[3] for row in values], chrom
'''


def run_probe(tmp_path, script_text, *args):
    script = tmp_path / "probe.py"
    script.write_text(script_text)
    proc = subprocess.run([sys.executable, str(script), *map(str, args)],
                          cwd=tmp_path, env=macs3_env(tmp_path),
                          capture_output=True, text=True, timeout=30)
    assert proc.returncode == 0, proc.stdout + proc.stderr


@pytest.mark.parametrize("encoding", ["plain", "gzip", "bgzf"])
@pytest.mark.parametrize("operation", ["build", "append"])
def test_frag_layout_and_interner_growth(tmp_path, encoding, operation):
    # ~1 s per probe. Exceeds the threaded byte threshold, local chromosome
    # and run bounds, initial chromosome reserve, hash load and arena reserve.
    # Distinct keys straddle the inline-key bound; later rows repeat them.
    barcodes = [(str(i) + "B" * (length - len(str(i)))).encode()
                for i, length in enumerate([35, 36, 37, 80] * 384)]
    data = b"".join(
        b"\t".join([f"chr{i % 73:02}".encode(), str(i * 3).encode(),
                    str(i * 3 + 40).encode(), barcodes[i % len(barcodes)],
                    str(i % 5 + 1).encode()]) + b"\n"
        for i in range(6144))
    reference = tmp_path / "reference.tsv"
    reference.write_bytes(data)
    path = write_frag(tmp_path / "fragments", data, encoding)
    run_probe(tmp_path, FRAG_SCRIPT, path, reference, operation)


TRACK_SCRIPT = r'''
import sys
from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI, PETrackII

kind = sys.argv[1]
track = (FWTrack(fw=40) if kind == 'fw' else
         PETrackI() if kind == 'pe' else PETrackII())
sizes = {b'below': 65535, b'at': 65536, b'above': 65537}
for chrom, size in sizes.items():
    for i in reversed(range(size)):
        left = i // 2
        if kind == 'fw':
            track.add_loc(chrom, left, 0)
        elif kind == 'pe':
            track.add_loc(chrom, left, left + 40)
        else:
            track.add_loc(chrom, left, left + 40, str(left % 17).encode(), 2)
track.finalize()
assert track.total == sum(sizes.values()) * (2 if kind == 'frag' else 1)
for chrom, size in sizes.items():
    expected = [i // 2 for i in range(size)]
    locs = track.get_locations_by_chr(chrom)
    if kind == 'fw':
        assert locs[0].tolist() == expected and not len(locs[1])
    else:
        assert locs.tolist() == [(left, left + 40) + ((2,) if kind == 'frag' else ())
                                 for left in expected]
        if kind == 'frag':
            codes = track.barcode_dict
            assert track.barcodes[chrom].tolist() == [codes[str(left % 17).encode()]
                                                      for left in expected]
if kind != 'frag':
    track.filter_dup(1)
    assert track.total == sum((size + 1) // 2 for size in sizes.values())
    for chrom, size in sizes.items():
        locs = track.get_locations_by_chr(chrom)
        expected = list(range((size + 1) // 2))
        if kind == 'fw':
            assert locs[0].tolist() == expected
        else:
            assert locs.tolist() == [(left, left + 40) for left in expected]
'''


@pytest.mark.parametrize("kind", ["fw", "pe", "frag"])
def test_tracks_at_parallel_threshold(tmp_path, kind):
    # ~1 s. Two chromosomes meet the threshold, allowing threaded work.
    run_probe(tmp_path, TRACK_SCRIPT, kind)


@pytest.mark.parametrize("method", ["write_to_bed", "write_to_summit_bed",
                                   "write_to_narrowPeak"])
def test_peak_writer_growth_preserves_utf8_and_formatting(method):
    from MACS3.IO.PeakIO import PeakIO

    peaks = PeakIO()
    chrom = "chrΩ"
    prefix = "試験" * 80
    expected = []
    for i in range(1500):
        left = i * 50
        peaks.add(chrom.encode(), left, left + 40, summit=left + 20,
                  peak_score=2.5, pileup=7.0, pscore=3.25,
                  fold_change=4.5, qscore=2.5)
        name = f"{prefix}{i + 1}"
        if method == "write_to_bed":
            line = f"{chrom}\t{left}\t{left + 40}\t{name}\t2.5\n"
        elif method == "write_to_summit_bed":
            line = f"{chrom}\t{left + 20}\t{left + 21}\t{name}\t2.5\n"
        else:
            line = (f"{chrom}\t{left}\t{left + 40}\t{name}\t25\t.\t"
                    "4.5\t3.25\t2.5\t20\n")
        expected.append(line)
    output = io.StringIO()
    getattr(peaks, method)(output, name_prefix=prefix.encode(), trackline=False)
    assert output.getvalue() == "".join(expected)
    assert len(output.getvalue().encode()) > 65536

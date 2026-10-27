"""Public-API and CLI regressions for #800, part of #757."""

import gzip
from pathlib import Path
import random
import subprocess
import sys

import MACS3
import pytest

from conftest import macs3_env, no_core_dumps


HEADER = "@SQ\tSN:chr1\tLN:1000000\n"


def record(name="rev", flag=16, cigar="50M", pos=5000):
    return (f"{name}\t{flag}\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t"
            + "A" * 50 + "\t" + "I" * 50 + "\n")


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gzip"])
@pytest.mark.parametrize("flag", [16, 83], ids=["single", "first-mate"])
def test_sam_reverse_append(tmp_path, gz, flag):
    path = tmp_path / ("reads.sam.gz" if gz else "reads.sam")
    text = HEADER + record(flag=flag, cigar="5H3S10M2I3D100N4=2X1P29M")
    with (gzip.open if gz else open)(path, "wt") as stream:
        stream.write(text)
    run_probe(tmp_path, path, "append", '''
from MACS3.Signal.FixWidthTrack import FWTrack
track = FWTrack()
track.add_loc(b'chr1', 100, 0)
assert Parser.SAMParser(sys.argv[1]).append_fwtrack(track) is track
track.finalize()
plus, minus = track.get_locations_by_chr(b'chr1')
assert plus.tolist() == [100]
# POS - 1 + M + D + N + = + X + M; clips, insertion and pad add zero.
assert minus.tolist() == [5147]
''')


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gzip"])
@pytest.mark.parametrize("case", ["empty", "header", "blank", "truncated",
                                 "missing-cigar", "invalid-cigar"])
def test_sam_edge_inputs(tmp_path, gz, case):
    text = {"empty": "", "header": HEADER, "blank": HEADER + "\n",
            "truncated": HEADER + "r\t16\tchr1\t5000\n",
            "missing-cigar": HEADER + record(cigar="*"),
            "invalid-cigar": HEADER + record(cigar="garbage")}[case]
    path = tmp_path / ("reads.sam.gz" if gz else "reads.sam")
    with (gzip.open if gz else open)(path, "wt") as stream:
        stream.write(text)
    run_probe(tmp_path, path, case, '''
parser = Parser.SAMParser(sys.argv[1])
if sys.argv[2] == 'truncated':
    try:
        parser.build_fwtrack()
    except IndexError:
        pass
    else:
        raise AssertionError('Truncated SAM record must raise IndexError')
else:
    track = parser.build_fwtrack()
    track.finalize()
    if sys.argv[2] in ('missing-cigar', 'invalid-cigar'):
        # Preserve existing permissive behavior; no reference operations.
        plus, minus = track.get_locations_by_chr(b'chr1')
        assert plus.tolist() == [] and minus.tolist() == [4999]
    else:
        assert track.get_chr_names() == set()
''')


def run_probe(tmp_path, path, case, source):
    script = tmp_path / "probe.py"
    script.write_text('''
from pathlib import Path
import sys
from MACS3.IO import Parser
assert Path(Parser.__file__).resolve().is_relative_to(Path(sys.argv[3]))
assert Path(Parser.__file__).suffix != '.py'
''' + source)
    root = Path(MACS3.__file__).resolve().parent.parent
    result = subprocess.run(
        [sys.executable, str(script), str(path), case, str(root)],
        cwd=tmp_path, env=macs3_env(tmp_path), capture_output=True,
        text=True, timeout=30, preexec_fn=no_core_dumps)
    assert result.returncode == 0, result.stdout + result.stderr


def test_sam_filterdup_original_reproduction(tmp_path, run_macs3):
    rng = random.Random(0)
    text = HEADER + "".join(
        record(name=f"r{i}", flag=0, pos=rng.randint(1000, 900000))
        for i in range(2000))
    sam = tmp_path / "reads.sam"
    sam.write_text(text + record())
    output = tmp_path / "out.bed"
    result = run_macs3(["filterdup", "-i", sam, "-f", "SAM", "-g", "1e6",
                        "-o", output], timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    rows = output.read_text().splitlines()
    assert [row for row in rows if row.endswith("\t-")] == [
        "chr1\t4999\t5049\t.\t.\t-"]
    sam.write_text(text)
    result = run_macs3(["filterdup", "-i", sam, "-f", "SAM", "-g", "1e6",
                        "-o", output], timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert output.read_text().splitlines() == [
        row for row in rows if row.endswith("\t+")]

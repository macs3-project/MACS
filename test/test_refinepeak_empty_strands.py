"""Public API/CLI regressions for strand-empty chromosomes (#812, #757)."""

import subprocess
import sys

import pytest

from conftest import macs3_env, no_core_dumps


REGION_PROBE = '''
from pathlib import Path
import sys
from MACS3.Signal import FixWidthTrack
from MACS3.IO.PeakIO import PeakIO

assert Path(FixWidthTrack.__file__).resolve() == Path(sys.argv[1]).resolve()
assert Path(FixWidthTrack.__file__).suffix != '.py'

def check(track, chromosomes):
    peaks = PeakIO()
    expected = []
    for chrom in sorted(chromosomes):
        plus, minus = track.get_locations_by_chr(chrom)
        for start, end in ((100, 200), (500, 600), (5000, 5100)):
            name = b'peak_' + str(start).encode()
            peaks.add(chrom, start, end, name=name)
            expected.append((chrom,
                             [p for p in plus.tolist() if start-50 <= p <= end+50],
                             [p for p in minus.tolist() if start-50 <= p <= end+50],
                             start-50, end+50, name, 50, 2.5))
    def collect(chrom, plus, minus, start, end, name, window_size, cutoff):
        assert plus.dtype.name == minus.dtype.name == 'int32'
        return (chrom, plus.tolist(), minus.tolist(), start, end, name,
                window_size, cutoff)
    result = track.compute_region_tags_from_peaks(
        peaks, collect, window_size=50, cutoff=2.5)
    assert result == expected, (chromosomes, result, expected)
    assert track.compute_region_tags_from_peaks(PeakIO(), collect) == []

for strand in (0, 1):
    for preceding in (False, True):
        for positions in ([150], [50, 100, 150, 250, 450, 550, 650, 900]):
            track = FixWidthTrack.FWTrack(buffer_size=2)
            chromosomes = [b'chrB', b'chrC']
            if preceding:
                chromosomes.append(b'chrA')
                for side in (0, 1):
                    for pos in (50, 100, 150, 250, 450, 550, 650, 900):
                        track.add_loc(b'chrA', pos, side)
            # Alternate which strand is empty across adjacent chromosomes.
            for pos in positions:
                track.add_loc(b'chrB', pos, strand)
                track.add_loc(b'chrC', pos, 1-strand)
            track.finalize()
            check(track, chromosomes)

# Sampling retains the chromosome keys even when both strands become empty.
for percent in (0.0, 0.25):
    track = FixWidthTrack.FWTrack(buffer_size=2)
    for side in (0, 1):
        for pos in range(50, 950, 10):
            track.add_loc(b'chrA', pos, side)
        track.add_loc(b'chrB', 150, side)
    track.finalize()
    track.sample_percent(percent, seed=17)
    assert all(len(x) == 0 for x in track.get_locations_by_chr(b'chrB'))
    check(track, [b'chrA', b'chrB'])

track = FixWidthTrack.FWTrack()
track.finalize()
check(track, [])
'''


def test_region_tags_empty_strands(tmp_path):
    from MACS3.Signal import FixWidthTrack
    script = tmp_path / 'region_probe.py'
    script.write_text(REGION_PROBE)
    result = subprocess.run(
        [sys.executable, str(script), FixWidthTrack.__file__], cwd=tmp_path,
        env=macs3_env(tmp_path), capture_output=True, text=True,
        timeout=30, preexec_fn=no_core_dumps)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('preceding', [False, True])
def test_refinepeak_single_strand_chromosome(tmp_path, run_macs3,
                                            strand, preceding):
    # Original reproduction and its mirror, also with chrB first (#812).
    reads = []
    peaks = []
    if preceding:
        for i in range(200):
            start = 1000 + 20*i
            reads.extend([('chrA', start, '+'), ('chrA', start+100, '-')])
        peaks.append(('chrA', 1500, 1800))
    reads.extend(('chrB', 1000 + 20*i, strand) for i in range(100))
    peaks.extend([('chrB', 1200, 1500), ('chrB', 2200, 2500)])
    bed = tmp_path / 'reads.bed'
    bed.write_text(''.join(f'{c}\t{s}\t{s+1}\tr\t0\t{st}\n'
                           for c, s, st in reads))
    narrow = tmp_path / 'peaks.narrowPeak'
    narrow.write_text(''.join(
        f'{c}\t{s}\t{e}\tpeak_{i}\t100\t.\t5.0\t10.0\t8.0\t150\n'
        for i, (c, s, e) in enumerate(peaks)))
    output = tmp_path / 'refined.bed'
    result = run_macs3(['refinepeak', '-b', narrow, '-i', bed,
                       '-f', 'BED', '-o', output], timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    rows = [line.split('\t') for line in output.read_text().splitlines()]
    assert len(rows) == len(peaks)
    # Independent tag collection checks scoring and valid both-strand output.
    from MACS3.Commands.refinepeak_cmd import find_summit
    expected = []
    for i, (chrom, start, end) in enumerate(peaks):
        plus = [s for c, s, st in reads
                if c == chrom and st == '+' and start-200 <= s <= end+200]
        minus = [s+1 for c, s, st in reads
                 if c == chrom and st == '-' and start-200 <= s+1 <= end+200]
        row = find_summit(chrom.encode(), plus, minus, start-200, end+200,
                          name=f'peak_{i}'.encode(), window_size=200)
        expected.append([row[0].decode(), str(row[1]), str(row[2]),
                         row[3].decode(), f'{row[4]:.2f}'])
    assert rows == expected
    assert all(row[3].endswith('_F') and float(row[4]) <= 0
               for row in rows if row[0] == 'chrB')

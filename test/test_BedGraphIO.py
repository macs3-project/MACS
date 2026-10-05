import gzip
import sys
import types

import numpy as np
import pytest

from MACS3.IO.BedGraphIO import bedGraphIO
from MACS3.Signal.BedGraph import bedGraphTrackI


# Provide lightweight Cython stubs so the module imports without compiled extensions.
cython_stub = sys.modules.get("cython")
if cython_stub is None:
    cython_stub = types.ModuleType("cython")
    sys.modules["cython"] = cython_stub


def _identity_decorator(*args, **kwargs):
    def decorate(func):
        return func
    if args and callable(args[0]) and len(args) == 1 and not kwargs:
        return args[0]
    return decorate


cython_stub.cfunc = getattr(cython_stub, "cfunc", _identity_decorator)
cython_stub.ccall = getattr(cython_stub, "ccall", _identity_decorator)
cython_stub.cclass = getattr(cython_stub, "cclass", _identity_decorator)
cython_stub.locals = getattr(cython_stub, "locals", _identity_decorator)
cython_stub.inline = getattr(cython_stub, "inline", _identity_decorator)
cython_stub.returns = getattr(cython_stub, "returns", _identity_decorator)
cython_stub.declare = getattr(cython_stub, "declare", lambda *a, **k: None)
cython_stub.short = getattr(cython_stub, "short", int)
cython_stub.float = getattr(cython_stub, "float", float)
cython_stub.double = getattr(cython_stub, "double", float)
cython_stub.int = getattr(cython_stub, "int", int)
cython_stub.long = getattr(cython_stub, "long", int)
cython_stub.bint = getattr(cython_stub, "bint", bool)

cimports_mod = sys.modules.get("cython.cimports")
if cimports_mod is None:
    cimports_mod = types.ModuleType("cython.cimports")
    sys.modules["cython.cimports"] = cimports_mod
cython_stub.cimports = cimports_mod

cpython_mod = sys.modules.get("cython.cimports.cpython")
if cpython_mod is None:
    cpython_mod = types.ModuleType("cython.cimports.cpython")
    cpython_mod.bool = bool
    sys.modules["cython.cimports.cpython"] = cpython_mod
cimports_mod.cpython = cpython_mod

numpy_mod = sys.modules.get("cython.cimports.numpy")
if numpy_mod is None:
    numpy_mod = types.ModuleType("cython.cimports.numpy")
    numpy_mod.ndarray = np.ndarray
    sys.modules["cython.cimports.numpy"] = numpy_mod
cimports_mod.numpy = numpy_mod


def to_python_arrays(track, chrom):
    data = track.get_data_by_chr(chrom)
    if not data:
        return [], []
    positions, values = data
    return list(positions), list(values)


def test_read_bedgraph_skips_headers_and_preserves_values(tmp_path):
    bedgraph_content = b"""track type=bedGraph name=\"demo\"\n# comment line\nbrowse position chr1\nchr1 10 20 1.50\nchr1 20 30 0.40\nchr2 0 5 -2.0\n"""
    input_path = tmp_path / "input.bdg"
    input_path.write_bytes(bedgraph_content)

    reader = bedGraphIO(str(input_path))
    result_track = reader.read_bedGraph(baseline_value=-1.0)

    if not hasattr(result_track, "get_chr_names"):
        pytest.skip("bedGraphTrackI implementation does not expose get_chr_names")

    assert result_track is reader.data
    assert result_track.get_chr_names() == {b"chr1", b"chr2"}

    if not hasattr(result_track, "get_data_by_chr"):
        pytest.skip("bedGraphTrackI implementation does not expose get_data_by_chr")

    chr1_pos, chr1_val = to_python_arrays(result_track, b"chr1")
    chr2_pos, chr2_val = to_python_arrays(result_track, b"chr2")

    assert chr1_pos == [10, 20, 30]
    assert chr1_val[0] == pytest.approx(-1.0)
    assert chr1_val[1] == pytest.approx(1.50)
    assert chr1_val[2] == pytest.approx(0.40)

    assert chr2_pos == [5]
    assert chr2_val == pytest.approx([-2.0])


def test_write_bedgraph_emits_trackline_and_sorted_regions(tmp_path):
    data = bedGraphTrackI()
    data.add_loc(b"chr2", 0, 4, 5.0)
    data.add_loc(b"chr1", 0, 3, 1.23456)
    data.add_loc(b"chr1", 3, 7, 0.789)

    output_path = tmp_path / "output.bdg"
    writer = bedGraphIO(str(output_path), data=data)
    writer.write_bedGraph(name='quot"ed', description='desc "here"', trackline=True)

    lines = output_path.read_text().splitlines()

    assert lines[0] == 'track type=bedGraph name="quot\\"ed" description="desc \\"here\\"" visibility=2 alwaysZero=on'
    assert lines[1] == "chr1\t0\t3\t1.23456"
    assert lines[2] == "chr1\t3\t7\t0.78900"
    assert lines[3] == "chr2\t0\t4\t5.00000"


# ------------------------------------
# Helpers for the tests below
# ------------------------------------

def _f32(x):
    """``x`` after a round trip through a C float (bedGraph values are float32)."""
    return float(np.float32(x))


def _read(tmp_path, text, baseline=0.0, name="in.bdg"):
    path = tmp_path / name
    path.write_text(text)
    return bedGraphIO(str(path)).read_bedGraph(baseline_value=baseline)


def _write(tmp_path, track, fname="out.bdg", **kwargs):
    path = tmp_path / fname
    bedGraphIO(str(path), data=track).write_bedGraph(**kwargs)
    return path.read_text()


def _roundtrip(tmp_path, text, baseline=0.0):
    """Read ``text`` as a bedGraph and write it back without a trackline."""
    return _write(tmp_path, _read(tmp_path, text, baseline), trackline=False)


EMPTY_TRACKLINE = ('track type=bedGraph name="" description="" '
                   'visibility=2 alwaysZero=on\n')


# ------------------------------------
# bedGraphIO.__init__
# ------------------------------------

def test_init_defaults(tmp_path):
    bio = bedGraphIO(str(tmp_path / "x.bdg"))
    assert bio.bedGraph_filename == str(tmp_path / "x.bdg")
    assert isinstance(bio.data, bedGraphTrackI)
    assert bio.data.get_chr_names() == set()


def test_init_uses_given_track(tmp_path):
    track = bedGraphTrackI()
    assert bedGraphIO(str(tmp_path / "x.bdg"), data=track).data is track


def test_init_rejects_other_data(tmp_path):
    with pytest.raises(AssertionError):
        bedGraphIO(str(tmp_path / "x.bdg"), data={"chr1": []})


# ------------------------------------
# bedGraphIO.read_bedGraph
# ------------------------------------

def test_read_contiguous_round_trip_exact(tmp_path):
    text = "chr1\t0\t10\t1.5\nchr1\t10\t25\t2.25\nchr2\t0\t5\t-3\n"
    assert _roundtrip(tmp_path, text) == (
        "chr1\t0\t10\t1.50000\nchr1\t10\t25\t2.25000\nchr2\t0\t5\t-3.00000\n")


def test_read_track_content(tmp_path):
    track = _read(tmp_path, "chr1\t0\t10\t0.1\nchr1\t10\t30\t2\n")
    p, v = track.get_data_by_chr(b"chr1")
    assert list(p) == [10, 30]
    assert list(v) == [_f32(0.1), 2.0]
    assert track.get_chr_names() == {b"chr1"}
    assert track.get_data_by_chr(b"chr2") == []


def test_read_space_separated_columns(tmp_path):
    assert _roundtrip(tmp_path, "chr1 0 10 1.5\nchr1  10 20\t2\n") == \
        "chr1\t0\t10\t1.50000\nchr1\t10\t20\t2.00000\n"


@pytest.mark.parametrize("header", [
    "track type=bedGraph name=x",
    "#comment",
    "# comment with space",
    "browser position chr1:1-100",
    "browse hide all",
])
def test_read_skips_header_lines_anywhere(tmp_path, header):
    text = (header + "\nchr1\t0\t10\t1\n" + header + "\nchr1\t10\t20\t2\n")
    assert _roundtrip(tmp_path, text) == \
        "chr1\t0\t10\t1.00000\nchr1\t10\t20\t2.00000\n"


@pytest.mark.parametrize("baseline,text", [
    (0.0, "0.00000"),
    (-1.5, "-1.50000"),
    (2.5, "2.50000"),
])
def test_read_leading_gap_gets_baseline(tmp_path, baseline, text):
    assert _roundtrip(tmp_path, "chr1\t10\t20\t1\n", baseline=baseline) == \
        "chr1\t0\t10\t%s\nchr1\t10\t20\t1.00000\n" % text


def test_read_merges_equal_adjacent_values(tmp_path):
    text = ("chr1\t0\t10\t1\nchr1\t10\t20\t1\nchr1\t20\t30\t2\n"
            "chr1\t30\t40\t1\n")
    assert _roundtrip(tmp_path, text) == (
        "chr1\t0\t20\t1.00000\nchr1\t20\t30\t2.00000\nchr1\t30\t40\t1.00000\n")


def test_read_unsorted_input_is_stored_in_file_order(tmp_path):
    # Unsorted input is outside the documented contract: add_loc says "The
    # caller is responsible for providing non-overlapping, sorted regions"
    # and the bdgcmp/bdgpeakcall docs say regions on a chromosome should be
    # continuous. read_bedGraph neither sorts nor rejects, so each end is
    # appended in file order (writing this track back gives start > end).
    # By hand: (20,30,2) -> leading block [0,20)=0 then end 30 value 2;
    # (0,10,1) -> end 10 value 1; (10,20,3) -> end 20 value 3.
    shuffled = "chr1\t20\t30\t2\nchr1\t0\t10\t1\nchr1\t10\t20\t3\n"
    track = _read(tmp_path, shuffled)
    assert to_python_arrays(track, b"chr1") == ([20, 30, 10, 20],
                                                [0.0, 2.0, 1.0, 3.0])


def test_read_many_chromosomes_written_in_bytes_order(tmp_path):
    text = "chr2\t0\t5\t2\nchr10\t0\t5\t10\nchr1\t0\t5\t1\nchrX\t0\t5\t3\n"
    track = _read(tmp_path, text)
    assert track.get_chr_names() == {b"chr1", b"chr2", b"chr10", b"chrX"}
    out = _write(tmp_path, track, trackline=False)
    assert out == ("chr1\t0\t5\t1.00000\nchr10\t0\t5\t10.00000\n"
                   "chr2\t0\t5\t2.00000\nchrX\t0\t5\t3.00000\n")


@pytest.mark.parametrize("value,text", [
    ("1e-3", "0.00100"),
    ("-0.5", "-0.50000"),
    ("7", "7.00000"),
    ("+2.5", "2.50000"),
    ("0.1", "0.10000"),
])
def test_read_value_parsing(tmp_path, value, text):
    assert _roundtrip(tmp_path, "chr1\t0\t10\t%s\n" % value) == \
        "chr1\t0\t10\t%s\n" % text


def test_read_negative_start_clipped_and_empty_interval_skipped(tmp_path):
    # add_loc clips start < 0 to 0 and ignores intervals ending at <= 0
    text = "chr1\t0\t0\t9\nchr1\t-5\t10\t1\n"
    assert _roundtrip(tmp_path, text) == "chr1\t0\t10\t1.00000\n"


def test_read_int32_max_end(tmp_path):
    assert _roundtrip(tmp_path, "chr1\t0\t2147483647\t1\n") == \
        "chr1\t0\t2147483647\t1.00000\n"


def test_read_float32_max_value(tmp_path):
    fmax = float(np.finfo(np.float32).max)
    out = _roundtrip(tmp_path, "chr1\t0\t10\t%r\n" % fmax)
    assert out == "chr1\t0\t10\t%.5f\n" % fmax


def test_read_too_few_columns_raises(tmp_path):
    with pytest.raises(IndexError):
        _read(tmp_path, "chr1\t0\t10\n")


def test_read_gzipped_bedgraph_is_not_decompressed(tmp_path):
    """Pins the current output.

    read_bedGraph opens the file with open(..., "rb") and has no gzip
    detection; gzip input is documented only for callpeak's tag files, not
    for bedGraph inputs, so the compressed bytes are parsed as text. With
    a fixed gzip header (mtime=0, no file name) the whole file is one
    whitespace-free token, so the column lookup raises IndexError. No
    specification says what reading compressed bytes should give.
    """
    path = tmp_path / "x.bdg.gz"
    with open(path, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="",
                           mtime=0) as fh:
            fh.write(b"chr1\t0\t10\t1.5\nchr1\t10\t20\t2\n")
    assert len(path.read_bytes().split()) == 1
    with pytest.raises(IndexError, match="list index out of range"):
        bedGraphIO(str(path)).read_bedGraph()


def test_read_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        bedGraphIO(str(tmp_path / "missing.bdg")).read_bedGraph()


def test_read_empty_file(tmp_path):
    track = _read(tmp_path, "")
    assert track.get_chr_names() == set()
    assert _write(tmp_path, track, trackline=False) == ""


# ------------------------------------
# bedGraphIO.write_bedGraph
# ------------------------------------

def _track(*rows, baseline=0.0):
    t = bedGraphTrackI(baseline_value=baseline)
    for chrom, s, e, v in rows:
        t.add_loc(chrom, s, e, v)
    return t


def test_write_default_trackline(tmp_path):
    out = _write(tmp_path, _track((b"chr1", 0, 10, 1.0)))
    assert out == EMPTY_TRACKLINE + "chr1\t0\t10\t1.00000\n"


def test_write_trackline_name_description(tmp_path):
    out = _write(tmp_path, _track((b"chr1", 0, 10, 1.0)), name="treat",
                 description='a "b" c')
    assert out.splitlines()[0] == ('track type=bedGraph name="treat" '
                                   'description="a \\"b\\" c" visibility=2 '
                                   'alwaysZero=on')


def test_write_without_trackline(tmp_path):
    out = _write(tmp_path, _track((b"chr1", 0, 10, 1.0)), name="x",
                 trackline=False)
    assert out == "chr1\t0\t10\t1.00000\n"


def test_write_empty_track(tmp_path):
    assert _write(tmp_path, bedGraphTrackI(), trackline=False) == ""
    assert _write(tmp_path, bedGraphTrackI(), fname="e2.bdg") == EMPTY_TRACKLINE


@pytest.mark.parametrize("value,text", [
    pytest.param(1.0 / 3.0, "0.33333", id="third"),
    pytest.param(1e-6, "0.00000", id="tiny"),
    pytest.param(5e-6, "0.00000", id="float32_below_half"),  # 4.99999987e-06
    pytest.param(-2.5, "-2.50000", id="negative"),
    pytest.param(123456.789, "123456.78906", id="float32_large"),  # .7890625
    pytest.param(0.0, "0.00000", id="zero"),
])
def test_write_value_format(tmp_path, value, text):
    out = _write(tmp_path, _track((b"chr1", 0, 10, value)), trackline=False)
    assert out == "chr1\t0\t10\t%s\n" % text


def test_write_leading_gap_uses_track_baseline(tmp_path):
    out = _write(tmp_path, _track((b"chr1", 5, 10, 1.0), baseline=0.5),
                 trackline=False)
    assert out == "chr1\t0\t5\t0.50000\nchr1\t5\t10\t1.00000\n"


def test_write_does_not_merge_unmerged_track(tmp_path):
    t = bedGraphTrackI()
    t.add_loc_wo_merge(b"chr1", 0, 10, 1.0)
    t.add_loc_wo_merge(b"chr1", 10, 20, 1.0)
    assert _write(tmp_path, t, trackline=False) == \
        "chr1\t0\t10\t1.00000\nchr1\t10\t20\t1.00000\n"


def test_write_overwrites_existing_file(tmp_path):
    path = tmp_path / "out.bdg"
    path.write_text("old content\n" * 5)
    bedGraphIO(str(path), data=_track((b"chr1", 0, 3, 2.0))).write_bedGraph(
        trackline=False)
    assert path.read_text() == "chr1\t0\t3\t2.00000\n"


def test_write_read_write_is_stable(tmp_path):
    t = _track((b"chr2", 0, 7, 0.25), (b"chr1", 3, 9, 1.75), (b"chr1", 9, 12, 4.0))
    first = _write(tmp_path, t, fname="a.bdg")
    again = _write(tmp_path, bedGraphIO(str(tmp_path / "a.bdg")).read_bedGraph(),
                   fname="b.bdg")
    assert first == again == EMPTY_TRACKLINE + (
        "chr1\t0\t3\t0.00000\nchr1\t3\t9\t1.75000\nchr1\t9\t12\t4.00000\n"
        "chr2\t0\t7\t0.25000\n")

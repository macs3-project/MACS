import sys
import types

import numpy as np
import pytest

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

from MACS3.IO.BedGraphIO import bedGraphIO
from MACS3.Signal.BedGraph import bedGraphTrackI


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

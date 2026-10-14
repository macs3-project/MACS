"""Native allocation failures through PeakModel.build(), isolated from pytest.

The session build takes roughly 15 seconds; each small child run takes about
one second. Run explicitly with --run-slow. Cython and a C compiler are needed.
"""

import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from conftest import macs3_env


pytestmark = pytest.mark.slow
TEST_DIR = Path(__file__).resolve().parent
SOURCE = TEST_DIR.parent / "MACS3" / "Signal" / "PeakModel.py"
ALLOC_IMPORT = "from cython.cimports.libc.stdlib import malloc, realloc, free"
TEST_IMPORT = """from cython.cimports.libc.stdlib import malloc
from cython.cimports.peakmodel_alloc import (
    pm_realloc as realloc, pm_free as free, pm_reset, pm_stats)
"""
CONTROLS = """

def allocation_reset(target, worker_only=False):
    pm_reset(target, worker_only)

def allocation_stats():
    stats: cython.size_t[8]
    pm_stats(cython.address(stats[0]))
    return (stats[0], stats[1], stats[2], stats[3],
            stats[4], stats[5], stats[6], stats[7])
"""
BUILD = """
from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

setup(name="peakmodel-allocation-test",
      ext_modules=cythonize([
          Extension("peakmodel_fault", ["peakmodel_fault.py"],
                    include_dirs=[numpy.get_include(), "."],
                    extra_compile_args=["-O2", "-pthread"],
                    extra_link_args=["-pthread"])
      ], compiler_directives={"language_level": 3}))
"""


@pytest.fixture(scope="session")
def allocation_module(tmp_path_factory):
    if importlib.util.find_spec("Cython") is None:
        pytest.skip("native fault injection requires Cython and a C compiler")
    directory = tmp_path_factory.mktemp("peakmodel-allocation")
    source = SOURCE.read_text()
    assert source.count(ALLOC_IMPORT) == 1
    (directory / "peakmodel_fault.py").write_text(
        source.replace(ALLOC_IMPORT, TEST_IMPORT) + CONTROLS)
    for name in ("peakmodel_alloc.pxd", "peakmodel_alloc.h"):
        shutil.copyfile(TEST_DIR / name, directory / name)
    (directory / "setup.py").write_text(BUILD)
    proc = subprocess.run([sys.executable, "setup.py", "build_ext", "--inplace"],
                          cwd=directory, capture_output=True, text=True,
                          timeout=180)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return directory


def run_model(directory, kind, target=0, worker_only=False,
              module="peakmodel_fault"):
    proc = subprocess.run(
        [sys.executable, str(TEST_DIR / "peakmodel_alloc_runner.py"),
         str(directory), module, kind, str(target), str(int(worker_only))],
        env=macs3_env(), capture_output=True, text=True, timeout=30)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return json.loads(proc.stdout)


@pytest.mark.parametrize("kind,target,size,growth", [
    ("last", 1, 64 * 4, False),
    ("between", 1, 64 * 4, False),
    ("ties", 2, 128 * 4, True),
    ("last", 2, 256 * 4, False),
    ("last", 3, 256 * 4, False),
    ("between", 2, 256 * 4, False),
    ("between", 3, 256 * 4, False),
    ("growth_last", 4, 512 * 4, True),
    ("growth_last", 5, 512 * 4, True),
    ("growth_between", 4, 512 * 4, True),
    ("growth_between", 5, 512 * 4, True),
])
def test_build_allocation_failure(allocation_module, kind, target, size, growth):
    result = run_model(allocation_module, kind, target)
    assert result["stats"][0] == target
    assert result["stats"][4:6] == [size, int(growth)]
    assert result["recovered"] == "not_enough_pairs"


def test_build_worker_allocation_failure(allocation_module):
    result = run_model(allocation_module, "normal", 1, worker_only=True)
    assert result["recovered"] == "built"


@pytest.mark.parametrize("kind", ["normal", "normal_ties", "ties"])
def test_build_instrumented_matches_production(allocation_module, kind):
    expected = run_model(allocation_module, kind, module="MACS3.Signal.PeakModel")
    if expected["outcome"] == "unavailable":
        pytest.skip("PeakModel public attributes unavailable in current build")
    assert run_model(allocation_module, kind) == expected

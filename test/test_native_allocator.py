"""Allocator lifecycle through public tracks, in bounded fresh processes.

The Linux shim builds in <1 s; each growth/sort/sample probe takes ~1–3 s
and <150 MiB. Use --run-slow. No MACS3 native symbol is patched or called.
"""

import os
from pathlib import Path
import shlex
import subprocess
import sys
import sysconfig

import pytest

from test_native_streams import native_probe

pytestmark = pytest.mark.slow


@pytest.fixture(scope='session')
def sysfs_shim(tmp_path_factory):
    directory = tmp_path_factory.mktemp('native-sysfs')
    library = directory / 'sysfs.so'
    compiler = shlex.split(sysconfig.get_config_var('CC') or 'cc')
    result = subprocess.run([*compiler, '-shared', '-fPIC', '-O2',
                             str(Path(__file__).with_name('native_sysfs.c')),
                             '-o', str(library), '-ldl'],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    return library


@pytest.mark.skipif(sys.platform != 'linux', reason='Linux allocator discovery')
@pytest.mark.parametrize('mode', ['enabled', 'disabled', 'missing', 'invalid-size'])
def test_linux_track_allocator_lifecycle(tmp_path, sysfs_shim, mode):
    fixture = tmp_path / 'sysfs'
    fixture.mkdir()
    if mode != 'missing':
        (fixture / 'enabled').write_text('always [madvise] never' if mode != 'disabled'
                                         else 'always madvise [never]')
        (fixture / 'hpage_pmd_size').write_text('2097152' if mode != 'invalid-size' else '3')
    preload = str(sysfs_shim)
    if os.environ.get('LD_PRELOAD'):
        preload += os.pathsep + os.environ['LD_PRELOAD']
    expected = 'macs3_huge_pages' if mode == 'enabled' else 'default_allocator'
    results = native_probe(tmp_path, 'allocator', expected,
                           extra_env={'LD_PRELOAD': preload, 'MACS3_TEST_SYSFS': str(fixture)})
    assert len(results) == 6
    assert {row[2] for row in results[:3]} == {expected}


@pytest.mark.skipif(sys.platform != 'darwin', reason='macOS safe allocator fallback')
def test_macos_track_allocator_fallback(tmp_path):
    results = native_probe(tmp_path, 'allocator', 'default_allocator')
    assert len(results) == 6
    assert {row[2] for row in results} == {'default_allocator'}

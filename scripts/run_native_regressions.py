"""Run native compatibility tests with explicit parent installation guards.

Run from a source checkout. --installed requires that checkout to have no
in-place extensions, imports the installed package before pytest collection,
and verifies that every subprocess test derives that same installation.
"""

import argparse
import importlib
from pathlib import Path
import sys


SUITES = ['native_compatibility', 'native_streams', 'native_allocator', 'native_cli',
          'Execution_paths', 'Performance', 'Parser_imports',
          'Parser_bam_truncation', 'Parser_frag_eof',
          'chromosome_batching', 'CallPeakUnit_backend', 'CallPeakUnit_concurrency',
          'thread_lifecycle', 'pileup_lifecycle', 'caller_ownership', 'PeakModel_allocation']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--installed', action='store_true')
    parser.add_argument('--run-slow', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.installed:
        assert not list((root / 'MACS3').rglob('*.so'))
        assert not list((root / 'MACS3').rglob('*.pyd'))
    else:
        sys.path.insert(0, str(root))
    import MACS3
    package = Path(MACS3.__file__).resolve().parent
    assert package.is_relative_to(root) != args.installed, package
    for name in ('IO.Parser', 'IO.BAMStream', 'Signal.PairedEndTrack',
                 'Signal.CallPeakUnit', 'Signal.CallerStorage', 'Signal.PeakModel'):
        module = importlib.import_module('MACS3.' + name)
        path = Path(module.__file__).resolve()
        assert path.is_relative_to(package) and path.suffix in ('.so', '.pyd'), path
        print('Selected compiled module:', path, flush=True)
    import pytest
    tests = [str(root / 'test' / ('test_' + suite + '.py')) for suite in SUITES]
    tests.append(str(root / 'test' / 'test_callpeak_cmd.py') +
                 '::test_run_in_threaded_host_avoids_peak_pools')
    raise SystemExit(pytest.main(['-q', *tests, *(['--run-slow'] if args.run_slow else [])]))


if __name__ == '__main__':
    main()

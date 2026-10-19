"""Verify native build inputs and third-party notices in release artifacts.

Usage: python scripts/check_native_distribution.py dist/*.tar.gz dist/*.whl
Run against the same source revision that produced the artifacts.
"""

import argparse
from pathlib import Path
import tarfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
LICENSES = ['LICENSE', 'MACS3/IO/libdeflate/COPYING',
            'MACS3/Signal/NUMPY_LICENSE.txt', 'MACS3/fermi-lite/LICENSE.txt',
            'MACS3/fermi-lite/lib/COPYING']
NATIVE_INPUTS = ['MACS3/IO/bgzf_mt.c', 'MACS3/IO/bgzf_mt.h',
                 'MACS3/IO/BAMStream.py', 'MACS3/IO/BAMStream.pxd',
                 'MACS3/IO/cbgzf_mt.pxd', 'MACS3/IO/czlib.pxd',
                 'MACS3/Signal/chugepages.pxd', 'MACS3/Signal/cmalloc.pxd',
                 'MACS3/Signal/cPosValCalculation.c', 'MACS3/Signal/cPosValCalculation.h',
                 'MACS3/Signal/swalign.c', 'MACS3/Signal/swalign.h']
EXTENSIONS = ['IO/Parser', 'IO/BAMStream', 'Signal/PairedEndTrack',
              'Signal/CallPeakUnit', 'Signal/PeakModel', 'Signal/PileupV2']


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def check_sdist(path):
    with tarfile.open(path) as bundle:
        names = bundle.getnames()
        prefix = names[0].split('/')[0] + '/'
        relative = {name.removeprefix(prefix) for name in names}
        required = set(LICENSES + NATIVE_INPUTS)
        for directory in ('MACS3/IO/libdeflate', 'MACS3/fermi-lite'):
            required.update(str(p.relative_to(ROOT)) for p in (ROOT / directory).rglob('*')
                            if p.suffix in ('.c', '.h') and p.is_file())
        for name in sorted(required):
            require(name in relative, f'{path}: missing native input/license {name}')
            require(bundle.extractfile(prefix + name).read() == (ROOT / name).read_bytes(),
                    f'{path}: stale native input/license {name}')
        require(not any(name.startswith('test/') for name in relative), 'sdist unexpectedly includes tests')
        require(not any(name.endswith(('.so', '.pyd')) for name in relative), 'sdist contains in-place binaries')


def check_wheel(path):
    with zipfile.ZipFile(path) as bundle:
        names = bundle.namelist()
        for name in ('MACS3/IO/BAMStream.py', 'MACS3/IO/BAMStream.pxd'):
            require(name in names, f'{path}: missing internal stream source {name}')
            require(bundle.read(name) == (ROOT / name).read_bytes(),
                    f'{path}: stale internal stream source {name}')
        for module in EXTENSIONS:
            require(any(name.startswith('MACS3/' + module + '.') and name.endswith(('.so', '.pyd'))
                        for name in names), f'{path}: missing compiled {module}')
        for license in LICENSES:
            candidates = [name for name in names if '.dist-info/licenses/' in name
                          and name.split('.dist-info/licenses/', 1)[1] == license]
            require(len(candidates) == 1, f'{path}: missing/ambiguous binary notice {license}')
            require(bundle.read(candidates[0]) == (ROOT / license).read_bytes(),
                    f'{path}: stale binary notice {license}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('artifacts', type=Path, nargs='+')
    args = parser.parse_args()
    for path in args.artifacts:
        if path.name.endswith('.tar.gz'):
            check_sdist(path)
        elif path.suffix == '.whl':
            check_wheel(path)
        else:
            parser.error(f'unsupported artifact: {path}')
        print(f'Verified {path}')


if __name__ == '__main__':
    main()

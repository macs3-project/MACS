"""Public parser imports in fresh interpreters, including stream-first loading."""

import os
import subprocess
import sys

import pytest

from conftest import macs3_env


@pytest.mark.parametrize('first', ['MACS3.IO.Parser', 'MACS3.IO.BAMStream'])
def test_public_parser_imports_in_fresh_process(tmp_path, first):
    script = tmp_path / 'imports.py'
    script.write_text('''
import importlib
from pathlib import Path
import sys
import MACS3

importlib.import_module(sys.argv[1])
from MACS3.IO import Parser
from MACS3.IO.Parser import (
    guess_parser, StrandFormatError, GenericParser, BEDParser, BEDPEParser,
    ELANDResultParser, ELANDMultiParser, ELANDExportParser, SAMParser,
    BAMParser, BAMPEParser, BowtieParser, FragParser)

assert str(Path(MACS3.__file__).resolve().parent.parent) == sys.argv[2]
for value in (guess_parser, StrandFormatError, GenericParser, BEDParser,
              BEDPEParser, ELANDResultParser, ELANDMultiParser,
              ELANDExportParser, SAMParser, BAMParser, BAMPEParser,
              BowtieParser, FragParser):
    assert value.__module__ == 'MACS3.IO.Parser'
    assert getattr(Parser, value.__name__) is value
assert issubclass(BAMPEParser, BAMParser)
assert issubclass(FragParser, GenericParser)
''')
    env = macs3_env(tmp_path)
    result = subprocess.run([sys.executable, str(script), first,
                             env['PYTHONPATH'].split(os.pathsep)[0]],
                            cwd=tmp_path, env=env, capture_output=True,
                            text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr

#!/usr/bin/env python

"""Module Description: Test the memory-reporting logger in
MACS3/Utilities/Logger.py.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import importlib
import logging
import os
import re
import resource
import subprocess
import sys
import time
from collections import namedtuple

import pytest

import MACS3.Utilities.Logger as Logger
from MACS3.Utilities.Logger import (MemoryLogger)

# ------------------------------------
# helpers
# ------------------------------------

FAKE_MB = 1234
PREFIX = re.compile(r"^\[(\d+) MB\] ")


class ListHandler(logging.Handler):
    """Keep every record handled."""

    def __init__(self):
        super().__init__(level=logging.NOTSET)
        self.records = []

    def emit(self, record):
        self.records.append(record)


@pytest.fixture
def memlogger(request):
    """A MemoryLogger outside the logging manager with a ListHandler and
    no propagation, so the root logger and pytest never see it."""
    lg = MemoryLogger("macs3test." + request.node.name, logging.DEBUG)
    lg.propagate = False
    handler = ListHandler()
    lg.addHandler(handler)
    yield lg, handler.records
    lg.removeHandler(handler)


@pytest.fixture
def fake_memory(monkeypatch):
    """Make get_memory_usage return FAKE_MB and count its calls."""
    calls = []

    def fake():
        calls.append(1)
        return FAKE_MB
    monkeypatch.setattr(MemoryLogger, "get_memory_usage", staticmethod(fake))
    return calls


def run_python(code, timeout=60):
    """Run ``code`` in a fresh interpreter with the same MACS3 build."""
    import MACS3
    env = os.environ.copy()
    pkg_parent = os.path.dirname(os.path.dirname(
        os.path.abspath(MACS3.__file__)))
    env["PYTHONPATH"] = pkg_parent
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[var] = "1"
    return subprocess.run([sys.executable, "-c", code], env=env,
                          capture_output=True, text=True, timeout=timeout)


# ------------------------------------
# MemoryLogger
# ------------------------------------

def test_memorylogger_is_a_logging_logger():
    assert issubclass(MemoryLogger, logging.Logger)
    lg = MemoryLogger("macs3test.ctor")
    assert lg.name == "macs3test.ctor"
    assert lg.level == logging.NOTSET


@pytest.mark.parametrize("level", [logging.DEBUG, logging.WARNING, 5])
def test_memorylogger_constructor_level(level):
    assert MemoryLogger("macs3test.level", level).level == level


def test_import_sets_memorylogger_as_logger_class():
    assert logging.getLoggerClass() is MemoryLogger
    lg = logging.getLogger("macs3test.created_after_import")
    assert isinstance(lg, MemoryLogger)


@pytest.mark.parametrize("module", [
    "MACS3.Utilities.OptValidator",
    "MACS3.IO.Parser",
    "MACS3.Signal.CallPeakUnit",
    "MACS3.Signal.PairedEndTrack",
    "MACS3.Signal.HMMR_EM",
    "MACS3.Signal.HMMR_Signal_Processing",
    "MACS3.Commands.bdgpeakcall_cmd",
    "MACS3.Commands.bdgbroadcall_cmd",
    "MACS3.Commands.bdgdiff_cmd",
])
def test_macs3_module_loggers_are_memoryloggers(module):
    # these modules import MACS3.Utilities.Logger before getLogger
    mod = importlib.import_module(module)
    assert isinstance(mod.logger, MemoryLogger)
    assert mod.logger.name == module


def test_message_gets_memory_prefix(memlogger, fake_memory):
    lg, records = memlogger
    lg.info("hello %s", "world")
    assert len(records) == 1
    assert records[0].getMessage() == "[1234 MB] hello world"
    assert records[0].msg == "[1234 MB] hello %s"
    assert records[0].args == ("world",)
    assert fake_memory == [1]


def test_message_prefix_uses_real_memory_usage(memlogger):
    lg, records = memlogger
    before = MemoryLogger.get_memory_usage()
    lg.warning("x")
    after = MemoryLogger.get_memory_usage()
    m = PREFIX.match(records[0].getMessage())
    assert m is not None
    assert before <= int(m.group(1)) <= after
    assert records[0].getMessage()[m.end():] == "x"


@pytest.mark.parametrize("method, levelno", [
    ("debug", logging.DEBUG),
    ("info", logging.INFO),
    ("warning", logging.WARNING),
    ("error", logging.ERROR),
    ("critical", logging.CRITICAL),
])
def test_every_level_gets_the_prefix(memlogger, fake_memory, method, levelno):
    lg, records = memlogger
    getattr(lg, method)("msg %d", 7)
    assert records[0].levelno == levelno
    assert records[0].getMessage() == "[1234 MB] msg 7"


def test_log_with_explicit_level(memlogger, fake_memory):
    lg, records = memlogger
    lg.log(25, "custom")
    assert records[0].levelno == 25
    assert records[0].getMessage() == "[1234 MB] custom"


def test_filtered_message_skips_memory_lookup(memlogger, fake_memory):
    lg, records = memlogger
    lg.setLevel(logging.WARNING)
    lg.info("dropped")
    lg.debug("dropped")
    assert records == []
    assert fake_memory == []
    lg.warning("kept")
    assert [r.getMessage() for r in records] == ["[1234 MB] kept"]


def test_non_string_message_is_formatted_into_prefix(memlogger, fake_memory):
    lg, records = memlogger
    lg.info(42)
    lg.info({"a": 1})
    assert [r.getMessage() for r in records] == ["[1234 MB] 42",
                                                 "[1234 MB] {'a': 1}"]


def test_exc_info_passes_through(memlogger, fake_memory):
    lg, records = memlogger
    try:
        raise RuntimeError("boom")
    except RuntimeError:
        lg.exception("failed")
    rec = records[0]
    assert rec.levelno == logging.ERROR
    assert rec.getMessage() == "[1234 MB] failed"
    assert rec.exc_info[0] is RuntimeError


def test_extra_and_stack_info_pass_through(memlogger, fake_memory):
    lg, records = memlogger
    lg.info("x", extra={"macs3_field": 5}, stack_info=True)
    assert records[0].macs3_field == 5
    assert records[0].stack_info.startswith("Stack (most recent call last)")


# ------------------------------------
# MemoryLogger.get_memory_usage
# ------------------------------------

def test_get_memory_usage_is_int_megabytes_of_maxrss():
    # Linux reports ru_maxrss in KiB
    before = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    mb = MemoryLogger.get_memory_usage()
    after = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    assert type(mb) is int
    assert before <= mb <= after
    assert mb > 0


def test_get_memory_usage_callable_on_instance():
    lg = MemoryLogger("macs3test.instance")
    assert type(lg.get_memory_usage()) is int


FakeUsage = namedtuple("FakeUsage", "ru_maxrss")
FakeUname = namedtuple("FakeUname", "sysname")


@pytest.mark.parametrize("maxrss, expected", [
    (0, 0),
    (1023, 0),                  # truncates toward 0
    (1024, 1),
    (2048 * 1024 + 1023, 2048),
    (64 * 1024 * 1024, 65536),
])
def test_get_memory_usage_linux_kib_to_mib(monkeypatch, maxrss, expected):
    monkeypatch.setattr(resource, "getrusage",
                        lambda who: FakeUsage(maxrss))
    monkeypatch.setattr(os, "uname", lambda: FakeUname("Linux"))
    assert MemoryLogger.get_memory_usage() == expected


@pytest.mark.parametrize("maxrss, expected", [
    (1024 * 1024 - 1, 0),
    (1024 * 1024, 1),
    (3 * 1024**3 + 7, 3072),
])
def test_get_memory_usage_darwin_bytes_to_mib(monkeypatch, maxrss, expected):
    # macOS reports ru_maxrss in bytes
    monkeypatch.setattr(resource, "getrusage",
                        lambda who: FakeUsage(maxrss))
    monkeypatch.setattr(os, "uname", lambda: FakeUname("Darwin"))
    assert MemoryLogger.get_memory_usage() == expected


def test_get_memory_usage_asks_for_this_process(monkeypatch):
    seen = []

    def fake(who):
        seen.append(who)
        return FakeUsage(4096)
    monkeypatch.setattr(resource, "getrusage", fake)
    assert MemoryLogger.get_memory_usage() == 4
    assert seen == [resource.RUSAGE_SELF]


# ------------------------------------
# logging.basicConfig set at import
# ------------------------------------
# basicConfig only acts when the root logger has no handler, which is
# not the case inside pytest, so the format is checked in a fresh
# interpreter.

FORMAT_CODE = """
import logging
import MACS3.Utilities.Logger
lg = logging.getLogger("macs3test.format")
lg.debug("not shown")
lg.info("hello %s", "there")
lg.warning("warn %d", 3)
lg.error("err")
lg.critical("crit")
print(type(lg).__name__, logging.getLogger().level)
"""

LINE = re.compile(r"^(\w+) +@ (\d\d \w\w\w \d{4} \d\d:\d\d:\d\d): "
                  r"\[\d+ MB\] (.*) $")


@pytest.fixture(scope="module")
def format_run():
    return run_python(FORMAT_CODE)


def test_basicconfig_stdout_reports_class_and_root_level(format_run):
    assert format_run.returncode == 0, format_run.stderr
    assert format_run.stdout == "MemoryLogger 20\n"


def test_basicconfig_lines_and_level_padding(format_run):
    lines = format_run.stderr.splitlines()
    # levelname is left-justified to 5 characters, message is followed
    # by one space; INFO is the root level, so debug is not shown
    heads = [ln.split(" @ ")[0] for ln in lines]
    assert heads == ["INFO ", "WARNING", "ERROR", "CRITICAL"]
    messages = []
    for ln in lines:
        m = LINE.match(ln)
        assert m is not None, ln
        time.strptime(m.group(2), "%d %b %Y %H:%M:%S")
        messages.append((m.group(1), m.group(3)))
    assert messages == [("INFO", "hello there"), ("WARNING", "warn 3"),
                        ("ERROR", "err"), ("CRITICAL", "crit")]


def test_basicconfig_output_parsed_by_conftest_parser(format_run,
                                                      parse_log):
    assert parse_log(format_run.stderr) == [
        ("INFO", "hello there"), ("WARNING", "warn 3"), ("ERROR", "err"),
        ("CRITICAL", "crit")]


def test_basicconfig_writes_to_stderr_not_stdout(format_run):
    assert "hello" not in format_run.stdout
    assert "hello there" in format_run.stderr


def test_module_level_handler_targets_stderr_when_configured():
    # inside a fresh interpreter, the root handler is a StreamHandler on
    # sys.stderr with the documented format and date format
    code = ("import logging, sys; import MACS3.Utilities.Logger;"
            "h = logging.getLogger().handlers[0];"
            "print(type(h).__name__, h.stream is sys.stderr);"
            "print(repr(h.formatter._fmt)); print(repr(h.formatter.datefmt))")
    proc = run_python(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.splitlines() == [
        "StreamHandler True",
        repr('%(levelname)-5s @ %(asctime)s: %(message)s '),
        repr('%d %b %Y %H:%M:%S'),
    ]


def test_logger_module_exports_configured_logging():
    # MACS3 modules do `from MACS3.Utilities.Logger import logging`
    assert Logger.logging is logging

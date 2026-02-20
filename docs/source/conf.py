import os
import sys
from types import ModuleType
from unittest.mock import MagicMock

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# Ensure the repository root is discoverable for autodoc when building from
# a source checkout instead of an installed wheel.
sys.path.insert(0, os.path.abspath('../..'))

class _CythonType:
    """Minimal stub that behaves like a Cython type at runtime."""

    def __getitem__(self, _):
        return self

    def __call__(self, *args, **kwargs):
        return None


try:
    import cython  # noqa: F401
except ImportError:  # pragma: no cover - only used in docs builds
    cython = ModuleType("cython")
    sys.modules['cython'] = cython
    _decorator = lambda *args, **kwargs: (lambda obj: obj)
    for _name in ("cclass", "ccall", "cfunc", "locals", "returns"):
        setattr(cython, _name, _decorator)
    for _name in ("pointer", "const", "uchar", "ushort", "int", "long"):
        setattr(cython, _name, _CythonType())
    cython.declare = lambda *args, **kwargs: None

if not hasattr(cython, 'bytes'):
    cython.bytes = bytes
if not hasattr(cython, 'tuple'):
    cython.tuple = tuple
if not hasattr(cython, 'list'):
    cython.list = list
_cython_type_attrs = ("pointer", "const", "uchar", "ushort", "int", "long")
for _name in _cython_type_attrs:
    if not hasattr(cython, _name):
        setattr(cython, _name, _CythonType())

class _CImportStub(ModuleType):
    def __getattr__(self, name):
        if name.isupper():
            return 0.0
        return _fake_callable


def _fake_callable(*args, **kwargs):
    return None


_cimport_modules = [
    'cython.cimports',
    'cython.cimports.cpython',
    'cython.cimports.libc',
    'cython.cimports.libc.math',
    'cython.cimports.libc.stdint',
    'cython.cimports.libc.stdlib',
    'cython.cimports.libc.stdio',
    'cython.cimports.numpy',
    'cython.cimports.MACS3',
    'cython.cimports.MACS3.Signal',
    'cython.cimports.MACS3.Signal.cPosValCalculation',
    'cython.cimports.MACS3.IO',
    'cython.cimports.MACS3.Utilities',
    'cython.cimports.MACS3.Commands',
]

for _mod in _cimport_modules:
    module = _CImportStub(_mod)
    module.__path__ = []
    sys.modules[_mod] = module
    parent_name, _, child_name = _mod.rpartition('.')
    if parent_name:
        parent = sys.modules.setdefault(parent_name, _CImportStub(parent_name))
        setattr(parent, child_name, module)

sys.modules['cython'].cimports = sys.modules['cython.cimports']

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MACS3'
copyright = '2025, Tao Liu, Philippa Doherty'
author = 'Tao Liu, Philippa Doherty'
release = '3.0.4'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_parser',
 #   'sphinx.ext.autodoc',
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
  #  "sphinx_autodoc_typehints",
    #'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages'
    ]

autosummary_generate = True

autodoc_mock_imports = ["scipy", "cykhash", "scikit-learn", "hmmlearn"]
autodoc_typehints = "none"

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

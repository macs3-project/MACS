#!/usr/bin/env python3
"""Utility script to regenerate API docs and build HTML output with Sphinx.

The script wraps two common steps needed for documentation refreshes:

1. Run ``sphinx-apidoc`` so that the ``docs/source/api`` tree is rebuilt from
   the project's docstrings.
2. Invoke ``sphinx-build`` to render the documentation as HTML.

Usage
-----
Run the script from the project root::

    python scripts/build_docs.py

Optional arguments allow you to customise the build, see ``--help`` for details.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def run(cmd: list[str], cwd: Path) -> None:
    """Execute *cmd* in *cwd* and stream output to the terminal."""
    print("$", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract docstrings and build HTML docs with Sphinx.")
    parser.add_argument(
        "--apidoc",
        default="sphinx-apidoc",
        help="Executable used to generate the API sources (default: sphinx-apidoc)",
    )
    parser.add_argument(
        "--sphinx-build",
        default="sphinx-build",
        help="Executable used to render the documentation (default: sphinx-build)",
    )
    parser.add_argument(
        "--api-dir",
        default="docs/source/api",
        help="Destination directory for generated API sources (default: docs/source/api)",
    )
    parser.add_argument(
        "--builder",
        default="html",
        help="Sphinx builder name to invoke (default: html)",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove the API directory and builder output before regenerating",
    )
    parser.add_argument(
        "--exclude",
        nargs="*",
        default=["test", "docs", "build", "dist"],
        help="Additional paths to exclude from sphinx-apidoc generation",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = Path(__file__).resolve().parents[1]
    docs_root = project_root / "docs"
    source_dir = docs_root / "source"
    api_dir = (project_root / args.api_dir).resolve()
    build_dir = docs_root / "build" / args.builder

    if args.clean:
        if api_dir.exists():
            print(f"Cleaning API directory: {api_dir}")
            shutil.rmtree(api_dir)
        if build_dir.exists():
            print(f"Cleaning build directory: {build_dir}")
            shutil.rmtree(build_dir)

    api_dir.mkdir(parents=True, exist_ok=True)

    apidoc_cmd = [
        args.apidoc,
        "--force",
        "--module-first",
        "--output-dir",
        str(api_dir),
        str(project_root / "MACS3"),
    ]
    for item in args.exclude:
        apidoc_cmd.append(str(project_root / item))

    run(apidoc_cmd, cwd=project_root)

    build_cmd = [
        args.sphinx_build,
        "-b",
        args.builder,
        str(source_dir),
        str(build_dir),
    ]

    run(build_cmd, cwd=project_root)


if __name__ == "__main__":
    main()

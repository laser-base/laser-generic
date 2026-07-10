#!/usr/bin/env python3
"""
execute_notebooks.py — execute every .ipynb under docs/ into ``<exec_dir>``.

Mirrors each notebook's docs/-relative path under ``<exec_dir>`` so
``docs/concat_mkdocs.py`` can resolve mkdocs.yml nav entries by relative path
(e.g. ``tutorials/notebooks/01_SI...ipynb`` -> ``<exec_dir>/tutorials/notebooks/01_SI...ipynb``).

Always passes ``--allow-errors`` so nbconvert produces an output notebook for
each input even when cells raise — the error-gate (``docs/check_executed_nbs.py``)
is what decides whether to fail the pipeline.

Usage:
    python docs/execute_notebooks.py <exec_dir> [--timeout SECONDS] [--exclude SUBSTR[,SUBSTR...]]

Exit codes:
    0 — every notebook was processed by nbconvert (cell errors are tolerated)
    non-zero — nbconvert itself failed (e.g. kernel won't start) on some notebook
"""

import argparse
import subprocess
import sys
from pathlib import Path

from tqdm_strip import strip_notebook_file


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("exec_dir", help="output directory for executed notebooks")
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="per-cell execution timeout in seconds (default: 600)",
    )
    parser.add_argument(
        "--exclude",
        default="",
        help=(
            "comma-separated substrings; any notebook whose docs/-relative path "
            "contains one of these is skipped (intended for research-only "
            "notebooks that are maintained for manual exploration but aren't "
            "suitable for automated execution)."
        ),
    )
    args = parser.parse_args()

    docs = Path("docs")
    exec_dir = Path(args.exec_dir)
    exec_dir.mkdir(parents=True, exist_ok=True)

    excludes = [p.strip() for p in args.exclude.split(",") if p.strip()]

    all_notebooks = sorted(p for p in docs.rglob("*.ipynb") if ".ipynb_checkpoints" not in p.parts)
    if not all_notebooks:
        print(f"No notebooks found under {docs}/", file=sys.stderr)
        return 0

    notebooks = []
    skipped = []
    for nb in all_notebooks:
        rel_str = nb.relative_to(docs).as_posix()
        if any(pat in rel_str for pat in excludes):
            skipped.append(nb)
        else:
            notebooks.append(nb)

    if skipped:
        print(f"Skipping {len(skipped)} excluded notebook(s) (matched: {excludes}):")
        for nb in skipped:
            print(f"  - {nb}")

    print(f"Executing {len(notebooks)} notebook(s) -> {exec_dir}/ ...")
    for nb in notebooks:
        rel = nb.relative_to(docs)
        out_dir = exec_dir / rel.parent
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            sys.executable,
            "-m",
            "jupyter",
            "nbconvert",
            "--to",
            "notebook",
            "--execute",
            "--allow-errors",
            f"--ExecutePreprocessor.timeout={args.timeout}",
            "--output-dir",
            str(out_dir),
            "--output",
            nb.stem,
            str(nb),
        ]
        rc = subprocess.run(cmd).returncode
        if rc != 0:
            print(f"nbconvert exited {rc} for {nb}", file=sys.stderr)
            return rc

        # Strip tqdm progress-bar noise from cell outputs before the artifact
        # ships anywhere downstream. Doing this at execute time keeps the
        # canonical CI artifact clean by construction — MkDocs Deploy's
        # notebook overlay and Build Combined Doc's RAG concat both consume
        # it, so both get clean outputs without needing to filter separately.
        # See docs/tqdm_strip.py for the regex + shape it targets.
        out_path = out_dir / f"{nb.stem}.ipynb"
        n_removed = strip_notebook_file(out_path)
        if n_removed:
            print(f"  stripped {n_removed} tqdm-progress line(s) from {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

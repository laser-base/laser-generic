#!/usr/bin/env python3
"""
check_executed_nbs.py — report executed notebooks that contain error outputs.

Walks ``<exec_dir>`` recursively, opens each .ipynb, and lists any whose code
cells contain at least one cell output with ``output_type == "error"``.

Usage:
    python docs/check_executed_nbs.py <exec_dir> [--allow-errors]

Exit codes:
    0 — no error outputs found, OR --allow-errors was passed
    1 — at least one notebook has error outputs and --allow-errors was NOT passed

Prints a per-file summary either way so failures are always visible. When run
under GitHub Actions (``GITHUB_ACTIONS=true``), the summary uses ``::error::``
or ``::warning::`` annotations so failures surface on the run page.
"""

import argparse
import json
import os
import sys
from pathlib import Path


def has_error_cells(nb_path: Path) -> bool:
    with nb_path.open(encoding="utf-8") as f:
        nb = json.load(f)
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        if any(out.get("output_type") == "error" for out in cell.get("outputs", [])):
            return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("exec_dir", help="directory containing executed .ipynb files")
    parser.add_argument(
        "--allow-errors",
        action="store_true",
        help="report failures but exit 0 (useful for debug runs that still want the artifact)",
    )
    args = parser.parse_args()

    exec_dir = Path(args.exec_dir)
    if not exec_dir.is_dir():
        print(f"ERROR: {exec_dir} is not a directory")
        return 1

    executed = sorted(exec_dir.rglob("*.ipynb"))
    failed = [p for p in executed if has_error_cells(p)]
    print(f"Checked {len(executed)} executed notebook(s); {len(failed)} had errors.")

    if not failed:
        return 0

    in_ci = os.environ.get("GITHUB_ACTIONS", "").lower() == "true"
    annotation = ("::warning::" if args.allow_errors else "::error::") if in_ci else ""
    print(f"{annotation}{len(failed)} notebook(s) contain execution errors:")
    for p in failed:
        print(f"  - {p}")

    if args.allow_errors:
        return 0
    print("Fail-fast: pass --allow-errors to proceed anyway.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

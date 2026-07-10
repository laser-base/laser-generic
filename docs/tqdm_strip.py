#!/usr/bin/env python3
"""tqdm_strip.py — drop tqdm progress-bar noise from Jupyter notebook outputs.

laser-generic notebooks call ``model.run()`` with a tqdm progress display,
which nbconvert captures as many mid-run snapshots per cell:

    1,000,000 agents in 1 node(s):  37%|███▋      | 3712/10000 [00:14<00:22, 271.83it/s]

These lines add no semantic value, hurt RAG/retrieval quality (they dominate
the token budget in some sections), and produce ugly renderings on the docs
site. The regex here matches any line ending with the tqdm rate suffix
(``it/s]`` or its inverse ``s/it]``) — a tqdm-specific signature that rarely
appears in legitimate output.

Callers:

    - docs/execute_notebooks.py — invokes ``strip_notebook_file`` after each
      nbconvert run so the ``executed_nbs`` CI artifact ships clean.
    - docs/concat_mkdocs.py — imports ``TQDM_PROGRESS_RE`` for its own
      progress-run collapse logic (which produces "first \n ... \n last"
      summaries in the combined markdown for readability).

CLI:

    python docs/tqdm_strip.py path1.ipynb path2.ipynb ...       # strip in place
    python docs/tqdm_strip.py --check path1.ipynb path2.ipynb   # exit 1 if any file has tqdm lines
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Iterable, Tuple


# Matches lines ending with the tqdm rate marker. Catches:
#   - Named progress bars ("1,000,000 agents in 1 node(s): 37%|... | ... it/s]")
#   - Anonymous ones            ("100%|... 25/25 [02:15<00:00, 5.41s/it]")
#   - Sim-count descriptions    ("10/10 simulations: 45%|... | ... it/s]")
# The `(?<![A-Za-z])` lookbehind prevents matches inside words like "digits/second"
# ending in `s]`, and requires the rate marker to be at end-of-line to avoid
# clipping legitimate output that happens to mention "it/s" mid-sentence.
TQDM_PROGRESS_RE = re.compile(r"(?<![A-Za-z])(?:it/s|s/it)\]\s*$")


def _text_of(output: dict) -> Tuple[str, str]:
    """Return ``(field, joined_text)`` for whichever text field an output uses.

    ``stream`` outputs use ``text``; ``execute_result`` / ``display_data`` use
    ``data['text/plain']``. Both can be either a string or a list of strings
    (nbformat allows both). Returns ("", "") if the output has no text at all.
    """
    if "text" in output:
        text = output["text"]
        if isinstance(text, list):
            text = "".join(text)
        return "text", text
    data = output.get("data") or {}
    if "text/plain" in data:
        text = data["text/plain"]
        if isinstance(text, list):
            text = "".join(text)
        return "data.text/plain", text
    return "", ""


def _set_text(output: dict, field: str, new_text: str) -> None:
    """Write ``new_text`` back into the output under ``field``, splitting on
    newlines to match nbformat's preferred list-of-lines shape.

    nbformat lets you use either a plain string or a list of lines with
    embedded ``\\n``. To minimize noise in future diffs, we prefer the
    list-of-lines shape if the original was a list, else keep it a string.
    """
    if field == "text":
        original = output.get("text")
    else:
        original = (output.get("data") or {}).get("text/plain")

    if isinstance(original, list):
        # Preserve the split-on-newlines shape nbformat produces.
        lines = new_text.splitlines(keepends=True)
        if new_text and not new_text.endswith("\n"):
            # Preserve the final-line no-trailing-newline case
            pass
        value = lines
    else:
        value = new_text

    if field == "text":
        output["text"] = value
    else:
        output.setdefault("data", {})["text/plain"] = value


def strip_tqdm_from_notebook(nb: dict) -> int:
    """Mutate ``nb`` in place, dropping tqdm-progress lines from output cells.

    Returns the number of lines removed across all cells. Handles both
    ``stream`` outputs (``text`` field) and ``execute_result`` / ``display_data``
    outputs (``data.text/plain`` field). Never touches source cells, cell
    metadata, or image outputs.
    """
    removed = 0
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        for output in cell.get("outputs", []):
            field, text = _text_of(output)
            if not text:
                continue
            original_chunks = text.splitlines(keepends=True)
            kept_lines = [ln for ln in original_chunks if not TQDM_PROGRESS_RE.search(ln.rstrip("\r\n"))]
            n_removed = len(original_chunks) - len(kept_lines)
            if n_removed == 0:
                continue
            removed += n_removed
            _set_text(output, field, "".join(kept_lines))
    return removed


def strip_notebook_file(path: Path) -> int:
    """Strip tqdm crud from ``path`` in place. Returns the number of lines removed.

    No-ops (does not rewrite the file) when the notebook is already clean, so
    file mtimes and git diffs stay stable across repeated runs.
    """
    with path.open(encoding="utf-8") as f:
        nb = json.load(f)
    removed = strip_tqdm_from_notebook(nb)
    if removed > 0:
        with path.open("w", encoding="utf-8") as f:
            json.dump(nb, f, indent=1, ensure_ascii=False)
            f.write("\n")
    return removed


def count_tqdm_lines(path: Path) -> int:
    """Return the number of tqdm-progress lines in ``path`` (read-only)."""
    with path.open(encoding="utf-8") as f:
        nb = json.load(f)
    count = 0
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        for output in cell.get("outputs", []):
            _, text = _text_of(output)
            if not text:
                continue
            for line in text.splitlines():
                if TQDM_PROGRESS_RE.search(line):
                    count += 1
    return count


def _iter_paths(argv_paths: Iterable[str]) -> Iterable[Path]:
    for p in argv_paths:
        path = Path(p)
        if not path.exists():
            print(f"warning: {p} does not exist, skipping", file=sys.stderr)
            continue
        yield path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="+", help="notebook paths (.ipynb)")
    parser.add_argument(
        "--check",
        action="store_true",
        help="exit non-zero if any file has tqdm-progress lines (does not modify files)",
    )
    args = parser.parse_args()

    if args.check:
        offenders = []
        for path in _iter_paths(args.paths):
            n = count_tqdm_lines(path)
            if n > 0:
                offenders.append((path, n))
        if offenders:
            print("Notebooks with tqdm-progress noise in outputs:")
            for path, n in offenders:
                print(f"  {path}: {n} line(s)")
            print(
                "\nRun without --check to strip in place, or clear outputs and re-execute "
                "with `make docs-executed-nbs` (which now strips tqdm output automatically).",
                file=sys.stderr,
            )
            return 1
        print(f"All {len(args.paths)} notebook(s) clean.")
        return 0

    total_removed = 0
    changed = 0
    for path in _iter_paths(args.paths):
        n = strip_notebook_file(path)
        if n > 0:
            print(f"  stripped {n:>4} lines from {path}")
            total_removed += n
            changed += 1
    print(f"Done. Stripped {total_removed} line(s) across {changed} file(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Prototype: inject plt.figtext(...) captions into plot-generating code cells.

Reads tools/plot_descriptions/config.json (the fat manifest). For each entry:
  1. Extracts the takeaway sentence from the .md file using a "last bold
     containing a takeaway trigger" heuristic (same as tools/plot_descriptions.slim).
  2. Finds the target code cell in the notebook.
  3. Inserts a plt.figtext(...) line just before the final plt.show() in that
     cell (or appends at end if no show call).

Idempotent: no-ops if the exact figtext line is already present.
"""
import json
import re
import sys
from pathlib import Path


FIGTEXT_TRIGGERS = (
    "demonstrate", "takeaway", "confirms", "shows that",
    "validates", "verifies", "reproduces", "matches",
    "reveals", "consistent with", "recovers", "predicts",
    "successfully", "scales as", "agrees with", "sanity check",
)

MANUAL_TAKEAWAY = {
    "01_SI_nobirths_sanity.md":
        "The three curves overlay across the entire range, sanity-checking that "
        "the model preserves N = S + I at every tick and that cumulative "
        "incidence tracks running I_t minus the initial seed.",
}


def extract_takeaway(md_text: str, md_filename: str) -> str:
    if md_filename in MANUAL_TAKEAWAY:
        return MANUAL_TAKEAWAY[md_filename]
    bolds = list(re.finditer(r"\*\*([^*]+)\*\*", md_text))
    rich = [b.group(1).strip() for b in bolds if len(b.group(1).strip()) > 20]
    for b in reversed(rich):
        if any(t in b.lower() for t in FIGTEXT_TRIGGERS):
            return re.sub(r"\s+", " ", b).strip().lstrip("—:- ")
    m = re.search(r"\*\*(The takeaway:)\*\*\s*(.+?)(?:\n\n|\Z)", md_text, re.DOTALL)
    if m:
        return "The takeaway: " + re.sub(r"\s+", " ", m.group(2)).strip()
    if rich:
        return re.sub(r"\s+", " ", rich[-1]).strip().lstrip("—:- ")
    return None


def build_figtext_line(takeaway: str) -> str:
    text = takeaway.strip().strip("*").strip()
    py_repr = repr(text)
    return (
        f"plt.figtext(0.5, -0.05, {py_repr}, "
        f'ha="center", va="top", wrap=True, fontsize=8)'
    )


def inject_before_show(source_lines, figtext_line):
    src = "".join(source_lines)
    if figtext_line in src:
        return source_lines, False
    lines = src.splitlines(keepends=True)
    last_show_idx = None
    for i, ln in enumerate(lines):
        if re.match(r"^\s*plt\.show\s*\(", ln):
            last_show_idx = i
    if last_show_idx is None:
        if lines and not lines[-1].endswith("\n"):
            lines[-1] = lines[-1] + "\n"
        lines.append(figtext_line + "\n")
    else:
        indent = re.match(r"^(\s*)", lines[last_show_idx]).group(1)
        lines.insert(last_show_idx, f"{indent}{figtext_line}\n")
    return lines, True


def find_code_cell(cells, ins):
    if "after_cell_containing" in ins:
        needle = ins["after_cell_containing"]
        matches = [
            i for i, c in enumerate(cells)
            if c.get("cell_type") == "code" and needle in "".join(c.get("source", []))
        ]
    elif "after_cell_id" in ins:
        cid = ins["after_cell_id"]
        matches = [i for i, c in enumerate(cells) if c.get("id") == cid]
    else:
        raise ValueError("insert needs after_cell_containing or after_cell_id")
    if len(matches) != 1:
        raise LookupError(f"got {len(matches)} matches for {ins!r}")
    return matches[0]


def strip_description_cells(cells):
    """Remove markdown cells whose first line starts with '### Reading the'."""
    before = len(cells)
    cells[:] = [
        c for c in cells
        if not (
            c.get("cell_type") == "markdown"
            and "".join(c.get("source", [])).lstrip().startswith("### Reading the")
        )
    ]
    return before - len(cells)


def apply(manifest_path: Path):
    base = manifest_path.parent
    entries = json.loads(manifest_path.read_text(encoding="utf-8"))

    for entry in entries:
        nb_path = (base / entry["notebook"]).resolve()
        data = json.loads(nb_path.read_text(encoding="utf-8"))
        cells = data["cells"]
        n_removed = strip_description_cells(cells)
        changed = n_removed > 0
        for ins in entry["inserts"]:
            idx = find_code_cell(cells, ins)
            md_text = (base / ins["md_file"]).read_text(encoding="utf-8")
            takeaway = extract_takeaway(md_text, ins["md_file"])
            if takeaway is None:
                print(f"WARN {nb_path.name} — no takeaway extracted for {ins['md_file']}")
                continue
            figtext_line = build_figtext_line(takeaway)
            new_lines, did_change = inject_before_show(cells[idx]["source"], figtext_line)
            if did_change:
                cells[idx]["source"] = new_lines
                changed = True
        if changed:
            nb_path.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
            print(f"CHANGED  {nb_path.name}  (stripped {n_removed} md, added figtext)")
        else:
            print(f"ok       {nb_path.name}")


if __name__ == "__main__":
    apply(Path(sys.argv[1] if len(sys.argv) > 1 else "tools/plot_descriptions/config.json").resolve())

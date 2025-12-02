#!/usr/bin/env python3
import pathlib

root = pathlib.Path("sphinx/source/converted")
ignore = {"reference", "notebooks"}
out_lines = [
    ".. toctree::",
    "   :maxdepth: 2",
    "   :caption: Getting Started",
    "",
]
for rst in sorted(root.rglob("*.rst")):
    if any(part in ignore for part in rst.parts):
        continue
    rel = f"source/converted/{rst.relative_to(root).with_suffix('')}"
    out_lines.append(f"   {rel}")
out = "\n".join(out_lines)
print(out)

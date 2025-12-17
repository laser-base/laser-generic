from pathlib import Path
import re

DOCS_ROOT = Path("../docs")
RST_FILE = Path("index.rst")
CONVERTED_PREFIX = "source/converted"

EXCLUDE_DIRS = {"reference", "tutorials", "__pycache__"}
EXCLUDE_FILES = {".DS_Store", "api.md.gz"}


def collect_markdown_paths():
    paths = []
    for path in DOCS_ROOT.rglob("*.md"):
        if any(part in EXCLUDE_DIRS for part in path.parts):
            continue
        if path.name in EXCLUDE_FILES:
            continue
        rel = path.relative_to(DOCS_ROOT).with_suffix("")
        converted = f"{CONVERTED_PREFIX}/{rel.as_posix()}"
        paths.append(converted)
    return sorted(paths, key=lambda p: (not p.endswith("index"), p))


def build_toctree_block(paths):
    lines = [".. toctree::", "   :maxdepth: 2", "   :caption: Getting Started", ""]
    lines += [f"   {p}" for p in paths]
    return "\n".join(lines)


def replace_toctree_block(rst_path, new_block):
    content = rst_path.read_text()

    # Match from the "Getting Started" toctree through its content
    pattern = re.compile(r"\.\. toctree::\n\s+:maxdepth: 2\n\s+:caption: Getting Started\n\n(   .+\n)+", re.MULTILINE)

    # Ensure the pattern matches before performing the substitution to avoid silently
    # writing the file unchanged when the expected block is missing.
    if not pattern.search(content):
        print(f"⚠️ No matching 'Getting Started' toctree block found in {rst_path}; file left unchanged.")
        return
    new_content = pattern.sub(new_block + "\n\n", content)
    rst_path.write_text(new_content)
    print(f"✅ Updated: {rst_path}")


if __name__ == "__main__":
    md_paths = collect_markdown_paths()
    new_block = build_toctree_block(md_paths)
    replace_toctree_block(RST_FILE, new_block)

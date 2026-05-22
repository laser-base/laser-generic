#!/usr/bin/env python3
"""
concat_mkdocs.py — flatten the laser-generic MkDocs build into a single markdown file.

The combined output contains three kinds of content:
  1. Human-authored .md pages, in mkdocs.yml nav order, recovered from the rendered HTML site.
  2. API reference pages (mkdocstrings output), inserted at the nav position of ``reference/``.
  3. Notebook tutorials, converted directly from the executed .ipynb files for output fidelity.

Usage:
    python docs/concat_mkdocs.py <mkdocs_site_dir> <executed_notebooks_dir> <output_file>

Example:
    python docs/concat_mkdocs.py site dist/executed_nbs dist/combined_mkdocs.md
"""

import json
import re
import sys
from pathlib import Path

import markdownify
import yaml
from bs4 import BeautifulSoup

try:
    import laser.generic

    _LASER_GENERIC_VERSION = laser.generic.__version__
except Exception:
    _LASER_GENERIC_VERSION = "unknown"


MIN_SECTION_CHARS = 150  # drop near-empty placeholder/index pages
EXPECTED_MIN_MAIN_PAGES = 8  # fail loud if the site walk looks broken


def extract_markdown(html_path: Path) -> str:
    """Pull the main article from a MkDocs HTML page and convert it back to markdown."""
    text = html_path.read_text(encoding="utf-8", errors="replace")
    soup = BeautifulSoup(text, "html.parser")

    content = (
        soup.find("article")
        or soup.find("div", role="main")
        or soup.find("div", {"class": "md-content"})
        or soup.find("body")
    )
    if content is None:
        return ""

    for tag in content.find_all(["nav", "footer", "script", "style"]):
        tag.decompose()
    for tag in content.find_all(
        class_=["md-nav", "md-sidebar", "md-search", "md-header", "md-footer", "headerlink", "md-breadcrumb"]
    ):
        tag.decompose()

    # MkDocs Material renders syntax-highlighted code blocks as a two-column
    # table (linenos | code). markdownify mangles that into a garbled markdown
    # table. Replace each such table with just the <pre> from the code cell.
    for table in content.find_all("table", class_="highlighttable"):
        code_td = table.find("td", class_="code")
        if code_td:
            pre = code_td.find("pre")
            if pre:
                table.replace_with(pre)
            else:
                table.decompose()
        else:
            table.decompose()

    md = markdownify.markdownify(
        str(content),
        heading_style=markdownify.ATX,
        code_language="python",
        strip=["a"],
    )
    md = re.sub(r"\n{3,}", "\n\n", md)
    return md.strip()


def notebook_to_markdown(nb_path: Path) -> str:
    """Convert an executed Jupyter notebook to plain markdown.

    Markdown cells included as-is; code cells as fenced ``python`` blocks;
    text outputs appended after each code fence; image outputs skipped.
    """
    with nb_path.open(encoding="utf-8") as f:
        nb = json.load(f)

    parts = []
    for cell in nb["cells"]:
        if cell["cell_type"] == "markdown":
            src = "".join(cell["source"]).strip()
            if src:
                parts.append(src)
            continue

        if cell["cell_type"] != "code":
            continue

        src = "".join(cell["source"]).strip()
        if not src:
            continue
        parts.append(f"```python\n{src}\n```")

        for output in cell.get("outputs", []):
            otype = output.get("output_type", "")
            if otype not in ("stream", "execute_result", "display_data"):
                continue
            text = output.get("text") or output.get("data", {}).get("text/plain", [])
            if isinstance(text, list):
                text = "".join(text)
            if not text or not text.strip():
                continue
            cleaned = re.sub(r"<[^>]+>", "", text).strip()
            if cleaned:
                parts.append(f"```\n{cleaned}\n```")

    md = "\n\n".join(parts)
    return re.sub(r"\n{3,}", "\n\n", md).strip()


def iter_nav_pages(nav):
    """Yield each leaf entry from a parsed mkdocs.yml nav (strings only)."""
    if isinstance(nav, str):
        yield nav
    elif isinstance(nav, dict):
        for v in nav.values():
            yield from iter_nav_pages(v)
    elif isinstance(nav, list):
        for item in nav:
            yield from iter_nav_pages(item)


def docs_path_to_site_path(rel: str, site_dir: Path) -> Path:
    """Translate a docs/-relative .md or .ipynb path to its rendered site/<...>/index.html."""
    p = Path(rel)
    stem = p.with_suffix("")
    if stem.name == "index":
        return site_dir / stem.parent / "index.html"
    return site_dir / stem / "index.html"


def get_reference_pages(site_dir: Path):
    """All reference/ pages from the rendered site, sorted for stable ordering."""
    ref_dir = site_dir / "reference"
    if not ref_dir.exists():
        return []
    return sorted(ref_dir.rglob("index.html"))


def append_section(parts: list, label, md: str) -> bool:
    """Append an extracted section if it has enough content. Returns True if appended."""
    if md and len(md) >= MIN_SECTION_CHARS:
        parts.append(f"\n\n---\n<!-- {label} -->\n\n{md}")
        return True
    return False


def concat(mkdocs_dir: str, notebooks_dir: str, output_file: str):
    site_dir = Path(mkdocs_dir)
    nb_dir = Path(notebooks_dir)
    mkdocs_yml = Path("mkdocs.yml")

    if not mkdocs_yml.exists():
        raise RuntimeError(f"mkdocs.yml not found at {mkdocs_yml.resolve()} (run from repo root).")

    with mkdocs_yml.open(encoding="utf-8") as f:
        # mkdocs.yml uses !!python/name tags; safe_load rejects those, so use unsafe_load.
        # We control this file in the repo, so this is acceptable.
        config = yaml.unsafe_load(f)

    nav_entries = list(iter_nav_pages(config.get("nav", [])))
    parts = [f"# laser-generic documentation\n\n**laser-generic version: {_LASER_GENERIC_VERSION}**"]

    main_included = nb_included = ref_included = 0
    skipped = 0

    print(f"laser-generic version: {_LASER_GENERIC_VERSION}")
    print(f"=== Walking {len(nav_entries)} nav entries from mkdocs.yml ===")

    for entry in nav_entries:
        if entry.endswith(".md"):
            page = docs_path_to_site_path(entry, site_dir)
            label = page.relative_to(site_dir) if page.exists() else entry
            if not page.exists():
                print(f"  skip (not built): {entry}")
                skipped += 1
                continue
            md = extract_markdown(page)
            if append_section(parts, label, md):
                print(f"  ok (main): {label}")
                main_included += 1
            else:
                print(f"  skip (empty/too short): {label}")
                skipped += 1

        elif entry.endswith(".ipynb"):
            nb_path = nb_dir / entry  # preserve docs/-relative layout under executed dir
            label = f"{entry} (executed)"
            if nb_path.exists():
                md = notebook_to_markdown(nb_path)
            else:
                # Fall back to the rendered HTML so we don't silently drop a notebook page.
                fallback = docs_path_to_site_path(entry, site_dir)
                if not fallback.exists():
                    print(f"  skip (no executed nb, no rendered html): {entry}")
                    skipped += 1
                    continue
                md = extract_markdown(fallback)
                label = f"{fallback.relative_to(site_dir)} (html fallback)"

            if append_section(parts, label, md):
                print(f"  ok (notebook): {label}")
                nb_included += 1
            else:
                print(f"  skip (empty/too short): {label}")
                skipped += 1

        elif entry.rstrip("/") == "reference":
            # mkdocstrings + literate-nav drop generated API pages here.
            ref_pages = get_reference_pages(site_dir)
            print(f"  -- expanding reference/: {len(ref_pages)} pages")
            for path in ref_pages:
                rel = path.relative_to(site_dir)
                md = extract_markdown(path)
                if append_section(parts, rel, md):
                    print(f"    ok (reference): {rel}")
                    ref_included += 1
                else:
                    skipped += 1

        else:
            # Unknown shape (some other literate-nav dir, etc.) — leave for future work.
            print(f"  skip (unhandled nav entry): {entry}")
            skipped += 1

    if main_included < EXPECTED_MIN_MAIN_PAGES:
        raise RuntimeError(
            f"Only {main_included} main pages contributed content "
            f"(expected at least {EXPECTED_MIN_MAIN_PAGES}). Did `mkdocs build` finish cleanly?"
        )
    if ref_included == 0:
        print("WARNING: 0 reference pages included — was the API reference built?")
    if nb_included == 0:
        print("WARNING: 0 notebook pages included — were notebooks executed into the dir?")

    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(parts), encoding="utf-8")

    total = main_included + nb_included + ref_included
    print(f"\nWrote {total} sections ({skipped} skipped) -> {output_file}")
    print(f"  main: {main_included}, notebooks: {nb_included}, reference: {ref_included}")
    print(f"Output size: {output.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python docs/concat_mkdocs.py <mkdocs_site_dir> <executed_notebooks_dir> <output_file>")
        sys.exit(1)
    concat(sys.argv[1], sys.argv[2], sys.argv[3])

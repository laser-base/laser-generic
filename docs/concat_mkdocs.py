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

# laser-generic notebooks call tqdm with a description like
# "1,000,000 agents in 1 node(s)", which produces many "<n> agents in <m>
# node(s):  37%|...| ..." progress snapshots in cell outputs. They add no
# semantic value and hurt RAG/retrieval quality, so we drop matching lines
# at output-build time. Pattern matches the description + percent token.
_TQDM_PROGRESS_RE = re.compile(r"agents in \d+ node\(s\):\s*\d+%")


class _NavLoader(yaml.SafeLoader):
    """SafeLoader that tolerates mkdocs.yml ``!!python/name:`` tags.

    mkdocs.yml references Python callables (emoji generators, mermaid fence
    handlers, etc.) via ``!!python/name:`` tags. yaml.unsafe_load() resolves
    those — and would happily import arbitrary modules from an untrusted
    mkdocs.yml. We only walk ``nav:`` for string paths here, so the actual
    Python object is irrelevant; this loader replaces each such tag with an
    opaque sentinel string and lets the rest of the file parse safely.
    """


def _python_name_placeholder(loader, suffix, node):
    return f"<python/name:{suffix}>"


_NavLoader.add_multi_constructor("tag:yaml.org,2002:python/name:", _python_name_placeholder)


def extract_markdown(html_path: Path) -> str:
    """Pull the main article from a MkDocs HTML page and convert it back to markdown."""
    text = html_path.read_text(encoding="utf-8", errors="replace")
    soup = BeautifulSoup(text, "html.parser")

    content = soup.find("article") or soup.find("div", role="main") or soup.find("div", {"class": "md-content"}) or soup.find("body")
    if content is None:
        return ""

    for tag in content.find_all(["nav", "footer", "script", "style"]):
        tag.decompose()
    for tag in content.find_all(class_=["md-nav", "md-sidebar", "md-search", "md-header", "md-footer", "headerlink", "md-breadcrumb"]):
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

    # Preserve <a> tags so content cross-links and external URLs survive in
    # the combined markdown. Unwanted decorative anchors (headerlinks, nav
    # links, sidebar links) have already been decomposed above.
    md = markdownify.markdownify(
        str(content),
        heading_style=markdownify.ATX,
        code_language="python",
    )
    md = re.sub(r"\n{3,}", "\n\n", md)
    return md.strip()


def _collapse_progress_runs(blocks):
    """Collapse runs of consecutive ``progress`` blocks into one summarized block.

    A run of length N>=3 becomes a single block "<first>\\n...\\n<last>";
    a run of length 2 is joined as-is (no middle to ellide); a run of length
    1 is passed through unchanged. ``other`` blocks are emitted untouched.
    """
    result = []
    i = 0
    while i < len(blocks):
        if blocks[i][0] != "progress":
            result.append(blocks[i])
            i += 1
            continue
        j = i + 1
        while j < len(blocks) and blocks[j][0] == "progress":
            j += 1
        run = [b[1] for b in blocks[i:j]]
        if len(run) == 1:
            combined = run[0]
        elif len(run) == 2:
            combined = f"{run[0]}\n{run[1]}"
        else:
            combined = f"{run[0]}\n...\n{run[-1]}"
        result.append(("progress", combined))
        i = j
    return result


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

        # Classify each output as "progress" (purely tqdm progress lines) or
        # "other". A run of consecutive progress outputs is then collapsed to
        # "first \n ... \n last" so the reader keeps the start/end context
        # without the dozens of intermediate snapshots that dominate output.
        cell_blocks = []  # list of (kind, text) tuples
        for output in cell.get("outputs", []):
            otype = output.get("output_type", "")
            if otype not in ("stream", "execute_result", "display_data"):
                continue
            text = output.get("text") or output.get("data", {}).get("text/plain", [])
            if isinstance(text, list):
                text = "".join(text)
            if not text or not text.strip():
                continue

            # Note: deliberately no HTML-tag stripping here. Output is wrapped
            # in a triple-backtick fence, so any literal "<...>" content (e.g.
            # "<function foo at 0x...>", "<class 'Foo'>", "<lambda>") renders
            # as plain text — and stripping it would silently delete legitimate
            # Python reprs that carry information for the RAG corpus.
            non_empty_lines = [ln for ln in text.splitlines() if ln.strip()]
            if not non_empty_lines:
                continue

            if all(_TQDM_PROGRESS_RE.search(ln) for ln in non_empty_lines):
                cell_blocks.append(("progress", "\n".join(non_empty_lines).strip()))
            else:
                kept = "\n".join(ln for ln in text.splitlines() if not _TQDM_PROGRESS_RE.search(ln)).strip()
                if kept:
                    cell_blocks.append(("other", kept))

        for _, body in _collapse_progress_runs(cell_blocks):
            parts.append(f"```\n{body}\n```")

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


def _is_reference_entry(entry: str) -> bool:
    """True if a mkdocs.yml nav entry points anywhere under the ``reference/`` tree.

    Catches all of the common shapes — ``reference/``, ``reference/SUMMARY.md``
    (literate-nav idiom), ``reference/index.md``, ``reference/<anything>`` —
    so the API reference expansion isn't accidentally treated as a single .md
    page by the .md branch of the nav walk.
    """
    parts = Path(entry).parts
    return bool(parts) and parts[0] == "reference"


def concat(mkdocs_dir: str, notebooks_dir: str, output_file: str):
    site_dir = Path(mkdocs_dir)
    nb_dir = Path(notebooks_dir)
    mkdocs_yml = Path("mkdocs.yml")

    if not mkdocs_yml.exists():
        raise RuntimeError(f"mkdocs.yml not found at {mkdocs_yml.resolve()} (run from repo root).")

    with mkdocs_yml.open(encoding="utf-8") as f:
        config = yaml.load(f, Loader=_NavLoader)

    nav_entries = list(iter_nav_pages(config.get("nav", [])))
    parts = [f"# laser-generic documentation\n\n**laser-generic version: {_LASER_GENERIC_VERSION}**"]

    main_included = nb_included = ref_included = 0
    skipped = 0

    print(f"laser-generic version: {_LASER_GENERIC_VERSION}")
    print(f"=== Walking {len(nav_entries)} nav entries from mkdocs.yml ===")

    ref_expanded = False  # reference/ pages are expanded once, in nav position
    for entry in nav_entries:
        # Reference check goes first so reference/SUMMARY.md or reference/index.md
        # don't get swallowed by the .md branch — any entry under reference/
        # triggers the full mkdocstrings expansion regardless of its file shape.
        if _is_reference_entry(entry):
            if ref_expanded:
                # Multiple nav entries under reference/ all collapse into a single
                # expansion (the first one we hit) — no point re-emitting pages.
                continue
            ref_expanded = True
            ref_pages = get_reference_pages(site_dir)
            print(f"  -- expanding reference/ (triggered by nav entry '{entry}'): {len(ref_pages)} pages")
            for path in ref_pages:
                rel = path.relative_to(site_dir)
                md = extract_markdown(path)
                if append_section(parts, rel, md):
                    print(f"    ok (reference): {rel}")
                    ref_included += 1
                else:
                    skipped += 1

        elif entry.endswith(".md"):
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

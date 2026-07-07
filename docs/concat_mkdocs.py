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


def _env_placeholder(loader, node):
    """Handle mkdocs' ``!ENV [VAR_NAME, default]`` env-var directive.

    PR #220 introduced ``!ENV`` in mkdocs.yml to toggle notebook execution
    between local and CI builds. That tag is MkDocs-specific; plain PyYAML
    (which this loader is built on) doesn't know it and aborts the whole
    parse with a ConstructorError. Return the declared default (last element
    of the sequence) so the file parses; nothing downstream in
    concat_mkdocs.py consumes the value — it only walks ``nav:`` for
    string paths.
    """
    if isinstance(node, yaml.SequenceNode):
        seq = loader.construct_sequence(node)
        return seq[-1] if len(seq) >= 2 else ""
    return ""


_NavLoader.add_constructor("!ENV", _env_placeholder)


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

    # Notebook pages rendered to HTML (the .ipynb fallback path) carry Jupyter
    # cell chrome that pollutes the markdown: "In [N]:"/"Out [N]:" prompts, the
    # "Copied!" copy-button widget, and — worst — a hidden raw-source copy of
    # every code cell (clipboard-copy-txt) that duplicates the highlighted code
    # we actually keep. Strip all three so the fallback reads like a clean notebook.
    _NOTEBOOK_CHROME = ["jp-InputPrompt", "jp-OutputPrompt", "zeroclipboard-container", "clipboard-copy-txt"]
    for tag in content.find_all(class_=_NOTEBOOK_CHROME):
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
    text outputs appended after each code fence; image outputs become a
    visible caption placeholder.
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

            data = output.get("data", {})
            # Image outputs (matplotlib figures, etc.) can't live in a text
            # corpus, but silently dropping them erases the fact that the cell
            # renders a plot at this point in the narrative. Emit a visible
            # caption instead. Angle brackets are stripped from the figure's
            # text/plain repr so it isn't mistaken for an HTML tag and hidden.
            if any(key.startswith("image/") for key in data):
                repr_txt = "".join(data.get("text/plain", [])).strip()
                caption = repr_txt.replace("<", "").replace(">", "") or "image"
                cell_blocks.append(("figure", f"_[Figure output omitted: {caption}]_"))
                continue

            text = output.get("text") or data.get("text/plain", [])
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

        # Figure captions are plain markdown lines; text/progress outputs go in
        # fenced blocks. Keep them distinct so captions stay readable prose.
        for kind, body in _collapse_progress_runs(cell_blocks):
            if kind == "figure":
                parts.append(body)
            else:
                parts.append(f"```\n{body}\n```")

    md = "\n\n".join(parts)
    md = re.sub(r"\n{3,}", "\n\n", md).strip()

    # Demote all headings by one level so:
    #   - notebook titles (originally H1) become H2, sitting cleanly under the
    #     combined doc's single H1 ("# laser-generic documentation"),
    #   - notebook subsections (originally H2 like "## Larger test suite")
    #     become H3, which the RAG header-splitter catches structurally,
    #   - nested notebook sections stay proportionally deeper.
    # Levels are clamped at H6 (Markdown's max) to avoid emitting invalid H7+.
    #
    # CRITICAL: skip lines inside fenced code blocks — Python comments like
    # ``# %%capture`` or ``# TODO`` would otherwise be misread as headings and
    # demoted, corrupting the code AND emitting bogus H2 section boundaries.
    out_lines = []
    in_fence = False
    for line in md.split("\n"):
        stripped = line.lstrip()
        if stripped.startswith("```"):
            in_fence = not in_fence
            out_lines.append(line)
            continue
        if not in_fence:
            m = re.match(r"^(#{1,5})(?= )", line)
            if m:
                line = "#" * (len(m.group(1)) + 1) + line[len(m.group(1)) :]
        out_lines.append(line)
    return "\n".join(out_lines)


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


def published_url(site_rel: Path, site_url: str) -> str:
    """Map a site-relative ``.../index.html`` path back to its public URL.

    MkDocs serves ``foo/bar/index.html`` at ``<site_url>/foo/bar/``, so we drop
    the trailing ``index.html`` and keep the directory with a trailing slash.
    Returns "" when mkdocs.yml has no site_url (no link can be built).
    """
    if not site_url:
        return ""
    parts = site_rel.parts
    dir_parts = parts[:-1] if parts and parts[-1] == "index.html" else parts
    path = "/".join(dir_parts)
    return f"{site_url}/{path}/" if path else f"{site_url}/"


def get_reference_pages(site_dir: Path):
    """All reference/ pages from the rendered site, sorted for stable ordering."""
    ref_dir = site_dir / "reference"
    if not ref_dir.exists():
        return []
    return sorted(ref_dir.rglob("index.html"))


def append_section(parts: list, source, md: str, url: str = "") -> bool:
    """Append an extracted section if it has enough content. Returns True if appended.

    Provenance is emitted as visible markdown (a ``**Source:**`` line plus an
    optional link to the published page) rather than an HTML comment, because
    many RAG/markdown loaders strip HTML comments before chunking — which would
    otherwise discard the one breadcrumb tying each chunk back to its origin.
    """
    if not md or len(md) < MIN_SECTION_CHARS:
        return False
    source_path = source.as_posix() if isinstance(source, Path) else str(source)
    provenance = f"**Source:** `{source_path}`"
    if url:
        provenance += f" · [View page online]({url})"
    parts.append(f"\n\n---\n\n{provenance}\n\n{md}")
    return True


# Soft cap for section size beyond which the RAG chunker starts making lossy
# splits (fracturing fenced code examples across chunks). Not a hard error —
# just a warning that surfaces which content needs subheadings or cell splits.
_SECTION_SOFT_MAX = 5000


def _find_source_before(combined: str, pos: int) -> str:
    """Walk backward from ``pos`` to the nearest ``**Source:** `<path>`
    breadcrumb (emitted by :func:`append_section`) and return that path.
    Returns ``"?"`` if no breadcrumb is found — shouldn't happen given
    every section is prefixed with one, but degrade gracefully.
    """
    m = None
    for candidate in re.finditer(r"\*\*Source:\*\* `([^`]+)`", combined[:pos]):
        m = candidate
    return m.group(1) if m else "?"


def _warn_oversized_sections(combined: str) -> int:
    """Report H2+ sections that exceed the soft cap, with source-file provenance.

    Downstream (see laser-mcp/ingest.py) uses MarkdownHeaderTextSplitter on
    H1/H2/H3/H4 then a 1200-char character splitter for anything still too big.
    Once a section is much larger than a few chunks the character splitter
    falls through its code-block-aware separators to line-level splits — which
    tear fenced code examples across chunks and hurt retrieval quality. Warn
    so authors can add subheadings (H3/H4/H5) or split the underlying notebook
    cells into smaller logical sections.

    Each warning line names the source file (notebook or .md) that contributed
    the section, resolved via the ``**Source:**`` breadcrumbs that
    :func:`append_section` emits — so authors can go straight to the file
    that needs subheadings without grepping.
    """
    positions = [(m.start(), m.group(0).strip()) for m in re.finditer(r"^##+ [^\n]+", combined, re.MULTILINE)]
    positions.append((len(combined), None))
    warnings = []
    for (start, heading), (end, _) in zip(positions[:-1], positions[1:]):
        size = end - start
        if size > _SECTION_SOFT_MAX:
            source = _find_source_before(combined, start)
            warnings.append((size, heading, source))
    if warnings:
        warnings.sort(reverse=True)
        print(f"\n  WARNING: {len(warnings)} section(s) exceed {_SECTION_SOFT_MAX:,} chars (RAG chunking may fracture code examples):")
        for size, heading, source in warnings[:15]:
            print(f"    {size:>7,} chars  {heading[:60]:60}  {source}")
        if len(warnings) > 15:
            print(f"    ...{len(warnings) - 15} more")
        print(f"  Add subheadings (H3/H4/H5) or split cells so no section exceeds {_SECTION_SOFT_MAX:,} chars.")
    return len(warnings)


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
    site_url = (config.get("site_url") or "").rstrip("/")  # base for per-section "View page" links
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
                if append_section(parts, rel, md, published_url(rel, site_url)):
                    print(f"    ok (reference): {rel}")
                    ref_included += 1
                else:
                    skipped += 1

        elif entry.endswith(".md"):
            page = docs_path_to_site_path(entry, site_dir)
            if not page.exists():
                print(f"  skip (not built): {entry}")
                skipped += 1
                continue
            label = page.relative_to(site_dir)
            md = extract_markdown(page)
            if append_section(parts, label, md, published_url(label, site_url)):
                print(f"  ok (main): {label}")
                main_included += 1
            else:
                print(f"  skip (empty/too short): {label}")
                skipped += 1

        elif entry.endswith(".ipynb"):
            nb_path = nb_dir / entry  # preserve docs/-relative layout under executed dir
            label = f"{entry} (executed)"
            url = published_url(docs_path_to_site_path(entry, site_dir).relative_to(site_dir), site_url)
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
                fb_rel = fallback.relative_to(site_dir)
                label = f"{fb_rel.as_posix()} (html fallback)"
                url = published_url(fb_rel, site_url)

            if append_section(parts, label, md, url):
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

    combined = "\n".join(parts)
    _warn_oversized_sections(combined)

    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(combined, encoding="utf-8")

    total = main_included + nb_included + ref_included
    print(f"\nWrote {total} sections ({skipped} skipped) -> {output_file}")
    print(f"  main: {main_included}, notebooks: {nb_included}, reference: {ref_included}")
    print(f"Output size: {output.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python docs/concat_mkdocs.py <mkdocs_site_dir> <executed_notebooks_dir> <output_file>")
        sys.exit(1)
    concat(sys.argv[1], sys.argv[2], sys.argv[3])

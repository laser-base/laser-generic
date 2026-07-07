# laser-generic documentation pipeline.
#
# Lets you reproduce locally exactly what .github/workflows/build-combined-doc.yml
# does in CI. The "all-in-one" target is `docs-jenner`, which reads the
# committed notebook outputs, builds the MkDocs site, and flattens both into
# a single combined_mkdocs.md suitable for RAG / MCP ingestion. It does NOT
# re-execute the notebooks — that's the ~30-minute step and it's not needed
# to build the corpus, since every notebook is committed with its outputs.
#
# If you want to re-execute notebooks (either for validation, or because the
# committed outputs have drifted from the current API), use `docs-jenner-execute`
# instead — that runs the historical full pipeline (execute + error-gate + build
# + concat). Or run `docs-check-nbs` on its own for validation without a full
# doc build.

# All recipes are single-command invocations (echo or python script), so they
# work identically under bash, cmd.exe, and PowerShell — no SHELL override
# needed. Multi-step logic lives in docs/*.py helpers instead of inline shell.

.PHONY: help docs-install docs-build docs-executed-nbs docs-check-nbs docs-jenner docs-jenner-execute clean-docs

PYTHON           ?= python
SITE_DIR         ?= site
EXEC_DIR         ?= dist/executed_nbs
COMBINED         ?= dist/combined_mkdocs.md
NB_TIMEOUT       ?= 600
ALLOW_NB_ERRORS  ?= 0
# Comma-separated substrings of docs/-relative paths that should NOT be executed
# by the doc-build pipeline. EW_analysis is a research-only notebook maintained
# for manual exploration — its dependencies and runtime aren't suitable for the
# automated execute-and-check flow, so it's left out by default.
NB_EXCLUDE       ?= EW_analysis

help:
	@echo "laser-generic documentation targets"
	@echo "===================================="
	@echo ""
	@echo "  make docs-install         Install runtime + docs dependencies"
	@echo "  make docs-build           Build the MkDocs HTML site  -> $(SITE_DIR)/"
	@echo "  make docs-executed-nbs    Execute every docs/**/*.ipynb -> $(EXEC_DIR)/  (~30 min)"
	@echo "  make docs-check-nbs       Fail if any executed notebook contains errors"
	@echo "  make docs-jenner          Fast pipeline: use committed notebook outputs +"
	@echo "                            build + concat.  Output: $(COMBINED)"
	@echo "  make docs-jenner-execute  Full pipeline: execute notebooks + check + build +"
	@echo "                            concat.  Slower but validates notebooks still run"
	@echo "                            against the current API.  Output: $(COMBINED)"
	@echo "  make clean-docs           Remove $(SITE_DIR)/, $(EXEC_DIR)/, $(COMBINED)"
	@echo ""
	@echo "Tunable variables (override on the command line):"
	@echo "  PYTHON=$(PYTHON)"
	@echo "  SITE_DIR=$(SITE_DIR)"
	@echo "  EXEC_DIR=$(EXEC_DIR)"
	@echo "  COMBINED=$(COMBINED)"
	@echo "  NB_TIMEOUT=$(NB_TIMEOUT)        per-cell execution timeout (seconds; docs-jenner-execute)"
	@echo "  ALLOW_NB_ERRORS=$(ALLOW_NB_ERRORS)   set to 1 to publish combined doc even if notebooks errored"
	@echo "  NB_EXCLUDE=$(NB_EXCLUDE)   comma-separated substrings of paths to skip during execution"

# ── Install ───────────────────────────────────────────────────────────────────
# `uv venv` creates minimal venvs without pip; bootstrap it from the stdlib's
# bundled wheel if it's missing. No-op on standard `python -m venv` envs.
docs-install:
	@$(PYTHON) -c "import pip" 2>/dev/null || $(PYTHON) -m ensurepip --default-pip
	$(PYTHON) -m pip install -e .
	$(PYTHON) -m pip install -r docs/requirements.txt

# ── MkDocs HTML build ─────────────────────────────────────────────────────────
docs-build:
	$(PYTHON) -m mkdocs build --site-dir $(SITE_DIR)

# ── Notebook execution ────────────────────────────────────────────────────────
docs-executed-nbs:
	$(PYTHON) docs/execute_notebooks.py $(EXEC_DIR) --timeout $(NB_TIMEOUT) --exclude "$(NB_EXCLUDE)"

# ── Notebook error gate ───────────────────────────────────────────────────────
# Depends on docs-executed-nbs so that `make -j` (or `make docs-check-nbs` on
# its own) cannot scan an empty $(EXEC_DIR) and pass spuriously — the two are
# serialized through the prereq graph regardless of parallelism.
docs-check-nbs: docs-executed-nbs
	@$(PYTHON) docs/check_executed_nbs.py $(EXEC_DIR) \
	  $(if $(filter 1 true yes,$(ALLOW_NB_ERRORS)),--allow-errors,)

# ── Combined markdown pipeline (fast — uses committed notebook outputs) ──────
# The default doc-build path. Skips the ~30-minute notebook re-execution step
# and reads each notebook's *committed* outputs directly from docs/. That's
# sufficient for the RAG corpus: the outputs already in the .ipynb files
# capture every code cell's result, and re-executing them just to re-embed the
# same content in the corpus is wasted CI time. If a notebook committer
# forgets to re-run a cell after editing it, this fast build will not detect
# stale outputs; use `docs-check-nbs` / `docs-jenner-execute` to re-execute and
# fail on runtime errors against the current API.
docs-jenner: docs-build
	$(PYTHON) docs/concat_mkdocs.py $(SITE_DIR) docs $(COMBINED)

# ── Combined markdown pipeline with fresh notebook execution ─────────────────
# The historical full path: execute every notebook, error-gate the outputs,
# build the site, concat everything. Use when you want to validate that all
# notebooks still run cleanly against the current laser-generic API (typical
# for release-time builds) — the docs-check-nbs gate is what catches drift
# between notebook code and current APIs.
docs-jenner-execute: docs-check-nbs docs-build
	$(PYTHON) -c "from pathlib import Path; Path('$(COMBINED)').parent.mkdir(parents=True, exist_ok=True)"
	$(PYTHON) docs/concat_mkdocs.py $(SITE_DIR) $(EXEC_DIR) $(COMBINED)

# ── Clean ─────────────────────────────────────────────────────────────────────
# Shell-agnostic removal via Python so this works on cmd.exe / PowerShell / bash.
clean-docs:
	$(PYTHON) -c "import shutil; from pathlib import Path; [shutil.rmtree(p, ignore_errors=True) for p in ('$(SITE_DIR)', '$(EXEC_DIR)')]; Path('$(COMBINED)').unlink(missing_ok=True); print('Removed $(SITE_DIR)/, $(EXEC_DIR)/, $(COMBINED)')"

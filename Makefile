# laser-generic documentation pipeline.
#
# Lets you reproduce locally exactly what .github/workflows/build-combined-doc.yml
# does in CI. The "all-in-one" target is `docs-jenner`, which executes every
# notebook under docs/, builds the MkDocs site, then flattens both into a single
# combined_mkdocs.md suitable for RAG / MCP ingestion.

# All recipes are single-command invocations (echo or python script), so they
# work identically under bash, cmd.exe, and PowerShell — no SHELL override
# needed. Multi-step logic lives in docs/*.py helpers instead of inline shell.

.PHONY: help docs-install docs-build docs-executed-nbs docs-check-nbs docs-jenner docs-jenner-artifact clean-docs

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
	@echo "  make docs-install        Install runtime + docs dependencies"
	@echo "  make docs-build          Build the MkDocs HTML site  -> $(SITE_DIR)/"
	@echo "  make docs-executed-nbs   Execute every docs/**/*.ipynb -> $(EXEC_DIR)/"
	@echo "  make docs-check-nbs      Fail if any executed notebook contains errors"
	@echo "  make docs-jenner       Full pipeline (execute + check + build + concat)"
	@echo "                           Output: $(COMBINED)"
	@echo "  make docs-jenner-artifact  Build + concat only (assumes \$$(EXEC_DIR) pre-populated)"
	@echo "                           Used by the Execute Notebooks -> Build Combined Doc CI chain"
	@echo "  make clean-docs          Remove $(SITE_DIR)/, $(EXEC_DIR)/, $(COMBINED)"
	@echo ""
	@echo "Tunable variables (override on the command line):"
	@echo "  PYTHON=$(PYTHON)"
	@echo "  SITE_DIR=$(SITE_DIR)"
	@echo "  EXEC_DIR=$(EXEC_DIR)"
	@echo "  COMBINED=$(COMBINED)"
	@echo "  NB_TIMEOUT=$(NB_TIMEOUT)        per-cell execution timeout (seconds)"
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

# ── Full combined markdown pipeline ───────────────────────────────────────────
# docs-executed-nbs is reached transitively via docs-check-nbs. docs-build is
# independent of notebook execution (mkdocs-jupyter has execute:false and
# reads the source .ipynb files directly), so it's safe to keep in parallel.
docs-jenner: docs-check-nbs docs-build
	$(PYTHON) -c "from pathlib import Path; Path('$(COMBINED)').parent.mkdir(parents=True, exist_ok=True)"
	$(PYTHON) docs/concat_mkdocs.py $(SITE_DIR) $(EXEC_DIR) $(COMBINED)

# ── Combined markdown from a pre-built executed-notebooks tree ────────────────
# Same output as docs-jenner but skips both docs-executed-nbs and docs-check-nbs
# — assumes $(EXEC_DIR) is already populated (typically by the Execute
# Notebooks GitHub Action's artifact download). No execution, no error gate,
# both already guaranteed upstream. Local users generally want `docs-jenner`
# instead.
docs-jenner-artifact: docs-build
	$(PYTHON) -c "from pathlib import Path; Path('$(COMBINED)').parent.mkdir(parents=True, exist_ok=True)"
	$(PYTHON) docs/concat_mkdocs.py $(SITE_DIR) $(EXEC_DIR) $(COMBINED)

# ── Clean ─────────────────────────────────────────────────────────────────────
# Shell-agnostic removal via Python so this works on cmd.exe / PowerShell / bash.
clean-docs:
	$(PYTHON) -c "import shutil; from pathlib import Path; [shutil.rmtree(p, ignore_errors=True) for p in ('$(SITE_DIR)', '$(EXEC_DIR)')]; Path('$(COMBINED)').unlink(missing_ok=True); print('Removed $(SITE_DIR)/, $(EXEC_DIR)/, $(COMBINED)')"

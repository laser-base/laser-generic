# Research notebooks

Full-fidelity copies of tutorial notebooks that are too slow to execute on
every CI doc build. Each file here is the "research" counterpart to a
lighter CI-facing version at `docs/tutorials/notebooks/<name>.ipynb`.

## What "research" means here

The CI-facing notebook keeps the same pedagogical arc — imports, single-run
sanity check, parameter sweep, validation, regime robustness — but with
smaller `nticks`, fewer `nsims`, and/or trimmed sweep grids so it finishes
in a few minutes instead of tens of minutes. The research copy preserves the
original parameters: longer time series for cleaner peak-finding, more sims
for a smoother scatter plot, wider CBR sweep for stronger regime coverage.

Use the research copy when:

- You want the highest-fidelity numbers to write about or cite
- You're debugging a scientific issue and want more headroom before
  stochastic noise dominates
- You're validating a change to the underlying model against a stringent
  reference

Use the CI-facing copy when:

- You just want to see the plot and read the description
- You're running the doc build locally and want it done in a reasonable time

## How this folder is wired up

- `Makefile`: `NB_EXCLUDE` contains `research/`, so
  `docs/execute_notebooks.py` skips everything in this directory during the
  automated doc-build pipeline (`make docs-jenner`).
- `mkdocs.yml`: the `mkdocs-exclude` plugin's glob list contains
  `tutorials/notebooks/research/*`, so these notebooks don't appear in the
  rendered mkdocs site or navigation.

To run a research notebook locally: `jupyter nbconvert --to notebook
--execute docs/tutorials/notebooks/research/<name>.ipynb`, or just open it
in JupyterLab.

## Current inventory

- `06_SIR_wbirths_natural_periodicity.ipynb` — 10-simulation sweep at
  `nticks = 365*100` (100 y) with 6-CBR regime robustness sweep. The
  CI-facing version at
  `../06_SIR_wbirths_natural_periodicity.ipynb` uses `nsims = 5`,
  `nticks = 365*40` (40 y), and a trimmed 4-CBR regime sweep at
  `nticks = 365*40`. That change cuts nb06's slice of the doc build from
  ~16 min to ~3–4 min.

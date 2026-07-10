# Tutorial notebooks

This folder holds the Jupyter notebooks used as `laser-generic` tutorials. They cover single-patch and multi-patch SI / SIR / SEIR / SEIRS dynamics, births and deaths, seasonality, calibration, and case studies like the England-and-Wales measles pre-vaccine record.

## Where do I read them?

**Read the rendered tutorials on the docs site: [laser.idmod.org/laser-generic/tutorials/](https://laser.idmod.org/laser-generic/tutorials/).**

That's the freshly-executed, browsable, hosted-by-GitHub-Pages view. No download, no local Jupyter setup. Every notebook has its own page — for example:

- [SI model with no demographics](https://laser.idmod.org/laser-generic/tutorials/notebooks/01_SI_nobirths_logistic_growth/)
- [Intrinsic periodicity of the SIR system](https://laser.idmod.org/laser-generic/tutorials/notebooks/06_SIR_wbirths_natural_periodicity/)
- [Periodicity of measles in England and Wales](https://laser.idmod.org/laser-generic/tutorials/notebooks/10_England_and_Wales/)

The docs site is auto-rebuilt whenever notebook or library source changes — the [Execute Notebooks](../../../.github/workflows/execute-notebooks.yml) workflow executes every notebook, and the [Build Combined Doc](../../../.github/workflows/build-combined-doc.yml) workflow feeds those outputs into the site build. So the version you see there always corresponds to a specific committed `main` state.

If you want to read a notebook **inline on github.com** (as an `.ipynb` file rendered by GitHub's built-in notebook viewer), that also works, but the figures you see are whatever was committed with the notebook — they may not match the current source. See "Source vs outputs" below.

## Source vs outputs — what's authoritative?

**Source is authoritative. Rendered outputs live at [the docs site](https://laser.idmod.org/laser-generic/tutorials/). Committed `outputs` inside the `.ipynb` files are decorative.**

Concretely:

| Where the figures live | Authoritative? | How fresh? |
|---|---|---|
| Docs site (`laser.idmod.org/laser-generic/tutorials/`) | ✅ yes | Fresh on every push to `main` that affects notebook or library source |
| Committed `outputs` field in `docs/tutorials/notebooks/*.ipynb` | ❌ no — decorative, purely for github.com inline rendering | Whatever a contributor happened to commit |
| GitHub Actions `executed_nbs` artifact | ✅ yes (same source the docs site is built from) | Same as the docs site — one artifact per successful Execute Notebooks run |

Consequences:

- Contributors may commit executed notebooks (nice for github.com inline rendering) OR strip outputs (nice for lean diffs) — both work identically for the docs site.
- Committed outputs may drift from source over time; that's OK because the docs site never reads them.
- If you're demoing the model or citing a plot, link to the docs site, not to the github.com view of an `.ipynb` file.

## Downloading a specific commit's executed notebooks

**You don't need this for regular reading — go to the [docs site](https://laser.idmod.org/laser-generic/tutorials/) instead.**

This section is for the debugging / archival case: you need the raw `.ipynb` files (with outputs) that were executed against a specific commit — e.g. to reproduce a doc-quality regression, or to compare corpora across time.

1. Open the [Execute Notebooks workflow runs](https://github.com/laser-base/laser-generic/actions/workflows/execute-notebooks.yml).
2. Filter for the branch or commit you care about; the latest successful main-push run is [here](https://github.com/laser-base/laser-generic/actions/workflows/execute-notebooks.yml?query=is%3Asuccess+branch%3Amain).
3. Scroll to the *Artifacts* panel and click `executed_nbs` (retention: 90 days per repo-default; can be raised).
4. The zip contains every executed notebook plus a `manifest.json` recording provenance.

## Provenance (`manifest.json` fields)

Each `executed_nbs` artifact includes a `manifest.json` at its root with:

- `commit_sha`, `commit_ref` — the exact tree the outputs were built against.
- `source_hash` — the hash used as the cache key. Two artifacts with the same `source_hash` were built from byte-identical inputs.
- `python_version` — the Python that executed the cells.
- `run_id`, `run_number`, `workflow_url` — back-links to the CI run.
- `was_cache_hit`, `was_forced`, `was_allow_errors` — provenance flags.

If you're doing bisection or comparing corpora across time, `source_hash` and `commit_sha` are the two fields to trust.

## Regenerating locally (contributors only)

To execute the notebooks locally:

```
make docs-executed-nbs           # populates dist/executed_nbs/
make docs-check-nbs              # fail if any notebook errored
```

Or the full local doc pipeline:

```
make docs-jenner                 # execute + check + build + concat -> dist/combined_mkdocs.md
```

Set `GITHUB_ACTIONS=true` in the environment to force nb06's env-var-lite path (n_years=20, nsims=2) — matches what CI does. Without it, nb06 runs at full scale (n_years=100, nsims=10, ~15 min).

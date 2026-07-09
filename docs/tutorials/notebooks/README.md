# Tutorial notebooks

This folder holds the Jupyter notebooks used as `laser-generic` tutorials. They cover single-patch and multi-patch SI / SIR / SEIR / SEIRS dynamics, births and deaths, seasonality, calibration, and case studies like the England-and-Wales measles pre-vaccine record.

## Where do I read them?

- **Rendered on the docs site (recommended for browsing):** [laser.idmod.org/laser-generic → tutorials](https://laser.idmod.org/laser-generic/tutorials/).
- **Rendered inline on GitHub:** open any `*.ipynb` file in this folder. Note that committed outputs are decorative — they may be stale relative to the current source. See "Provenance" below.
- **Executed against a specific commit:** see the artifact section below.

## Source vs outputs — what's authoritative?

**Source is authoritative. Committed outputs are decorative.**

The tutorial notebooks are executed in CI by the [Execute Notebooks](../../../.github/workflows/execute-notebooks.yml) workflow. That workflow runs on every push to `main` that touches notebook or library source, publishes an `executed_nbs` artifact, and feeds that artifact into the doc-site build. The rendered figures you see at `laser.idmod.org/laser-generic` come from that artifact, not from the `outputs` field of the committed notebook files.

That means:
- Contributors may commit executed notebooks (nice for GitHub inline rendering) OR strip outputs (nice for lean diffs) — both work identically for the doc build.
- Committed outputs may drift from source over time. That's OK, because no downstream product reads them.
- The source-of-truth executed notebooks live in the CI artifact.

## Downloading the executed notebooks

If you want a specific commit's freshly-executed notebooks (e.g. to debug a doc-quality regression, or to preview what the docs will render against a change), grab them from the workflow's artifact panel:

1. Open the [Execute Notebooks workflow runs](https://github.com/laser-base/laser-generic/actions/workflows/execute-notebooks.yml).
2. Filter for the branch or commit you care about; the latest successful main-push run is [here](https://github.com/laser-base/laser-generic/actions/workflows/execute-notebooks.yml?query=is%3Asuccess+branch%3Amain).
3. Scroll to the *Artifacts* panel and click `executed_nbs`.
4. Retention is 400 days.

## Provenance

Each `executed_nbs` artifact includes a `manifest.json` at its root with:

- `commit_sha`, `commit_ref` — the exact tree the outputs were built against.
- `source_hash` — the hash used as the cache key. Two artifacts with the same `source_hash` were built from byte-identical inputs.
- `python_version` — the Python that executed the cells.
- `run_id`, `run_number`, `workflow_url` — back-links to the CI run.
- `was_cache_hit`, `was_forced`, `was_allow_errors` — provenance flags.

If you're doing bisection or comparing corpora across time, the `source_hash` and `commit_sha` are the two fields to trust.

## Regenerating locally

To execute the notebooks locally against a specific check-out:

```
make docs-executed-nbs           # populates dist/executed_nbs/
make docs-check-nbs              # fail if any notebook errored
```

Or run the full local doc pipeline:

```
make docs-jenner-execute         # (available once #237 lands) execute + check + build + concat
```

Set `GITHUB_ACTIONS=true` in the environment to force `nb06`'s env-var-lite path (n_years=20, nsims=2) — matches what CI does. Without it, `nb06` runs at full scale (n_years=100, nsims=10, ~15 min).

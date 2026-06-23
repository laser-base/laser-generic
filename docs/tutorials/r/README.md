# laser-generic from R (via reticulate)

Tutorials for driving the `laser.generic` Python package from R, using
[`reticulate`](https://rstudio.github.io/reticulate/).

## Contents

| File                   | What it covers                                                                                                                                                                                        |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `01-models.Rmd`        | Each of the six standard disease models: SI, SIS, SIR, SIRS, SEIR, SEIRS. One self-contained section per model with a runnable example and a plotted trajectory.                                      |
| `02-customization.Rmd` | Customizing an SEIR base model with the four vital-dynamics components: `ConstantPopVitalDynamics`, `BirthsByCBR`, `MortalityByCDR`, `MortalityByEstimator`.                                          |
| `03-spatial.Rmd`       | Multi-node scenarios and the migration network: a two-patch coupling sweep, a 3×3 grid with the default gravity network, and replacing `model$network` with custom topologies (none, chain, uniform). |
| `helpers.R`            | Shared R helpers (scenario factory, compartment-trajectory extractor, per-node extractor, ggplot helpers). Sourced by every `.Rmd` file.                                                              |

## Requirements

- **R ≥ 4.1** (uses the native `|>` pipe).
- The R packages below. Installed once from CRAN:

    ```r
    install.packages(c("reticulate", "rmarkdown", "ggplot2", "dplyr", "tidyr", "knitr"))
    ```

- The first time you knit a document, `reticulate` will provision an
  ephemeral Python environment (via `py_require()`) that pulls
  `laser-generic` from PyPI along with its dependencies. Subsequent
  knits reuse the cached environment and start instantly.

## Knitting

From R (or RStudio):

```r
# Render to BOTH html_document and github_document in one call
rmarkdown::render("01-models.Rmd",        output_format = "all")
rmarkdown::render("02-customization.Rmd", output_format = "all")
rmarkdown::render("03-spatial.Rmd",       output_format = "all")
```

From the shell:

```sh
for f in 01-models 02-customization 03-spatial; do
    Rscript -e "rmarkdown::render('${f}.Rmd', output_format = 'all')"
done
```

Each render produces two artifacts:

- `*.html` — the rich format with floating TOC; what you'd open locally.
- `*.md` — GitHub-friendly markdown that renders inline on the repo page.

## Notes

- The tutorials prefer single-node scenarios for clarity. The same
  patterns scale to grids by passing `M > 1` and/or `N > 1` to
  `make_scenario()` in `helpers.R`.
- All randomized output is reproducible: each chunk seeds both R and
  `laser.core.random.seed(...)` before constructing the model.

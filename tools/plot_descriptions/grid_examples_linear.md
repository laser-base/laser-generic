#### Reading the linear (1-D) grid plot

A single row of 9 patches stretched horizontally — x-axis runs 0 to 90,000 meters, y-axis runs 0 to 10,000 meters, so each patch is 10 km square but laid out edge-to-edge along a line. Red `x` centroids sit at (5000, 15000, …, 85000) on a single row at y = 5000. The `viridis` colorbar spans 0.2 x 10^6 to 1.0 x 10^6.

The central patch (centroid at 45,000 m) is bright yellow at 1,000,000 — the maximum. Color falls off symmetrically to either side following `linear_pop` (`1_000_000 / (|col - CX| + 1)`): the immediate neighbors are teal at ~500,000, then ~333,000, ~250,000, and the two endpoints are darkest purple at ~200,000. **The figure demonstrates that `grid()` with `ROWS = 1` collapses to a one-dimensional chain of patches, useful for stylized linear-corridor scenarios where the `population_fn` only needs to depend on the column index.**

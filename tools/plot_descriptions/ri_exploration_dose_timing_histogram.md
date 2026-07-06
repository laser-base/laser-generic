**Reading the dose-timing histogram.**
The `RoutineImmunizationEx` component was constructed with `track=True`, exposing the per-agent `initial_ri` array — the day-of-life at which each individual is scheduled to receive their RI dose. This plot keeps only non-zero entries (i.e. agents actually scheduled, not the 30% un-covered tail).

- **Sky-blue histogram** (unlabeled): non-zero `initial_ri` values, 30 bins spanning roughly 210–340 days. The shape is a clean bell with peak counts of ~2,000 near the center bins.
- **"Mean 273.5"** (red, dashed vertical line): at the empirical mean of the non-zero `initial_ri` values.
- **"Normal(274, 15) PDF"** (dark-blue, solid curve): the analytic PDF, rescaled by `nonzero_initial_ri.size * bin_width` so its peak matches the histogram. It hugs the bar tops across the full range.

The empirical mean of 273.5 days lines up almost exactly with the requested location parameter of 274, and the bell's width matches the requested $\sigma = 15$. **The takeaway: the `dose_timing_dist = distributions.normal(loc=274, scale=15)` argument flows through `RoutineImmunizationEx` faithfully — the per-agent dose ages drawn at birth reproduce the requested normal distribution, confirming that the 9-month-mean / 15-day-spread vaccination schedule is being applied as configured.**

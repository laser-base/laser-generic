#### Reading the constant-CBR population-growth plot

A linear-axis plot of total population versus time in days (0 to 3650, i.e. ten years), with a faint background grid:

- **"Simulated Total Population"** (blue, solid): $S+E+I+R$ aggregated across the single node.
- **"Expected Growth (CBR)"** (orange, dashed): the analytic compound-interest formula $N_0(1 + \text{CBR}/1000)^{t/365}$ with $\text{CBR}=20$ per 1000 per year.

The trajectory starts at the seed population of 59,116 and climbs smoothly to roughly 72,000 over the ten-year window — almost exactly a $(1.02)^{10}\approx 1.22\times$ expansion, matching the printed actual-vs-expected difference of just $-2$ agents ($-0.00\%$). The blue and orange curves are visually indistinguishable across the entire range; you can only see the simulated line at all because it is drawn under the dashed overlay. The slope steepens slightly as the population compounds, giving the curve its characteristic gentle exponential bow rather than a straight line. **The takeaway: `BirthsByCBR` driven by a scalar `ValuesMap.from_scalar` reproduces continuous-compounding demographic growth to within rounding error, validating the births component against the closed-form CBR formula in the simplest single-node, time-constant case.**

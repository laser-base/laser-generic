### Reading the SI baseline (no seasonality) plot

A single-node SI run on 1,000,000 agents over 730 days, $\beta = 0.03125$, seeded with 10 infections. Single linear y-axis (population, 0 to 1e6). Two curves:

- **Blue — Susceptible.** Starts at 1,000,000, holds essentially flat through day ~200, then rolls over through an inflection near day ~370 and asymptotes to 0 by day ~600.
- **Red — Infectious.** The mirror image: near zero through day ~200, rising steeply through the same inflection, and saturating at 1,000,000 by day ~600.

The two curves cross at exactly 500,000 near day ~370 — the classic SI logistic midpoint where $S = I = N/2$ and $\dot{I}$ is maximal. There is no recovery compartment, so every agent eventually moves to I and stays there.

**The plot demonstrates the canonical undriven SI take-off-and-saturate curve at the chosen $\beta$ — the control trajectory against which the next cell's "radical seasonality" (an intervention applied between days 300 and 400) will be compared.**

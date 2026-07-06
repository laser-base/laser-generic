**Reading the SI "radical seasonality" intervention plot.**
Same SI configuration as the previous baseline — 1,000,000 agents, $\beta = 0.03125$, 730 days, linear y-axis 0 to 1e6 — but now the seasonality array is used as a non-cyclical intervention: $\beta(t)$ is multiplied by 1.0 through day 300, ramps linearly down to 0.0 between days 300 and 400, and stays at 0.0 thereafter.

The blue Susceptible curve is identical to the baseline through day ~300 (flat near 1,000,000), begins falling at the same inflection near day ~310, but instead of plunging to zero it flattens out at ~710,000 by day ~400 and stays there. The red Infectious curve mirrors it: rising from zero near day ~250, climbing through the ramp-down window, and plateauing at ~290,000 by day ~400. Roughly 29% of the population is ever infected before the outbreak halts.

**The plot demonstrates that the seasonality input is general-purpose temporal modulation of $\beta(t)$ — not restricted to cyclical forcing — and can be used to model the population-level effect of a time-limited intervention that arrests an outbreak partway through its logistic take-off.**

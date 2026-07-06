**Reading the SI reference suite plots.**
Two vertically stacked linear-scale subplots over 0–365 days, ten overlapping seed-trajectories visible as a single ribbon at $10^6$ agents each. Same $\beta = 0.3$ as the SEI run, but with no exposed compartment.

- **Top — Susceptible** $S_t$. Stays near $10^6$ only briefly (through day ~15), then sweeps down sigmoidally and hits ~0 by day ~70.
- **Bottom — Infectious** $I_t$. The complementary sigmoid: near-zero seed, exponential take-off starting around day 15, and saturation at $10^6$ by day ~70.

Compared with the SEI suite in the previous cell, the SI epidemic reaches half-saturation roughly 40 days earlier and finishes about 30 days sooner. There is no middle subplot because the SI model has no latent class — newly infected agents become immediately infectious. **The figure demonstrates the SI baseline against which the SEI delay can be measured: with the same per-tick transmission rate $\beta$, removing the incubation period strips out the lag and shifts the entire epidemic curve to the left.**

### Reading the SEI single-node suite plots

Three vertically stacked linear-scale subplots sharing a 0–365 day x-axis, one curve per seed (10 seeds), all overlapping into a single visual ribbon — the ten realizations are effectively indistinguishable at this population size of $10^6$ agents per simulation.

- **Top — Susceptible** $S_t$. Holds flat at $10^6$ from day 0 through day ~40, then drops in a sigmoidal sweep, reaching near zero around day 100 and remaining flat thereafter.
- **Middle — Exposed** $E_t$. A single sharp pulse: nearly zero until day ~50, climbing to a peak of ~215,000 around day 80, then decaying back to zero by day ~120. This is the transient pool of agents in the latent period waiting for their `etimer` to expire.
- **Bottom — Infectious** $I_t$. Mirror image of $S_t$: near zero through day ~50, sigmoidal climb, and saturating at the full $10^6$ population by day ~100 where it stays for the rest of the run.

The exposed peak coincides with the steepest part of the $S \to I$ transition — exactly when transmission is at maximum throughput. **The figure demonstrates that the new `ExposedSEI` component produces the textbook SEI signature: a transient exposed bump sandwiched between a depleting susceptible pool and a saturating infectious absorbing state.**

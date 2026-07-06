#### Reading the SEIS single-node suite plots

Three vertically stacked linear-scale subplots, 0–365 days, with all ten seed realizations overlapping into a single visual ribbon at $10^6$ agents per simulation.

- **Top — Susceptible** $S_t$. Holds at $10^6$ through day ~50, drops sigmoidally to a minimum around day ~100 of ~200,000, then rebounds slightly and locks onto an equilibrium plateau near ~215,000 for the rest of the run.
- **Middle — Exposed** $E_t$. Climbs from zero starting around day 50, overshoots to a peak near ~225,000 around day 90, performs a small damped oscillation, and settles to a steady-state pool of ~205,000.
- **Bottom — Infectious** $I_t$. Climbs in parallel, peaks near ~600,000 around day 100, oscillates once, and converges to an endemic equilibrium near ~580,000.

The damped oscillation visible in all three panels around days 90–130 is the signature of the SEIS system relaxing toward its endemic fixed point — agents flowing $S \to E \to I \to S$ no longer produce a one-shot epidemic. **The figure demonstrates that the new `InfectiousSEIS` component, by recycling recovered agents back to $S$, transforms the one-shot SEI outbreak into an endemic equilibrium where the three compartments settle at fixed fractions of the population determined by the ratio of exposure duration, infectious duration, and $\beta$.**

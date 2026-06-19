### Reading the baseline (no seasonality) plot

A 10-year SIR trajectory on a single node of 1,000,000 agents, time on the x-axis from 0 to ~3,650 days. Two y-axes share the figure:

- **Left axis (0 to 1e6)** carries blue (Susceptible) and green (Recovered) lines.
- **Right axis (0 to ~3,000)** carries the red Infectious curve and a flat gray "Seasonality (scaled)" reference line that sits dead-level at the midpoint — confirming $\beta(t)$ has no temporal modulation in this run.

S sits near 70,000 and R near 940,000 throughout, both nearly flat. The red I curve, by contrast, oscillates wildly with peaks of ~2,500 spaced roughly every 500 days — an emergent, undriven cycle arising from the interplay between waning susceptibility (births via `ConstantPopVitalDynamics` with CBR=33) and depletion-driven outbreak crashes. Peak heights drift between ~2,200 and ~2,500 with no strict periodicity.

**The plot demonstrates that even without seasonal forcing, an SIR model with births produces damped, irregular endemic oscillations — the natural baseline that subsequent seasonal-forcing experiments must be compared against.**

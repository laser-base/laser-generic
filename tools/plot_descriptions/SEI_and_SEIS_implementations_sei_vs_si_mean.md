#### Reading the SEI vs SI mean-trajectory plot

A single linear-scale panel, 0–365 days on the x-axis, 0 to $10^6$ agents on the y-axis. Four curves showing seed-averaged means:

- **"SEI Mean Susceptible"** (blue, solid): flat at $10^6$ through day ~40, sigmoidal drop, crossing $5 \times 10^5$ near day ~75, and ~0 by day ~100.
- **"SEI Mean Infectious"** (orange, solid): the mirror sigmoid, crossing $5 \times 10^5$ near day ~75 and saturating at $10^6$ by day ~100.
- **"SI Mean Susceptible"** (blue, dashed): same shape but shifted left — crosses $5 \times 10^5$ near day ~40 and bottoms out by day ~60.
- **"SI Mean Infectious"** (orange, dashed): mirror of the dashed-blue, saturating at $10^6$ by day ~60.

The blue/orange crossover for SI sits near day 40 while the SEI crossover sits near day 75 — a roughly 35-day lag attributable entirely to the ~5-day mean latent period being compounded across many generations of transmission. **The figure demonstrates that adding an exposed compartment with mean incubation of 5 days delays the full SEI epidemic by tens of days relative to an SI model with identical $\beta$, even though the final attack rate is identical (both saturate at $N$).**

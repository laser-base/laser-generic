### Reading the MortalityByCDR single-run SEIR channels plot

Four colored line traces on a linear-axis figure of population vs. time in days (0 to ~3650, i.e. 10 years). Initial compartments are seeded at S=37,500, E=12,500, I=12,500, R=37,500 out of a 100,000-agent single-node population.

- **"Susceptible (S)"** (blue): starts at ~37,500 and decays roughly linearly to ~30,500 by day 3650.
- **"Exposed (E)"** (orange) and **"Infectious (I)"** (red): both spike briefly in the first ~80 days as the seeded exposures progress through the SEIR pipeline (red "I" peaks near 14,000), then collapse to zero once the initial cohort flushes through and the susceptible pool is no longer being drawn down by transmission.
- **"Recovered (R)"** (green): climbs sharply from 37,500 to a peak near 62,000 by day ~80 as the initial E/I cohort recovers, then declines roughly linearly to ~51,000 by day 3650.

The post-day-100 dynamics are pure non-disease mortality: with no births and a constant crude death rate (CDR = 20 per 1,000 per year), each compartment decays at a common per-capita rate while the relative S:R proportions are preserved. **This figure demonstrates that `MortalityByCDR` removes agents from every SEIR compartment at the configured crude death rate without distorting the relative compartment composition, producing the characteristic linear-on-this-scale exponential decay of a closed population under constant mortality.**

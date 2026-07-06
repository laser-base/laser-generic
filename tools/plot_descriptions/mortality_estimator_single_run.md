#### Reading the MortalityByEstimator SEIR channels plot

Same four-line layout as the earlier `MortalityByCDR` single-run plot — population vs. time in days (0 to ~3650), four colored traces for S, E, I, R — but the mortality engine is now `MortalityByEstimator`, which draws a date-of-death (DoD) per agent from a Kaplan-Meier survival curve fit to Nigeria 2020 life-table data, with initial ages sampled from a Nigeria 2020 population pyramid via `AliasedDistribution`.

- **"Susceptible (S)"** (blue): starts at ~37,500 and decays gently and roughly linearly to ~33,500 by day 3650, a shallower slope than under the constant-CDR run.
- **"Exposed (E)"** (orange) and **"Infectious (I)"** (red): same transient SEIR pipeline behavior: a brief spike in the first ~80 days (red "I" peaks near 14,000), then both fall to zero.
- **"Recovered (R)"** (green): rises to a peak near 62,000 by day ~80, then declines to ~55,800 by day 3650 — again, less steeply than the constant-CDR analog.

The whole-population observed mortality printed below the plot is ~11.5 per 1,000 per year, well below the CDR=20 case, reflecting Nigeria's age-specific survival profile applied to the seeded age distribution. **This figure demonstrates that `MortalityByEstimator` produces a closed-population decay shaped by an age-at-death distribution drawn from a Kaplan-Meier survival curve and a population pyramid, rather than a single crude rate — yielding age-structured non-disease mortality whose realized annual rate emerges from the underlying life table and pyramid.**

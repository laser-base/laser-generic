#### Reading the 70%-coverage routine-immunization SEIR plot

Same single-node SEIR setup as the baseline, but now with `RoutineImmunizationEx` added: 70% coverage (`distributions.constant_float(0.7)`) and a dose timing drawn from `Normal(loc=274, scale=15)` days (about 9 months) with a 6-month floor.

Four curves on twin y-axes over 730 days:

- **Blue solid — Susceptible** (left axis): starts near 31,000 and drifts *upward* throughout the run, ending around 36,900 (slightly higher than start as births minus RI minus infections net positive in the final stretch). The slope is much shallower than the baseline because RI is continuously draining new births out of S into R.
- **Green solid — Recovered** (left axis): climbs from ~169,000 to ~170,500 — almost flat, since most over-fives were already R.
- **Red solid — Infectious** (right axis): noisy outbreak peaking near 180 around day 350.
- **Orange solid — Exposed** (right axis): tracks I, peaking near 130.

Compared to the baseline plot's peak of ~410 infectious individuals, the I peak here is less than half. **The takeaway: routine immunization at 70% coverage with a 9-month mean dose age preserves the outbreak shape but cuts the peak roughly in half and drops cumulative infections from ~14,000 to ~6,500.**

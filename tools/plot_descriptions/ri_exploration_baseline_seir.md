#### Reading the no-RI baseline SEIR plot

The baseline run is a single-node SEIR over 730 days with a 200,000 starting population, 100 seed infections, and $R_0 = 7$. After `initialize_susceptibility` reclassifies everyone over five as recovered, only ~31,000 individuals are left in S — so this is effectively an outbreak in the under-five compartment.

Four curves share the time axis (0–730 days) on twin y-axes:

- **Blue solid — Susceptible** (left axis): starts near 31,000, drifts gently down to ~27,000 around day 500, then climbs back toward ~28,800 as births refill the under-five pool.
- **Green solid — Recovered** (left axis): rises monotonically from ~169,000 to ~178,000.
- **Red solid — Infectious** (right axis): noisy outbreak rising from 100, peaking near 410 around day 350, then decaying to ~30 by day 730.
- **Orange solid — Exposed** (right axis): tracks I roughly 7 days earlier (the exposure mean), peaking near 300.

**The takeaway: in the no-RI baseline with $R_0 = 7$ acting on the ~31k under-five susceptibles, the outbreak runs its full course over ~700 days, ultimately recruiting roughly 14,000 cumulative infections (printed below the figure).**

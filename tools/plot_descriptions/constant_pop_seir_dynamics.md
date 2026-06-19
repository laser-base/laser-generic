### Reading the constant-population SEIR plot

A single wide panel with twin y-axes versus time in days (0 to ~3650, i.e. 10 years), six curves drawn from a 25-node SEIR run with `ConstantPopVitalDynamics` and a CBR of 35 per 1000 per year:

- **Left y-axis (Population, 0 to $1.4\times10^7$)** — **black solid** Total Population pinned flat across the top; **blue solid** Susceptible (S) flat at ~$1.4\times10^6$; **green solid** Recovered (R) flat at ~$1.28\times10^7$ after a brief settle from the $1/R_0$ initialization.
- **Right y-axis (Counts, 0 to ~11000)** — **red dashed** Infectious (I) and **orange dashed** Exposed (E), both noisy oscillations with the I curve riding ~3000 above the E curve; **purple dotted** Births and **brown dotted** Deaths sitting on top of each other at ~1500/day.

I and E show damped multi-annual waves with peaks around day 500, 1200, 2100, 2700, and 3200, never settling to a smooth endemic plateau but staying within a roughly $\pm 25\%$ envelope around the equilibrium prevalence of 9000/12M used to seed I. **Births and deaths trace the same horizontal band**, which is the point of `ConstantPopVitalDynamics`: the recycle process replaces each death with a new susceptible birth so $N = S + E + I + R$ never drifts. **The figure demonstrates that under constant-population vital dynamics the SEIR model sustains endemic transmission with stochastic multi-year cycles while births exactly balance deaths and total population stays fixed.**

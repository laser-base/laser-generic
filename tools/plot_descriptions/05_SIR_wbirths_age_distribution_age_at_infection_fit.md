### Reading the age-at-infection fit plot

A density-normalized histogram (blue) of simulated ages at infection (in years) for infections occurring after day 5000, overlaid with two curves:

- **"Expected exponential distribution - A = 0.98 years"** (orange, solid, thick): the analytic prediction $A = \frac{1}{R_0 \mu}$ once mortality censoring is folded in.
- **"Best fit age of infection, A = 0.98 years"** (black, dashed): the exponential fit to the simulated ages at infection.

The histogram itself is labeled **"Ages from simulation"** in the legend.

The x-axis spans 0 to 15 years, y-axis 0 to ~1.03. The histogram peaks at density ~1.03 in the 0-0.5-year bin and decays monotonically, reaching the floor by ~6 years. Both fit curves visually overlay the histogram across the entire range, and the analytic and fitted means agree to two decimal places. The KS statistic returned alongside the plot quantifies the goodness of fit.

**The figure demonstrates that, at endemic equilibrium, the simulated age at infection follows an exponential distribution with mean $\frac{1}{R_0 \mu}$ — i.e. the observed mean age of infection, properly censored by background mortality, matches Keeling/Rohani's textbook result.**

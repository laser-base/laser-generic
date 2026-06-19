### Reading the $\bar{A}_I$ expected-vs-observed plot

A scatter plot from a 25-run sweep over random pairs $(R_0 \in [5, 15], \text{cbr} \in [70, 100])$:

- **X-axis — Average_Iage_expected**, the analytic $\frac{1}{R_0 \mu}$ converted to years, ranging from ~0.73 to ~2.33.
- **Y-axis — Average_Iage_observed**, the mean of an exponential fit to simulated ages at infection (after day $365 \times 40$), spanning the same range.
- **Blue dots** — one per simulation; **red dashed** — the identity line.

All 25 points sit tightly on the identity line across the full 0.73-2.33 year range. The notebook layout reserves space for two additional panels (susceptibility-age and population-age expected-vs-observed) which are commented out and therefore appear empty below the populated top panel. The accompanying printout shows the average fractional deviation, max deviation, and counts of runs exceeding 5% and 10% error — the verbal pass/fail criteria for the test.

**The figure demonstrates that the predicted mean age at infection $\frac{1}{R_0 \mu}$ is recovered to within a few percent across a wide sweep of $R_0$ and birth-rate values, validating the SIR-with-demography implementation as a scientific instrument for endemic-equilibrium analysis.**

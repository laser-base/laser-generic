### Reading the CCS calibration-grid plot

A 6-by-2 grid of small-multiple panels, one per parameter set, showing the **critical community size (CCS)** signature: proportion of weeks with zero measles cases versus $\log_{10}(\text{population})$ across the England-and-Wales places.

Each panel overlays:

- **"Observed"** (blue dots): observed proportion zero per place (from the historical case records).
- **"Simulated"** (orange dots): simulated proportion zero from that sim's last 20 years.
- **"Obs fit"** (blue, solid): logistic fit to observed, transitioning from $\sim 1$ at $\log_{10}(\text{pop}) \approx 3$ down to $\sim 0$ near $\log_{10}(\text{pop}) \approx 5.5$.
- **"Sim fit"** (orange, solid): logistic fit to simulated.

Panels are sorted by the `similarity_CCS` score in the title (0.789 best, top-left, up to 1.200 in this top-12 slice). Panel titles also report the calibration draws ($\beta$, seasonal amplitude, gravity $k$, $b$, $c$). The best panels have the orange and blue logistic curves nearly overlying through the transition region around $\log_{10}(\text{pop}) \approx 4$–5; worse panels show the simulated curve dropping more steeply or shifted left of the observed.

**This figure demonstrates which corners of the gravity-and-seasonality parameter space reproduce the population-size-dependent fadeout structure that defines measles CCS in pre-vaccine England and Wales.**

#### Reading the best-combined-rank plots

The cell prints the best simulation index and parameters, then emits two figures side by side in the rendered output. Both panels share the same parameter-set title line: `beta=3.87, amplitude=1.31, k=0.0446, b=0.251, c=1.56`.

- **Top — CCS Plot.** Proportion of weeks with zero cases versus $\log_{10}(\text{population})$, x-axis 2.8 to 6.5. Four legend entries: **"Observed"** (blue dots) — the observed places; **"Simulated"** (orange dots) — simulated places from this best sim; **"Obs fit"** (blue, solid) and **"Sim fit"** (orange, solid) — logistic curves descending from $\sim 1$ at low population to $\sim 0$ near $\log_{10}(\text{pop}) = 6$. The simulated cloud sits very close to the observed cloud through the transition region around $\log_{10}(\text{pop}) \approx 4$–5, with the simulated logistic dropping slightly more steeply than the observed one between $\log_{10}(\text{pop}) = 4.5$ and 5.5.

- **Bottom — Wavelet Phase Plot.** Phase difference (degrees) versus distance from London, x-axis 0 to 30, y-axis $\sim -55$ to $+5$. Two legend entries: **"Observed"** (black dots) — observed phase lags; **"Simulated"** (blue dots) — simulated lags for this same best sim. Both clouds show the characteristic downward trend, with the blue simulated points interleaving with the black observed points across the full distance range.

**This figure demonstrates that a single parameter set ($\beta \approx 3.87$, seasonal amplitude $\approx 1.31$, gravity $k \approx 0.045$, $b \approx 0.25$, $c \approx 1.56$) simultaneously reproduces the CCS fadeout signature and the travelling-wave phase-lag signature in pre-vaccine England-and-Wales measles — the headline calibration result of the notebook.**

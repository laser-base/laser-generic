### Reading the CCS-vs-alpha log-log fit plot

A single log-log scatter overlaying the fitted power law on the simulation cloud:

- **x-axis — $\alpha$** on a log scale from ~50 to ~2000.
- **y-axis — CCS_est** on a log scale from ~$10^3$ to ~$10^6$.
- **"Simulations"** (filled markers, colored by $R_0$ on a viridis scale — dark ~2, yellow ~15): the 200 sims (after filtering the floor/ceiling cases).
- **"Fit R0=1.5"**, **"Fit R0=3"**, **"Fit R0=8"**, **"Fit R0=16"** (four straight diagonal lines): the fitted relationship evaluated at each of those $R_0$ values. All four have the same slope (~1.83, the recovered $\alpha$ exponent) and are simply offset vertically by the $(R_0/(R_0-1))^{1.67}$ factor — the blue $R_0=1.5$ line sits highest, then orange ($R_0=3$), then green ($R_0=8$) and red ($R_0=16$) nearly on top of each other near the bottom.

The coloured points slot into the corresponding bands — dark (low-$R_0$) points cluster above the orange line, yellow (high-$R_0$) points hug the red line, and the trend across nearly three decades of $\alpha$ is clean and roughly straight. **The figure demonstrates that the CCS scales as a near-$\alpha^{3/2}$ power law and that increasing $R_0$ lowers CCS only modestly once $R_0$ is comfortably above 1 — confirming the weak-$R_0$, strong-$\alpha$ behaviour predicted by the Nasell-style analytic formulas.**

#### Reading the SIR-trajectory and age-at-infection plots

The first cell emits two figures from a single 30-year simulation of an SIR model with births and deaths at a crude birth rate of 90 per 1000:

- **Top — compartment trajectories** versus tick (0 to ~11000 days). Four curves labeled by compartment: **"S"** (blue) crashes from ~300000 down through a deep oscillation, **"I"** (red) spikes briefly, **"R"** (green) shoots up to ~275000, and **"N"** (black) plateaus near 290000. After damped oscillations across the first ~2000 days, all four lines settle to a flat endemic equilibrium: $S^* \approx 25000$, $I^* \approx 3000$, $R^* \approx 258000$.
- **Bottom — histogram of age at infection** (in days, `doi - dob`) restricted to infections after tick $365 \times 28$. The shape is a clean monotonic exponential decay, peaking at ~11000 counts in the 0-100 day bin and tailing off into negligible counts past ~2000 days.

The trajectory plot confirms $N = S + I + R$ is preserved while the system reaches its predicted endemic equilibrium; the histogram is a first visual hint that age at infection at equilibrium is exponentially distributed. **Together the figures demonstrate that the SIR-with-demography model settles into the analytic endemic equilibrium and produces an exponentially distributed age at infection driven by the constant force-of-infection hazard.**

#### Reading the top-10 wavelet-phase-similarity plots

Ten single-axes scatter plots, one per figure, each emitted by the same loop and titled `Simulation <idx> Wavelet Phase Similarity`. Across all ten the axes are identical: x-axis is distance from London (0 to ~30), y-axis is phase difference in degrees (roughly $-90$ up to $+5$).

Two overlaid point clouds in every panel:

- **"Observed"** (black, filled): phase lags from the previous London-wavelet plot, the same downward-trending cloud from $\sim 0^\circ$ near London out to $\sim -45^\circ$ at distance 30.
- **"Sim {idx} (sim=...)"** (light-blue dots): phase lags for the named simulation; the title's `sim=` value reports that simulation's squared-difference score against the observed cloud.

Simulations are ranked by ascending similarity (best match first); scores in this top-10 slice range from $\sim 8093$ (Sim 33) up to $\sim 14803$ (Sim 121). In the better-ranked panels the blue cloud's central tendency closely tracks the black observed cloud's downward slope, with some sims showing a few outlier points reaching $-60^\circ$ to $-80^\circ$. Lower-ranked panels (e.g. Sim 46, Sim 121) display more diffuse blue clouds with weaker downward trend.

**Collectively these panels demonstrate that a subset of the 150 calibration draws reproduces the travelling-wave phase-lag signature with distance from London — and the wavelet-phase similarity metric correctly ranks them by how tightly the simulated lag-vs-distance cloud overlays the observed one.**

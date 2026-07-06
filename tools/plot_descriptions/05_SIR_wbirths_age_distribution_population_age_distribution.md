**Reading the population-age-distribution plot.**
Density-normalized histogram of agent ages at the final tick (in years, 0 to 60 on the x-axis, 0 to ~0.09 on the y-axis), overlaid with:

- **"Expected exponential distribution - A = 11.60 years"** (orange, solid, thick): the analytic prediction with mean $A = \frac{1}{\mu}$, where $\mu$ is the per-day mortality rate matching the 90-per-1000 crude birth rate.
- **"Best fit age distribution, A = 11.47 years"** (black, dashed): the exponential fit to the simulated age distribution.

The histogram itself is labeled **"Simulation output"** in the legend.

The histogram peaks near density 0.087 in the youngest bin and decays smoothly to the noise floor by ~50 years; the two overlaid curves visually overlay the histogram and each other across the entire range. The analytic and fitted means agree to within ~1%.

This is a demographic check rather than an epidemiological one: at constant equal birth and death rates, the stable population age distribution should be exponential with rate $\mu$. **The figure demonstrates that the `BirthsByCBR` plus `MortalityByEstimator` components produce the textbook stable exponential age distribution with mean $\frac{1}{\mu}$, confirming the demographic plumbing underneath the epidemiology is correct.**

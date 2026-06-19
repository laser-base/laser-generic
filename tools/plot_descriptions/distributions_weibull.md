### Reading the Weibull plot

Overlaid density histograms on $[0, 2.5]$ (samples clipped at 2.5), y-axis density from 0 to about 7.3. Four shape values at fixed scale $\lambda = 1$, sampled at $N = 100{,}000$:

- **Blue — Weibull(k = 0.5, λ = 1)** monotone decay with a huge spike at the origin reaching density ~7.3, the sub-exponential regime ($k < 1$) where the hazard rate decreases over time.
- **Red — Weibull(k = 1.0, λ = 1)** pure exponential decay (Weibull reduces to Exponential when $k = 1$), intercept ~1.0 at the origin.
- **Purple — Weibull(k = 1.5, λ = 1)** unimodal with mode near $x \approx 0.5$ at density ~0.7, the increasing-hazard regime.
- **Green — Weibull(k = 5.0, λ = 1)** sharply peaked symmetric hump centered near $x \approx 0.95$ with peak density ~1.9 — high shape values make Weibull look almost Gaussian around the scale parameter.

Each curve passes through the classic Weibull transition from infinite-density-at-zero through exponential to peaked. **The figure demonstrates that `distributions.weibull` reproduces the full Weibull family across decreasing, constant, and increasing hazard regimes, with the mode and tail behavior controlled by the shape parameter $k$.**

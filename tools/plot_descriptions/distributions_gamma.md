#### Reading the Gamma plot

Overlaid density histograms on $[0, 15]$, y-axis density from 0 to about 1.6. Seven shape–scale pairs sampled at $N = 100{,}000$:

- **Purple — Gamma(0.5, 1.0)** monotone decay from density ~1.6 at the origin; shape $\alpha < 1$ produces the characteristic spike at zero.
- **Red — Gamma(1.0, 2.0)** pure exponential decay (shape = 1 collapses to Exponential(mean = 2)), intercept ~0.5.
- **Orange — Gamma(2.0, 2.0)** mode near $x = 2$ at density ~0.2.
- **Yellow — Gamma(3.0, 2.0)** broader unimodal hump near $x = 4$.
- **Green — Gamma(5.0, 1.0)** mode near $x = 4$ at density ~0.2.
- **Black — Gamma(9.0, 0.5)** sharp symmetric hump centered on $x = 4$ peaking at density ~0.27 (smallest scale → tightest curve at fixed mean).
- **Blue — Gamma(7.5, 1.0)** broad mode near $x = 6.5$ with the longest tail.

The mode sits at $(\alpha - 1)\theta$ and the mean at $\alpha\theta$, matching where each histogram peaks. **The figure demonstrates that the `distributions.gamma` sampler covers the full shape spectrum from monotone-decreasing ($\alpha < 1$), through exponential ($\alpha = 1$), to increasingly symmetric Gaussian-like humps ($\alpha \gg 1$), with the mean tracking $\alpha\theta$.**

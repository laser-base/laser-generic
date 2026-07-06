**Reading the LogNormal plot.**
Overlaid density histograms on $[0, 3]$ (samples clipped at 3), y-axis density from 0 to about 1.65. Three sigma values at fixed $\mu = 0$ (so the median sits at $e^0 = 1$), sampled at $N = 100{,}000$:

- **Blue — Lognormal(0, 1)** broadest curve, mode near $x \approx 0.37 = e^{-1}$ with peak density ~0.8, long right tail extending past $x = 3$.
- **Green — Lognormal(0, 0.5)** mode near $x \approx 0.78 = e^{-0.25}$ with peak density ~0.95, moderate right skew.
- **Red — Lognormal(0, 0.25)** tightest curve, near-symmetric narrow peak at $x \approx 0.94 = e^{-0.0625}$ with peak density ~1.65.

All three share the same median at $x = 1$. As $\sigma$ shrinks the mode rises toward 1 and the skewness vanishes; as $\sigma$ grows the mode slides toward 0 and the right tail thickens. **The figure demonstrates that `distributions.lognormal` correctly produces the positively skewed lognormal PDF and that its variance parameter controls both the spread and the right-tail heaviness.**

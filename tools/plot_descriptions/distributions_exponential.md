#### Reading the Exponential plot

Overlaid density histograms on $[0, 8]$ (samples clipped at 8), y-axis density from 0 to about 1.45. Three rate parameters sampled at $N = 100{,}000$, each constructed by passing `scale = 1/λ`:

- **Red — Exponential(scale = 1/0.5) = Exponential(mean = 2)** the most slowly decaying curve, intercept near 0.5 at $x = 0$ and a long tail extending past $x = 6$.
- **Green — Exponential(scale = 1/1.0) = Exponential(mean = 1)** intercept near 1.0 at $x = 0$, the canonical unit-rate exponential.
- **Light blue — Exponential(scale = 1/1.5)** steepest curve, intercept near 1.45 at $x = 0$, decaying fastest.

Each curve is a clean monotone decay whose intercept equals the rate $\lambda$ — exactly the analytic PDF $\lambda e^{-\lambda x}$. Larger $\lambda$ means a higher density at zero and a shorter mean waiting time $1/\lambda$. **The figure demonstrates that the `distributions.exponential` float sampler reproduces the memoryless exponential PDF and that its rate parameter is correctly mapped through the scale convention.**

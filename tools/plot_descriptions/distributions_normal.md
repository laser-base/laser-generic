#### Reading the Normal plot

Overlaid density histograms on $[-5, 5]$ (samples clipped to this range), y-axis density from 0 to about 0.85. Four mean–variance pairs sampled at $N = 100{,}000$, with the scale parameter passed as $\sqrt{\sigma^2}$:

- **Blue — Normal(0, σ ≈ 0.447) (variance 0.2)** narrow symmetric peak at $x = 0$ reaching density ~0.85.
- **Red — Normal(0, σ = 1.0)** standard normal, peak at $x = 0$ with density ~0.4.
- **Orange — Normal(0, σ ≈ 2.236) (variance 5.0)** widest curve, broad peak at $x = 0$ with density ~0.18, tails extending past $\pm 4$.
- **Green — Normal(μ = −2, σ ≈ 0.707) (variance 0.5)** peak shifted to $x = -2$ with density ~0.56.

Each histogram traces the symmetric Gaussian bell. The peak heights satisfy the analytic $1/(\sigma\sqrt{2\pi})$ relationship (e.g. $1/\sqrt{2\pi} \approx 0.399$ for the unit normal), and the green curve shows the location parameter $\mu$ shifting the entire bell along $x$. **The figure demonstrates that `distributions.normal` reproduces the Gaussian PDF with correct mean and standard deviation across narrow and wide regimes.**

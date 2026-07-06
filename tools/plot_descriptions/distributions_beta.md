**Reading the Beta plot.**
Overlaid density histograms on $[0, 1]$, y-axis density from 0 to about 6.5. Five parameter pairs sampled at $N = 100{,}000$:

- **Red — Beta(0.5, 0.5)** U-shaped with tall spikes at both 0 and 1 reaching density ~6.5, the bathtub of the Jeffreys prior.
- **Blue — Beta(5.0, 1.0)** monotonically rising, peaking at 1.0 with density ~5, the power-law shape that occurs when one pseudo-count outweighs the other.
- **Green — Beta(1.0, 3.0)** monotonically decaying from density ~3 at 0 down to 0 at 1.
- **Purple — Beta(2.0, 2.0)** symmetric unimodal bell centered at 0.5 with peak density ~1.5.
- **Orange — Beta(2.0, 5.0)** right-skewed mode near $x \approx 0.2$ with peak density ~2.3.

Each histogram traces the canonical Beta PDF shape for its $(\alpha, \beta)$ pair. **The figure demonstrates that the Numba-compiled `distributions.beta` sampler reproduces the full Beta family — U-shaped, monotone, symmetric, and skewed regimes — within Monte-Carlo noise.**

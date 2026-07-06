**Reading the Uniform plot.**
Overlaid density histograms on $[-1, 4.2]$, y-axis density from 0 to about 2.5. Six $[low, high]$ pairs sampled at $N = 100{,}000$ — each histogram is a flat rectangle of height $1/(high - low)$:

- **Red — Uniform(0.0, 1.0)** rectangle on $[0, 1]$ at density 1.0.
- **Orange — Uniform(0.25, 1.25)** rectangle on $[0.25, 1.25]$ at density 1.0 (width 1).
- **Green — Uniform(0.0, 2.0)** rectangle on $[0, 2]$ at density 0.5 (width 2).
- **Blue — Uniform(−1.0, 1.0)** rectangle on $[-1, 1]$ at density 0.5 (width 2).
- **Indigo — Uniform(2.71828, 3.14159)** narrow tall rectangle on $[e, \pi]$ at density ~2.4 (width $\pi - e \approx 0.424$).
- **Violet — Uniform(1.30, 4.20)** wide low rectangle on $[1.3, 4.2]$ at density ~0.34 (width 2.9).

Every rectangle sits squarely on its $[low, high]$ support with zero density elsewhere, and the heights match $1/(high - low)$ to within Monte-Carlo noise. **The figure demonstrates that `distributions.uniform` produces flat densities on arbitrary intervals — including offset, negative, and irrational endpoints — with the correct inverse-width normalization.**

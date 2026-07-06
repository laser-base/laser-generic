**Reading the Logistic plot.**
Overlaid density histograms on $[-5, 25]$ (samples clipped to this range), y-axis density from 0 to about 0.25. Five location–scale pairs sampled at $N = 100{,}000$:

- **Purple — Logistic(μ = 2, s = 1)** tallest curve, sharp symmetric peak at $x = 2$ with density ~0.25.
- **Blue — Logistic(μ = 5, s = 2)** peak at $x = 5$ with density ~0.125.
- **Light blue — Logistic(μ = 6, s = 2)** peak at $x = 6$ with density ~0.125, same width as the blue.
- **Green — Logistic(μ = 9, s = 3)** peak at $x = 9$ with density ~0.08, broader heavier tails.
- **Red — Logistic(μ = 9, s = 4)** flattest, peak at $x = 9$ with density ~0.06, widest tails reaching past $x = 20$.

Each curve is the bell shape of the logistic PDF — visually similar to a Gaussian but with noticeably heavier symmetric tails. The peak height scales as $1/(4s)$ and the location sits exactly at $\mu$. **The figure demonstrates that `distributions.logistic` reproduces the analytic logistic PDF with correct location $\mu$ and scale $s$ across both narrow and broad regimes.**

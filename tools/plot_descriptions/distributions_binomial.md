#### Reading the Binomial plot

Discrete spike histogram with integer support from 0 to ~31 on the x-axis, density on the y-axis from 0 to about 1.3. Three parameter pairs sampled at $N = 100{,}000$:

- **Blue — Binomial(n = 20, p = 0.5)** symmetric envelope centered on $np = 10$, peak density ~1.0 at $k = 10$, support roughly 3–17.
- **Green — Binomial(n = 20, p = 0.7)** left-skewed envelope centered on $np = 14$, peak density ~1.3 at $k = 14$.
- **Red — Binomial(n = 40, p = 0.5)** broader envelope centered on $np = 20$, peak density ~0.5, with visibly wider variance ($np(1-p) = 10$) than the $n = 20$ cases ($np(1-p) = 5$).

Each spike sits at an integer, and the spike heights trace the binomial PMF. The mean shifts with $np$ and the spread grows with $\sqrt{np(1-p)}$, exactly as expected. **The figure demonstrates that the `distributions.binomial` integer sampler matches the analytic binomial PMF across small-$n$, skewed-$p$, and larger-$n$ regimes.**

### Reading the Poisson plot

Discrete spike histogram on integer support from 0 to ~20, y-axis density from 0 to about 5.4. Three rate parameters sampled at $N = 100{,}000$:

- **Orange — Poisson(λ = 1)** tall narrow spike at $k = 1$ reaching density ~5.4 (because the bin width is sub-unit), with secondary spikes at $k = 0$ and $k = 2$ near density ~0.5–1.0.
- **Purple — Poisson(λ = 4)** classic right-skewed cluster peaking jointly at $k = 3$ and $k = 4$ with density ~1.4, support mostly between 0 and 10.
- **Light blue — Poisson(λ = 10)** broad nearly-symmetric envelope centered at $k = 10$ with peak density ~0.6, tail extending out past $k = 18$.

Mean and variance both equal $\lambda$ for each curve, and the mode sits at $\lfloor\lambda\rfloor$. At $\lambda = 10$ the distribution already starts to look bell-shaped (the central-limit approach to the normal). **The figure demonstrates that `distributions.poisson` reproduces the analytic Poisson PMF, with the mean tracking $\lambda$ and the spread growing as $\sqrt{\lambda}$.**

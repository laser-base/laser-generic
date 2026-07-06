#### Reading the logistic-fit-with-offset plot

Same log-scale axes as the sanity-check plot — currently-infected versus day, 1 to $10^6$ on the y-axis. Three curves:

- **"Model output"** (blue, solid, thick): the stochastic SI trajectory from this realization.
- **"Logistic growth with known inputs, t0=0"** (orange, dashed): the analytic logistic curve with the known $\beta = 0.05$ and no time offset.
- **"Logistic growth with known inputs, best-fit t0 = -14.0"** (black, dashed): the same analytic curve with a best-fit time offset of $t_0 \approx -14$ days.

The orange curve sits visibly to the left of the model in the exponential-growth region — i.e. the textbook logistic with no offset predicts the outbreak takes off about two weeks sooner than this realization actually did. The black curve, which only adjusts $t_0$, lines up almost perfectly with the blue model output across the entire trajectory. **The takeaway: stochastic delays among the very first infections shift the entire downstream trajectory in time without distorting its shape.** Once the early lag is absorbed into a single offset parameter, the deterministic logistic equation with the known $\beta$ and population reproduces the discrete-time stochastic model — motivating the next section's strategy of fitting $\beta$ jointly with $t_0$ rather than against the wall-clock day count.

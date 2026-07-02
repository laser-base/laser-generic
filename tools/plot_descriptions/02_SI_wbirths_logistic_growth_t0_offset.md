### Reading the logistic-fit-with-offset plot

Same log-scale axes as the sanity-check plot — currently-infected versus day, 1 to $10^6$ on the y-axis. Three curves:

- **"Model output"** (blue, solid, thick): the stochastic SI-with-births trajectory from this realization.
- **"Logistic growth with known inputs, t0=0"** (orange, solid): the analytic generalized-logistic curve evaluated with the known $\beta = 0.04$, $\mathrm{CBR} = 40$, and no time offset.
- **"Logistic growth with known inputs, best-fit t0 = 5.2"** (red, dotted): the same analytic curve with a best-fit time offset of $t_0 \approx 5.2$ days.

All three curves visually overlay across the entire trajectory; the t0=0 and best-fit-t0 lines are nearly indistinguishable, and both track the stochastic blue trajectory closely from the exponential-growth phase through the plateau near $N = 10^6$ around day ~400. The fitted $t_0$ here is small and positive (about a 5-day delay), in contrast to the $\sim -14$-day offset seen in the births-free notebook. **The takeaway: the generalized-logistic solution with the demographic correction factor $x = 1 - \mu/\beta$ reproduces the discrete-time stochastic SI-with-vital-dynamics trajectory using the known $\beta$ and $\mathrm{CBR}$**, confirming both that the analytic form derived in the lead-in markdown is correct and that the residual stochastic timing jitter is small enough to absorb into a single offset parameter — setting up the joint $(\beta, \mathrm{CBR}, t_0)$ fit in the next section.

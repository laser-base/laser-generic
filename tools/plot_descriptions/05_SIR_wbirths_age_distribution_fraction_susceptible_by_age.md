### Reading the fraction-susceptible-by-age plot

Fraction of agents still in the $S$ compartment plotted against age in years (0 to 15 on the x-axis, 0 to 1 on the y-axis), forced through $(0, 1)$. Three series:

- **Blue dots — simulation** fraction susceptible in each age bin (180 bins across 0-15 years), restricted to agents born after day 5000.
- **Orange solid (thick) — expected exponential** $e^{-a \mu (R_0-1)}$ with analytic mean $A = \frac{1}{\mu(R_0-1)} \approx 1.071$ years.
- **Black dashed — best-fit exponential** to the simulation, fitted mean $A = 1.016$ years.

All three series trace essentially the same curve: 1.0 at birth, ~0.5 by age 1, ~0.1 by age 2.5, indistinguishable from zero past age 5. The slight gap between the analytic curve (1.071 years) and the best fit (1.016 years) is the residual stochastic noise from a single realization; the RMSE printed below the figure quantifies it. Unlike the previous plot, this quantity *conditions on survival to age $a$*, so it is not censored by mortality and recovers the uncensored mean $\frac{1}{\mu(R_0-1)}$.

**The figure demonstrates that the fraction of agents still susceptible at age $a$ falls off as an exponential with rate $\mu(R_0-1)$, giving a mortality-uncensored route to estimating $R_0$ from age-stratified susceptibility data.**

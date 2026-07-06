**Reading the SIR compartment plot.**
Four solid/dashed curves plus a horizontal reference line on linear axes — number of people (0 to $10^5$) versus day (0 to ~730), titled "SIR Model with No Births or Deaths":

- **"S"** (blue, solid): susceptible count $S_t$. Starts at $10^5$, holds flat until ~day 40, drops steeply through the epidemic, and settles at $S_\infty \approx 6{,}000$.
- **"I"** (red, solid): infected count $I_t$. Rises from the seed, peaks near $\approx 30{,}000$ around day 80, and decays back to zero by ~day 200.
- **"R"** (green, solid): recovered count $R_t$. Mirrors $S$ as it climbs from 0 to a plateau just below $10^5$.
- **"recoveries"** (orange, dashed): the `newly_recovered` flux $\Delta R_t$. A small bump tracking the recovery flow; peaks well below the $I_t$ peak.
- **"Est. R(∞)"** (gray, dashed, horizontal): analytic $R(\infty)$ from the Kermack-McKendrick / Lambert-W solution, sitting almost exactly on top of the green plateau.

With $R_0 = 3.0$, the textbook attack fraction is just under 80%, and the simulation lands there: $R$ plateaus at $\approx 94{,}000$ and $S$ at $\approx 6{,}000$, both flush with the analytic prediction line. **This plot demonstrates that the stochastic agent-based SIR model reproduces both the dynamic shape and the final-size prediction of the deterministic Kermack-McKendrick equations for a single $R_0 = 3$ outbreak.**

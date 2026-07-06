**Reading the sanity-check plots.**
Two side-by-side panels on log y-axes, both running from time 0 to ~3000 ticks:

- **Left — Susceptible over time.** Y-axis spans $10^5$ to $3\times10^5$. The susceptible count $S_t$ (blue, thick) sits at $N = 3\times10^5$ through day ~150, drops sharply between days ~180 and ~280, and settles into a noisy plateau just below $10^5$ for the remainder. The orange curve overlays the reconstruction $S_{t-1} + \text{recovered}_{t-1} - \text{newly\_infected}_{t-1}$.
- **Right — Currently infected over time.** Y-axis from 1 to ~$3\times10^5$. The infected count $I_t$ (blue, thick) climbs about five orders of magnitude from the single seed up through day ~280, then locks into the endemic plateau near $2\times10^5$. The orange curve overlays the reconstruction $I_{t-1} - \text{recovered}_{t-1} + \text{newly\_infected}_{t-1}$.

In both panels the two curves visually overlay everywhere, and the printed `np.allclose` checks both return `True`. **The takeaway: the SIS compartments are bookkept correctly tick-by-tick — every recovery that subtracts from $I$ adds back to $S$, and every new infection moves the opposite direction — and unlike the SI model the trajectories settle into a finite endemic equilibrium rather than absorbing the whole population into $I$.**

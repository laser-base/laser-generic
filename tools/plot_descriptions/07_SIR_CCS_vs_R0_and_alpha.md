### Reading the CCS-vs-parameter scatter plots

Two stacked scatter panels, one per simulation in the 200-sim parameter sweep, with the smallest persisting patch population (the empirical CCS estimate) on a log y-axis spanning $\sim 10^3$ to $10^6$:

- **Top — CCS_est vs $R_0$.** x-axis from ~1 to ~16, points coloured by $\log_{10}(\alpha)$ (viridis, ~1.75 dark to ~3.5 yellow, where $\alpha = (\gamma+\mu)/\mu$). The cloud is broad and shows only a weak downward trend with $R_0$; the dominant visual structure is the colour gradient — dark (low-$\alpha$) points sit at the bottom near $10^3$–$10^4$, yellow (high-$\alpha$) points pile up against the $10^6$ ceiling.
- **Bottom — CCS_est vs $\log_{10}(\alpha)$.** x-axis from ~1.75 to ~3.6, points coloured by $R_0$ (plasma, ~2 dark to ~15 yellow). Here a clean diagonal trend appears: CCS rises roughly two orders of magnitude as $\log_{10}(\alpha)$ moves from 1.75 to 3.5, with the colour gradient showing that within any vertical slice, higher-$R_0$ points sit lower.

**The figure demonstrates that $\alpha = (\gamma+\mu)/\mu$ — essentially the ratio of infectious period to demographic susceptible-replenishment time — is the dominant driver of CCS, while $R_0$ exerts only a secondary downward influence, exactly as the Nasell formulae predict.** The flat $10^6$ ceiling and $10^3$ floor are simulation artifacts: those are the largest and smallest patches in the grid, so points pegged there are upper/lower bounded rather than informative, motivating the filter applied in the next fitting cell.

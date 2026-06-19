### Reading the 3D best-fit surface plot

A three-dimensional view of the fitted power-law $\text{CCS} = c \cdot \alpha^{a} \cdot (R_0/(R_0-1))^{b}$ against the simulation points:

- **x-axis — $\alpha$** from 0 to ~2000 (linear).
- **y-axis — $R_0$** from ~2 to ~16 (linear).
- **z-axis — $\log_{10}(\text{CCS})$** from ~3 to ~7.
- **Viridis surface** — the fitted function evaluated on a meshgrid, sweeping up steeply in $\alpha$ and curling sharply upward at small $R_0$ where the $R_0/(R_0-1)$ factor blows up.
- **Red dots** — the 200 individual simulation outcomes (filtered to remove the $10^3$ floor and $10^6$ ceiling cases).
- **Title text** prints the fit: $y \approx 0.72 \cdot \alpha^{1.83} \cdot (R_0/(R_0-1))^{1.67}$.

The red cloud sits broadly on or just below the coloured surface, hugging it most tightly along the high-$\alpha$, low-$R_0$ ridge and showing more scatter in the bulk middle. **The figure demonstrates that a two-term power law in $\alpha$ and $R_0/(R_0-1)$ captures the simulated CCS surface across two orders of magnitude in population, with the recovered exponent on $\alpha$ ($\approx 1.83$) close to the $\alpha^{3/2}$ scaling of Nasell's $N_{\text{crit},1}$ formula.**

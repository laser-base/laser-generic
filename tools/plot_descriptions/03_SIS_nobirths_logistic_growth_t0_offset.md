### Reading the SIS logistic-fit-with-offset plot

Log y-axis from $10^{-1}$ to ~$3\times10^5$, time on the x-axis from 0 to 3000 ticks. Three curves:

- **Blue solid (thick) — model output**, the stochastic SIS trajectory from this single realization.
- **Orange solid — analytic SIS logistic** $\frac{Nx}{1+(Nx/I_0-1)e^{-\beta x t}}$ with $x = 1-\gamma/\beta$, evaluated at the known $\beta = 0.1$, $\gamma = 1/32$, and $t_0 = 0$.
- **Red dotted — same analytic curve** with the best-fit time offset $t_0 \approx 36$ days.

The orange curve takes off noticeably earlier than the blue model output during the exponential phase, then all three curves merge onto the same endemic plateau near $2\times10^5$ — the analytic prediction $Nx = N(1-\gamma/\beta) \approx 0.6875 \times 3\times10^5$. The red curve, shifted right by 36 days to absorb the stochastic delay in the very first transmissions, lies on top of the blue trajectory for the entire growth phase. **The takeaway: the same single-parameter time-offset trick that worked for the SI model also collapses the SIS realization onto the deterministic logistic, confirming that the discrete-time stochastic SIS reproduces the analytic shape — and the shared endemic plateau confirms the model is hitting the correct $1-\gamma/\beta$ equilibrium fraction.**

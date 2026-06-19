### Reading the case-trace-and-fit plot

A single linear-scale plot titled "Case Trace and Logistic Fit for Row 0" showing the last sweep row (index 9: $\beta = 0.055$, $\mathrm{CBR} = 30$) over 1825 days (5 years):

- **Blue solid — Case Trace.** The model's $I_t$ from `model.nodes.I[:, 0]` for this run.
- **Orange dashed — Logistic Fit.** The generalized-logistic curve evaluated with the four fitted parameters $\hat{\beta}$, $\hat{N}$, $\hat{\mathrm{CBR}}$, $\hat{t_0}$ recovered for this row.

Both curves trace a smooth S-shape: flat near zero through day ~100, a steep rise between day ~150 and day ~250 to the carrying capacity, then a long plateau at $\approx 100{,}000$ cases that holds out to day 1825. **The two curves visually overlay across the entire 5-year window** — no residual gap is visible even on the linear scale that exaggerates plateau-level discrepancies. **The takeaway: despite the joint $(\beta, \mathrm{CBR})$ fit showing >5% relative error for some sweep rows in the previous plots, the fitted logistic curve still reproduces the model's case trajectory to visual precision**, which is why the lead-in markdown notes the consistently-negative fitted $t_0$ (here $\hat{t_0} \approx -28$) as a quirk of seeding with multiple infections rather than a sign of a bad fit.

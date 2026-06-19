### Reading the sanity-check plot

Three overlapping curves on a shared log y-axis (1 to $10^6$) versus time in days (0 to ~730):

- **Blue solid (thick) — currently infected** $I_t$ from `model.nodes.I`.
- **Orange dashed — population minus susceptibles** $N_t - S_t$, the same quantity reconstructed from the susceptible compartment.
- **Black dotted — initial seed plus cumulative incidence** $I_0 + \sum_{t'} \Delta I_{t'}$, the same quantity reconstructed from the per-tick new-infection records.

The shape is the classic SI take-off-and-saturate: roughly four orders of magnitude of exponential growth from day 0 through day ~250, then a sharp transition to the plateau at $N = 10^6$ around day 300. **All three curves visually overlay across the entire range.** That overlap is the sanity check passing: the model preserves $N = S + I$ at every tick, and cumulative incidence accumulated from `newly_infected` equals the running $I_t$ minus the initial seed — exactly the two identities the lead-in markdown asserts.

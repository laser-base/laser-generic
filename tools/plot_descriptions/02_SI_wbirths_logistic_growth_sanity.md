### Reading the sanity-check plot

Two nearly-overlapping curves on a shared log y-axis (1 to $10^6$) versus time in days (0 to ~730):

- **Blue solid (thick) — currently infected** $I_t$ from `model.nodes.I`.
- **Orange dashed — population minus susceptibles** $N_t - S_t$, the same quantity reconstructed from the susceptible compartment.

(The legend still lists a third "cumulative incidence" entry from the births-free notebook, but only two curves are plotted here — with `ConstantPopVitalDynamics` recycling agents, cumulative incidence from `newly_infected` no longer equals $I_t - I_0$, since some early infectives die and are replaced by newborn susceptibles.) The trajectory shows the same SI take-off-and-saturate shape as the no-births notebook: roughly six orders of magnitude of exponential growth from day 0 through day ~400, then a plateau at $N = 10^6$. **The blue and orange curves visually overlay across the entire range**, and the cell prints `S = N-I: True`. **That overlap is the surviving sanity check: even with births and deaths constantly shuffling individuals between compartments, the model preserves $N = S + I$ at every tick**, which is exactly the constraint the lead-in markdown promises `ConstantPopVitalDynamics` will enforce.

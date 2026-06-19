### Reading the time-varying-CBR population-growth plot

A linear-axis plot of total population versus time in days (0 to 3650), with the same faint grid styling as the previous plot:

- **Blue solid — simulated total population** $S+E+I+R$ from the temporal-CBR model.
- **Orange dashed — expected growth** computed by stepping the analytic per-day compounding factor $(1 + \text{CBR}_t/1000)^{1/365}$ forward through the same 3650-day CBR ramp.

The trajectory starts at exactly 100,000 and grows to roughly 122,000 by day 3650. Unlike the constant-CBR plot, the slope visibly decreases over time: early years (CBR near 30 per 1000) add population fast, while later years (CBR sliding toward 10 per 1000) add it more slowly, giving the curve a concave-downward bend rather than the constant-CBR's slight exponential bow. The two curves again overlay tightly across the full range, consistent with the printed difference of well under a percent. **The takeaway: `BirthsByCBR` driven by a `ValuesMap.from_timeseries` correctly tracks a CBR schedule that varies day-by-day, matching the time-integrated analytic growth and confirming that the births component reads its per-tick rate from the time axis rather than freezing the initial value.**

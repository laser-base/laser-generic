### Reading the SEI vs SEIS mean-trajectory plot

A single linear-scale panel, 0–365 days on the x-axis, 0 to $10^6$ agents on the y-axis. Six curves comparing seed-averaged means — solid for SEI, dashed for SEIS — sharing the same color per compartment.

- **Blue (S):** SEI solid drops from $10^6$ to ~0 by day ~100. SEIS dashed tracks it closely but pulls up to an endemic plateau near ~215,000 from day ~120 onward.
- **Green (E):** SEI solid produces a transient pulse peaking ~215,000 around day 80, returning to zero by day ~120. SEIS dashed peaks slightly later (~225,000) and instead of decaying settles to ~205,000.
- **Orange (I):** SEI solid sigmoidally climbs to the absorbing $10^6$ saturation by day ~100. SEIS dashed peaks near ~600,000, then dips slightly and settles at the endemic ~580,000.

The SEIS curves also start their rise a few days later than SEI — recovery back to $S$ continuously erodes the force of infection, slowing the early exponential climb. **The figure demonstrates the qualitative difference between an absorbing-state SEI (everyone ends infectious) and a recycling SEIS (system relaxes to a non-trivial endemic equilibrium where $S$, $E$, and $I$ coexist indefinitely), driven entirely by whether the $I$ compartment leaks back to $S$ via the infectious-duration timer.**

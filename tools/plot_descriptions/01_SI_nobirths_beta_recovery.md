**Reading the β-recovery plots.**
The cell emits two scatter plots side by side from a 10-seed sweep over $\beta \in \{0.02, 0.03, \ldots, 0.11\}$:

- **Left — Fitted β vs True β.** Ten points lying very close to the identity line on a 0.02–0.11 axis. Recovery is near-perfect at low $\beta$ and drifts slightly below the diagonal as $\beta$ grows.
- **Right — Relative error in fitted β vs True β.** Same x-axis; y-axis runs from ~0.001 (0.1%) at the smallest $\beta$ to ~0.024 (2.4%) at the largest. The relationship is roughly linear with some scatter from seed-to-seed stochasticity.

**The takeaway:** the discrete-time finite-difference approximation systematically under-shoots the continuous-time $\beta$, with the bias growing roughly linearly as the per-step transmission probability gets larger — exactly the qualitative bias predicted in the notebook's lead-in markdown. Every point sits well inside the 5% tolerance, which is what makes the next cell's pass/fail assertion return `True`.

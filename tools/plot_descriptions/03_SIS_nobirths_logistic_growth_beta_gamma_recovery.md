### Reading the β- and γ-recovery plots

The cell emits four scatter plots from a 10-seed sweep with $\beta \sim U(0.02, 0.1)$ and $\gamma \sim U(1/300, 1/100)$:

- **Panel 1 — Fitted β vs True β.** Square frame 0.02–0.10 on both axes. Ten points hug the identity line; recovery is essentially unbiased across the range.
- **Panel 2 — Relative β error $1-\beta/\hat\beta$ vs True β.** X-axis 0.02–0.10, y-axis ±0.25. All ten points sit just below zero (roughly −0.01 to −0.03), confirming the same small, systematic discrete-time under-shoot the SI model exhibited.
- **Panel 3 — Fitted γ vs True γ.** Square frame ~0.003–0.011 on both axes. Points cluster around the identity line but with visibly more scatter than the β panel — a few points are noticeably below the diagonal.
- **Panel 4 — Relative γ error $1-\gamma/\hat\gamma$ vs True γ.** X-axis 0.003–0.011, y-axis ±0.25. Points sit between roughly −0.02 and −0.08, all negative, indicating the fitted $\gamma$ is systematically smaller than the true $\gamma$ but well inside the 20% tolerance.

**The takeaway: the joint fit recovers both $\beta$ and $\gamma$ from a single stochastic SIS realization with bounded bias — β within ~3%, γ within ~8% — which is what makes the next cell's `< 0.10` / `< 0.20` pass-fail asserts return `True` and validates the SIS implementation against its analytic logistic solution.**

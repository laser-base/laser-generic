#### Reading the β- and CBR-recovery plots

The cell emits four scatter plots from a 10-seed sweep over $\beta \in \{0.010, 0.015, \ldots, 0.055\}$ paired with randomly drawn crude birth rates $\mathrm{CBR} \in [15, 50)$:

- **Plot 1 — Fitted β vs True β.** Axes 0.00–0.06 on both sides. Ten points cluster tightly along the identity line. Notably, two points near $\beta = 0.040$ sit visibly *above* the diagonal (fitted ≈ 0.047) — an over-shoot rather than the systematic under-shoot seen in the births-free notebook.
- **Plot 2 — Relative error in fitted β vs True β.** Same x-axis; y-axis runs from $-0.1$ to $+0.1$. Most points sit between 0 and $+0.025$ (i.e. true is slightly higher than fitted, a small under-shoot); the cluster is roughly flat across the $\beta$ range rather than fanning out.
- **Plot 3 — Fitted CBR vs True CBR.** Axes 15–50 on both sides. Points lie along the identity line with visible scatter; the cluster near true CBR = 43 has one point that jumps up to fitted ≈ 49.
- **Plot 4 — Relative error in fitted CBR vs True CBR.** Same x-axis; y-axis $-0.2$ to $+0.2$. Most points are within $\pm 0.025$, but one outlier near true CBR = 43 sits at roughly $-0.13$ (about a 13% over-shoot in the fit).

**The takeaway: joint recovery of $\beta$ and $\mathrm{CBR}$ from the generalized-logistic fit is mostly tight but no longer uniformly inside the 5% tolerance — the added birth/death degree of freedom and short demographic transients introduce outliers**, which is why the subsequent assertion cell prints `False` for both the 10% β and 10% CBR pass/fail checks.

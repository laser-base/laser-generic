**Reading the MortalityByCDR observed-CDR sweep plot.**
A 2x2 grid of subplots, each sharing the same axes: x-axis is run index 1-11, y-axis is observed CDR (deaths per 1,000 per year) from 0 to ~42. Each panel corresponds to one configured CDR target: top-left CDR=2, top-right CDR=10, bottom-left CDR=20, bottom-right CDR=40. In every panel blue x-markers show the 11 individual run observations and a green dashed line marks the across-run mean.

- **CDR=2 panel:** all 11 markers cluster tightly on the green line at mean 2.01.
- **CDR=10 panel:** markers sit on mean 9.97 with barely visible run-to-run scatter.
- **CDR=20 panel:** mean 19.87, very slight visible scatter (~0.2 spread).
- **CDR=40 panel:** mean 39.82 with the largest visible spread of the four (still under ~1 unit), consistent with Poisson-like variance growing with the rate.

In every panel the observed mean falls within ~1% of the configured CDR. **This figure demonstrates that `MortalityByCDR` is statistically unbiased: across two orders of magnitude in target rate (2 to 40 per 1,000 per year), the empirically observed annual mortality reproduces the input CDR with low run-to-run variance, validating the component as a calibrated non-disease mortality driver.**

### Reading the per-patch infection heatmap

A single heatmap from one of the 200 simulations showing $\log_{10}$ of the prevalence (infected fraction) in each patch over time:

- **x-axis — time in days** from 0 to ~14,600 (40 years).
- **y-axis — patch population** on a log scale, with 61 patches arrayed from $10^3$ (bottom) to $10^6$ (top), tick labels at 1000, 3163, 10000, 31623, 100000, 316228, 1000000.
- **Colorbar labeled "Cases"** (viridis): despite the label, encodes $\log_{10}(I/N)$, roughly $-5$ (dark) to $-2$ (yellow).

The top rows (large patches, $\gtrsim 10^5$) show continuous coloured streaks across the whole 40-year span — the disease persists. Below ~$10^4$ the picture is dominated by white (zero-infection) gaps punctuated by short coloured stripes that mark each importation pulse (recall the `Importation_EachNode` component reseeds every 180 days for the first 20 years). After day ~7,300 the importations stop, and the small-patch rows go essentially white through the end of the run, while the large-patch rows keep ticking. **The figure demonstrates the CCS phenomenon visually: above some threshold population the SIR dynamics sustain themselves, below it the disease fades out between importations and cannot survive once reseeding ends.**

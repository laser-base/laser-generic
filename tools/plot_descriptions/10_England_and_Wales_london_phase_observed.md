#### Reading the observed London wavelet-phase plot

A single scatter of blue dots: **wavelet phase difference (degrees) versus distance from London**. The x-axis runs 5 to 30 (distance units as in the dataset) and the y-axis runs from $-90$ to $0$, with the title "Phase difference of London wavelet transform".

Each point is one England-and-Wales place within 30 distance units of London. Phase is extracted from the complex Morlet wavelet transform of weekly cases in the 2–3 year (biennial) band, then the angle of the mean cross-spectrum against London is taken.

The cloud shows a clear **downward trend**: places close to London (distance $\sim 8$) sit near $0^\circ$ to $-10^\circ$, and points drift toward $-30^\circ$ to $-50^\circ$ by distance $\sim 28$–30. A linear lag of roughly $-1.5^\circ$ to $-2^\circ$ per distance unit is visually evident through the scatter.

**The plot demonstrates the Grenfell–Bjørnstad–Kappey "travelling wave" signature: biennial measles epidemics arrive later in places farther from London, producing a systematic phase lag that the calibrated spatial model must reproduce.** This observed pattern becomes the target for the wavelet-phase calibration metric used in the next set of cells.

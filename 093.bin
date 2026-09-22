**Provenance correction (22 September):** Natalie’s notebooks perform buffer subtraction on the HDF5 input profiles. Both AA_082024.ipynb and AA_052026.ipynb load AA_052026.hdf5 and report buffer frames24–70; the primary WT analysis reports220–236. Earlier claims that the needed pre-subtraction profiles were absent, or that our HDF5-based subtraction was necessarily only a residual correction, were incorrect. The HDF5-based calculations already read `profiles` and subtract their chosen buffer once, so this correction does not itself change their numerical output. Original exported DAT curves are a later processing stage and must not be subtracted again. The new controlled test is in [09_buffer_window_sensitivity.ipynb](../09_buffer_window_sensitivity.ipynb).

## What do pre, linear and post mean?
These are **three alternative residual-background subtractions**, applied to the input curves before fitting. They are not EFA components, different polymer models, or three measured protein states. The HDF5 provenance does not establish the upstream subtraction history, so these anchors must not be described as verified original buffer-only exposures.

| Sample | Pre anchor frames | Post anchor frames |
|---|---|---|
| WT |55–94|210–234|
| AA |35–74|225–254|

For each q, average the pre-anchor frames to get $B_{pre}(q)$, and the post-anchor frames to get $B_{post}(q)$. Each is an entire q-dependent curve, not a single constant offset. For a protein frame t, subtract:

$$\widehat B(q,t)=\begin{cases}B_{pre}(q)&\text{pre}\\B_{post}(q)&\text{post}\\(1-f_t)B_{pre}(q)+f_tB_{post}(q)&\text{linear}\end{cases}$$

where $f_t=(t-t_{pre})/(t_{post}-t_{pre})$, clipped to [0,1], and the anchor times are the midpoints of their frame ranges. The three corrected estimates are $I(q,t)=Y(q,t)-\widehat B(q,t)$, then averaged over the selected protein window.

If the true signal is $Y(q,t)=c_tP(q)+B(q,t)$, a baseline error leaves $\widehat I=c_tP+\delta B$. Its relative effect is $\delta B/(c_tP)$: it grows as the protein signal falls on the peak edges. Because Rg comes from the low-q slope, $R_g^2=-3\,d\ln I/d(q^2)$ at q approaching0, a q-dependent residual background changes Rg, not just intensity. Extended Guinier accounts for polymer curvature but cannot fix an incorrect subtraction.

**Read the bottom-right panels with their different y-axis scales in mind.** WT spans roughly26–33Å; AA spans only31.6–32.3Å. The WT pre-corrected leading/tail can resemble AA, but the WT center does not:

| Window | WT pre / linear / post Rg (Å) | AA pre / linear / post Rg (Å) |
|---|---|---|
| Leading |31.86 /29.80 /27.94|32.11 /31.80 /31.57|
| Center |28.17 /27.20 /26.56|32.31 /32.20 /32.14|
| Trailing |28.75 /27.20 /26.42|32.10 /31.96 /31.90|
| Late tail |32.77 /27.74 /25.79|32.23 /31.87 /31.75|

These values use extended fits over the same0.012–0.055Å⁻¹ range. WT windows are150–159/160–169/170–179/180–189; AA windows are160–169/170–179/180–189/190–199. They are direct averages, not EFA-purified profiles. The newer five-frame leading-edge test (WT161–165;AA165–169) is a separate selection.

The top-right and bottom-left panels of the original four-panel window plots use the **linear** baseline only. Their shape normalization is multiplicative; it does not change Rg. WT's leading window also overlaps the earlier shoulder and shows structured shape residuals. Its late tail is weaker and particularly baseline-sensitive. Neither observation assigns an oligomer.

These tested corrections preserve a smaller WT center than AA center, while leaving uncertainty in the absolute size and edge estimates. This does not validate linear subtraction as the true baseline: anchor quality, residual time dependence, and original sample/buffer exports must decide the background model—not which model produces the preferred Rg. The spread is systematic sensitivity, not a statistical confidence interval.

The existing anchor stationarity screen passes WT55–94 (maximum band slope magnitude2.97 standard errors; half-window profile discrepancy0.59 in its mean-square metric) and WT210–234 (1.45 and0.71). These tests assess time stability, not absence of protein. The pre anchor is low and flat on the full trace; the post anchor is elevated. The reason for that difference remains unresolved.

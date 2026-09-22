**Provenance correction (22 September):** Natalie’s notebooks perform buffer subtraction on the HDF5 input profiles. Both AA_082024.ipynb and AA_052026.ipynb load AA_052026.hdf5 and report buffer frames24–70; the primary WT analysis reports220–236. Earlier claims that the needed pre-subtraction profiles were absent, or that our HDF5-based subtraction was necessarily only a residual correction, were incorrect. The HDF5-based calculations already read `profiles` and subtract their chosen buffer once, so this correction does not itself change their numerical output. Original exported DAT curves are a later processing stage and must not be subtracted again. The new controlled test is in [09_buffer_window_sensitivity.ipynb](../09_buffer_window_sensitivity.ipynb).

# WT and AA windows, and visible EFA component curves

These figures are embedded so no external image loading is required. Frame numbering is zero-based. WT windows: 150–159, 160–169, 170–179, 180–189. AA windows: 160–169, 170–179, 180–189, 190–199. Same ten-frame lengths, baseline anchors and propagation as documented in the comparison scripts. The shared extended fit interval is 0.012–0.055 Å^-1.

| WT window | Conventional Rg (Å) | Extended Rg (Å) |
|---|---:|---:|
| Leading | 33.82 | 29.80 |
| Center | 28.57 | 27.20 |
| Trailing | 28.86 | 27.20 |
| Late tail | 32.41 | 27.74 |

The leading-minus-center extended Rg has a conditional paired 95% interval of +1.25 to +3.75 Å (300 perturbations). The difference remains positive in pre-only and post-only checks, with the latter marginal. Leading versus center shape discrepancy chi-square/dof=3.39 versus 0.91 for trailing and 0.67 for tail. The tail is weak and background sensitive. No component identity is established.

The WT shoulder EFA curve is positive and structured, unlike AA's secondary curve. The exploratory extended fit is 40.9 Å with qmax Rg<2, but reduced chi-square=2.44 and qmin sensitivity gives 37.6–43.3 Å. Its single-chain 120-residue assumption is unverified; no oligomer stoichiometry follows. Error bars for secondary curves condition on fixed elution coefficients, whereas main-component errors are bootstrap estimates. Curve intensity scales are arbitrary and not directly comparable as populations. Negative AA intensities remain visible on a symmetric-log axis.

Recreate figures with `component_visibility/plot_components.py`; reproduce direct-window results with the two `compare_windows.py` scripts. They currently reference the original Desktop source HDF5 paths; edit SOURCE for another checkout. The Rg model implementation is included beside each script.

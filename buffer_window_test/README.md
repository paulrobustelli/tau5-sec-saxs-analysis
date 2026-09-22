# Controlled buffer-window sensitivity
**Provenance correction (22 September):** Natalie’s notebooks perform buffer subtraction on the HDF5 input profiles. Both AA_082024.ipynb and AA_052026.ipynb load AA_052026.hdf5 and report buffer frames24–70; the primary WT analysis reports220–236. Earlier claims that the needed pre-subtraction profiles were absent, or that our HDF5-based subtraction was necessarily only a residual correction, were incorrect. The HDF5-based calculations already read `profiles` and subtract their chosen buffer once, so this correction does not itself change their numerical output. Original exported DAT curves are a later processing stage and must not be subtracted again. The new controlled test is in [09_buffer_window_sensitivity.ipynb](09_buffer_window_sensitivity.ipynb).

## Question and fixed choices
Subtract the mean of **25–35**, **60–70**, or **210–220** (11 frames each, inclusive and zero-based) directly from the input profiles. Include Natalie’s buffer windows as references. Do not add a subsequent baseline correction. The AA reference is therefore her buffer-subtraction step, not her entire later endpoint-baseline workflow.

Keep protein windows fixed: WT selected159–175, leading161–165, apex165–169; AA selected169–188, leading165–169, apex175–179. Keep EFA ROIs, two-component supports, and q range unchanged. Conventional fits use0.012–0.030Å⁻¹; extended fits use0.012–0.055Å⁻¹. Report qRg validity flags; additional adaptive extended fits are retained in JSON. EFA-other common-range fits may be unsupported because the range is too wide or the curve is nonphysical; they are not species assignments.

For each sample,100 identical input-frame perturbations are propagated through all four subtraction choices. Shared-buffer covariance and overlapping frames are retained algebraically. EFA is rerun on each draw. Failed rotations remain failures; intervals for EFA are conditional on successful rotations and fixed supports. BIFT is rerun for selected averages and successful EFA-main estimates with identical q/Dmax/regularization search settings and80 local parameter samples. This does not supply between-experiment uncertainty.

## Findings
The selected-average extended Rg is **WT28.39/28.17/26.60Å** and **AA32.09/32.25/31.98Å** for25–35/60–70/210–220. Thus AA remains larger, but WT≈26Å is not invariant to buffer selection. Conventional Guinier and P(r) show additional method sensitivity and must not be pooled as interchangeable Rg estimates.

AA210–220 fails stationarity (maximum band slope≈7 standard errors); it is a deliberately tested, questionable buffer choice, not a recommended clean buffer. The AA60–70 nominal two-component EFA rotation fails. AA25–35 converges on only63/100 perturbations, versus100/100 for WT. These failures argue against relying on AA forced two-component extraction as the primary size estimate.

See results.json, paired_differences.json and anchor_checks.json for all numerical results. The executed notebook includes all I(q), extended-fit, Kratky and BIFT comparisons. `analysis/test_buffer_windows.py` reruns the experiment.

**Provenance correction (22 September):** Natalie’s notebooks perform buffer subtraction on the HDF5 input profiles. Both AA_082024.ipynb and AA_052026.ipynb load AA_052026.hdf5 and report buffer frames24–70; the primary WT analysis reports220–236. Earlier claims that the needed pre-subtraction profiles were absent, or that our HDF5-based subtraction was necessarily only a residual correction, were incorrect. The HDF5-based calculations already read `profiles` and subtract their chosen buffer once, so this correction does not itself change their numerical output. Original exported DAT curves are a later processing stage and must not be subtracted again. The new controlled test is in [09_buffer_window_sensitivity.ipynb](../09_buffer_window_sensitivity.ipynb).

# Exported-curve comparison

Run `07_results_comparison.ipynb` from the repository root with the SAXS kernel (or a Python kernel with requirements installed). It reads the registry and exported curves; no EFA rerun is required.

`curve_registry.json` catalogs81 curves. `fit_results.json` records source SHA256, conventional and extended fits, conditional fit errors, fixed-range Monte Carlo intervals and qmin sensitivity. `paired_results.json` records correlated differences. `edge_selection.json` records the leading-edge selection rule and fitted contributions.

## What is fitted

- Conventional weighted log-intensity Guinier, widest tested q range satisfying qmaxRg≤1.1.
- Zheng–Best extended model with120 residues and119 peptide bonds, weighted intensity objective, widest tested range satisfying qmaxRg≤2.
- Extended model over the conventional range to separate fit-form from range effects.
- Saved0.012–0.030 conventional comparisons retained and flagged when outside the IDP range.
- qmin0.010,0.015 and0.020 sensitivity in addition to primary0.012Å⁻¹.

Range criteria are necessary screening rules, not evidence of model validity; inspect residuals. A converged component fit does not establish a molecular species. WT shoulder single-chain fits are exploratory; negative/noise-like components can yield unsupported fits.

## Errors

Each curve uses its exported sigma. For standard log-space fits, the slope covariance is propagated through Rg=sqrt(-3 slope). Extended local errors use the Jacobian and the derivative of Rg with respect to nu. These are conditional on fixed baseline, support and q range.

The notebook uses200 independent-q normal perturbations per curve, except existing100 paired EFA perturbations for main/average. Report-hypothesis components have fixed-C errors from their saved frame errors; unrecorded cross-frame covariance is unavailable. Signed intensities are retained; conventional log fits reject nonpositive points within the fit window instead of silently dropping them.

For direct-window differences,300 joint perturbations use Cov[I_j(q),I_k(q)]=sum_t sigma(q,t)^2 w_tj w_tk, including common buffer anchors and overlap. EFA average/main differences use existing matched perturbations. Fit windows are fixed for these intervals. Baseline, component support, q selection and selection of the leading window contribute additional uncertainty not represented by the intervals.

Existing P(r) curves are compared and their Rg integrals recomputed. Full P(r) covariance is not available, so no new P(r) Rg confidence interval is manufactured.

## Leading edges

For each sample, choose five consecutive pre-apex frames with main coefficient≥10% of its maximum, minimizing the worst other-component absolute scattering fraction over q0.012–0.030 and0.012–0.15. A5% threshold is descriptive. AA165–169 passes (~1.23%); WT161–165 does not (~10.04%). Main/apex averages and pre/linear/post baseline variants are exported for direct comparison. Support-constrained zeros do not prove absence of overlap.

`analysis/export_comparison_inputs.py` rebuilds derived inputs/covariances from the saved data and decompositions. `analysis/compare_all_iq.py` contains the reusable analysis. Notebook07 recomputes fits without overwriting input profiles. Original input subtraction provenance remains incomplete.

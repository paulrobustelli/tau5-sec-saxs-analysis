**Provenance correction (22 September):** Natalie’s notebooks perform buffer subtraction on the HDF5 input profiles. Both AA_082024.ipynb and AA_052026.ipynb load AA_052026.hdf5 and report buffer frames24–70; the primary WT analysis reports220–236. Earlier claims that the needed pre-subtraction profiles were absent, or that our HDF5-based subtraction was necessarily only a residual correction, were incorrect. The HDF5-based calculations already read `profiles` and subtract their chosen buffer once, so this correction does not itself change their numerical output. Original exported DAT curves are a later processing stage and must not be subtracted again. The new controlled test is in [09_buffer_window_sensitivity.ipynb](09_buffer_window_sensitivity.ipynb).

# Actual EFA component curves — 22 September 2026

The requested **component-level I(q) curves now exist**. These are a tested two-component hypothesis, not verified monomer-only scattering. Open `04_EFA_component_Iq.ipynb`, which has been executed, then `component_results/component_Iq_comparison.png`.

## Direct curve files

- `WT_092025_EFA_candidate_main.dat`
- `AA_052026_EFA_candidate_main.dat`

Both are under `component_results/`, with columns **q (Å^-1), candidate main-component I(q), conditional bootstrap sigma**. They retain the original full q grid; the EFA fit uses only q=0.012–0.25 Å^-1. High-q values are a projection using the fitted elution coefficients, not independently resolved components. The arbitrary intensity scale is the fitted component's contribution to the selected frame average. It is not an absolute concentration or a mole fraction.

The corresponding `*_frame_average_same_baseline.dat`, `*_EFA_other_component.dat`, `*_EFA_elution.csv`, and `*_rank1_dominant_not_purified.dat` make the comparison auditable. Rank-one denoising is explicitly not represented as purification.

## Same-baseline, same-q-range comparison

Guinier fits use q=0.012–0.030 Å^-1 throughout. Candidate fits use the same bootstrap uncertainty column exported in their DAT files; averages use propagated frame/background errors. Fixed-elution fits are retained separately in the audit JSON, and may differ slightly because their error weights differ.

| Sample | Frame-average Rg (Å) | EFA candidate Rg (Å) | Candidate bootstrap 95% (Å) | Paired change 95% (Å) |
|---|---:|---:|---|---|
| WT_092025 | 28.60 | 27.42 | 25.15–29.29 | -2.32 to -0.26 |
| AA_052026 | 32.15 | 31.80 | 30.33–33.28 | -1.50 to +0.87 |

These intervals condition on the baseline, component count and support windows; they do not encompass all model ambiguity.

**WT:** the earlier shoulder supplies a distinguishable second scattering pattern, and the candidate main-peak curve has a lower apparent size than the matched average. The earlier pattern peaks near frame 140 and the main candidate near 167. This supports a conditional shoulder/main separation, not identification of either species. Absolute WT size remains materially baseline dependent.

**AA:** the forced two-component solution has almost coincident elution peaks and an extra curve that is largely noise/negative. Its second singular value remains below the supplied-error rank-one null threshold, and the conditional candidate-minus-average Guinier interval includes zero. The new AA curve is an explicit candidate for inspection, **not demonstrated removal of oligomer contamination**.

## Baseline and support choices

The working residual baseline interpolates q-dependent anchor means: WT 55–94 / 210–234 and AA 35–74 / 225–254. All four pass the declared slope/split-profile stationarity screen, but this does not establish absence of sample scattering. The WT pre-anchor high-q drift is close to the screening threshold. AA 230–259 was rejected for drift. Pre-only and post-only corrections, anchor shifts, multiple elution supports, and qmax 0.15 versus 0.25 Å^-1 are provided as sensitivity checks. This is not a uniquely established optimal solvent subtraction. Original unsubtracted sample/buffer frames are still required to redo the facility subtraction.

The fiducial EFA intervals/supports are WT 120–205 with [120,175]/[140,205], and AA 145–215 with [145,185]/[160,215]. The curves are not selected to match a desired Rg. `rotation_grid.json` retains all 144 hypotheses, including failures. The sensitivity chart uses explicit numerical and elution-identity screens, not a claim of biological validity. AA's extra component does not become valid merely because the main curve passes these screens.

The baseline transform is linear and propagated with its shared covariance. Component curves are recovered as S=D(C^T)^+. We ran 100 full-data perturbations with fresh rotations per sample, 200 rank-one null simulations per sample, an exact synthetic two-component recovery check, an exact linear-baseline subtraction check and a Monte Carlo variance check. Unknown correlations from the original facility subtraction remain a limitation.

## P(r) control: a cleaner-looking tail is not evidence by itself

The BIFT comparison includes the unmodified average fitted with candidate error weights. For AA, this alone substantially shortens the weak tail and lowers P(r) Rg. Thus much of the apparent P(r) cleanup is attributable to uncertainty weighting/regularization, not demonstrated removal of an oligomer. The crossed error-weight controls are in `bift_weight_controls.json` and the dashed curve in the main figure. BIFT uses diagonal errors and does not incorporate full q covariance from EFA; its conditional uncertainties and Dmax must not be overinterpreted.

## What remains unresolved

AA direct-average fixed-range Rg stays around 31–32 Å across the leading/central/trailing main peak with larger tail errors; see `elution_window_check.json`. This is consistent with a persistent main-peak size but cannot exclude co-eluting contaminants. Some equally numerical AA decompositions yield substantially different candidate sizes when baseline/support choices change. The available data therefore do not prove either intrinsic AA expansion or an oligomer-contamination explanation.

The practical output is now the actual candidate I(q), the putative removed pattern, reconstruction residuals and uncertainty/sensitivity evidence. A **validated purified AA monomer I(q)** has not been established.

# Poster-style SEC-SAXS figures

Recommended averaging windows, zero indexed and inclusive:

| Sample | Frames | Number averaged | Guinier Rg | BIFT Rg | BIFT Dmax |
|---|---|---:|---:|---:|---:|
| WT | 159–175 | 17 | 25.8 ± 0.6 Å | 27.2 Å | ~94 Å |
| W397A/W433A | 169–188 | 20 | 32.3 ± 0.5 Å | 34.6 Å | ~159 Å |

Guinier uncertainties are formal one-standard-error estimates with propagated residual-background uncertainty. They exclude background choice and other systematic effects. BIFT conditional errors are approximately 0.2 / 0.4 Å for Rg and 4 / 11 Å for Dmax (WT / AA), and also exclude systematic effects. Treat Dmax as model dependent, particularly the weak AA tail.

## Files

- `poster_saxs_panel.png` and `.svg`: combined I(q), P(r), WT Guinier and AA Guinier panel. PNG is 300 dpi; SVG preserves editable text and vector geometry.
- `Iq_normalized.*`: I(q)/I(0) comparison, with Guinier I(0) normalization and BIFT fit lines.
- `Iq_measured_scale.*`: I(q) in the original arbitrary intensity scale. Differences in scale cannot be interpreted as structural differences without concentration calibration.
- `Pr_comparison.*`: P(r)/I(0), using each BIFT I(0). The x axis is **pair distance r**, not Rg.
- `WT_092025_Guinier.*`, `AA_052026_Guinier.*`: individual Guinier plots with standardized residuals.
- `SEC_selected_windows.*`: SEC curves, selected intervals, and five-frame moving apparent Rg estimates. Only frames above 12% of the local trace maximum with valid fits are shown for Rg; all SEC intensities are shown. These overlapping Rg averages are correlated.
- `fit_and_window_diagnostics.*`: BIFT residuals and nearby-window Rg sensitivity.
- `Kratky_companion.*`: optional dimensionless Kratky comparison matching the additional analysis type on the poster.
- `*_selected_average.dat`: q, arithmetic mean I(q), and propagated uncertainty, including shared residual-background uncertainty. The data retain the full original q range, including values not used in fits or figures.
- `*_pr_primary.dat`: primary P(r) values and RAW's conditional uncertainty.
- `window_candidates.csv`, `selected_windows.json`, `bift_sensitivity.json`, `bift_window_sensitivity.json`: numerical audit trail.

## Window selection

Search covered peak-containing contiguous intervals of at least seven frames within WT 150–195 and AA 155–215. All such windows with evaluable fits were retained in the audit table: 475 WT and 876 AA candidates. A candidate passed the declared presentation-quality filters when:

1. The apparent Rg span over qmin 0.010, 0.012, 0.015 and 0.018 Å⁻¹ was ≤5%.
2. The Guinier fit's weighted residual sum divided by n−2 was ≤1.5, using qmax Rg≤1 and at least eight positive points.
3. The mean squared standardized difference between the first and second halves, after one multiplicative intensity rescaling, was ≤1.5. Shared residual-background covariance was included in the half-comparison errors.

Among passing windows, maximize the norm of the mean-profile signal-to-error vector over q=0.012–0.25 Å⁻¹. This is a transparent signal-versus-consistency criterion, not a proof of purity or a search for a target Rg. There were 384 passing WT and 871 passing AA windows under the supplied error model. Adjacent candidates have almost the same score; for example WT 159–174 differs by only ~0.005% from 159–175. The recommended intervals sit on a broad plateau, so the precise final frame is not biologically meaningful. The selected windows are arithmetic averages, not EFA-rotated component curves.

Primary Guinier fit ranges are 0.012–0.038 Å⁻¹ (WT) and 0.012–0.030 Å⁻¹ (AA). Points outside these ranges are gray and retained for context, not silently discarded. Both fits have weighted residual mean squares below one, consistent with conservative supplied errors. The EFA report explains why smaller empirical errors could reveal weak additional clean-peak structure.

## P(r) method and stability

P(r) was calculated with the actual `bioxtasraw.BIFT.doBift` implementation from the recorded upstream commit, using a 40-point r grid, an initial alpha search from 150 to 1e10 with 12 values, and Dmax 50–220 Å with 15 initial values. RAW then optimizes the parameters; the initial grid limits are not hard final bounds. The primary q range is 0.012–0.245 Å⁻¹ (0.245 is the last recorded point below the requested 0.25 cutoff). The primary uncertainty sampling uses 150 parameter draws; sensitivity runs use 80. This is RAW's local alpha/Dmax uncertainty estimate, not a bootstrap of every experimental error source.

Across qmin 0.010–0.015 and qmax cutoffs 0.20–0.30 Å⁻¹, BIFT Rg spans 27.15–27.25 Å for WT and 34.42–34.69 Å for AA. Dmax spans ~93–94 Å and ~157–160 Å. Moving one window boundary by ±2 frames gives BIFT Rg ~27.10–27.21 Å for WT and ~34.30–34.64 Å for AA. Primary results lie within or close to these narrow ranges.

Guinier and P(r)-derived Rg are deliberately labeled separately. Their difference reflects finite-q/model sensitivity and must not be concealed by forcing agreement. Long-r features are not evidence of a separate species. Background choice is a larger systematic, especially for WT: the alternative early background used elsewhere in the original notebook produces larger low-q apparent sizes. These recommended plots are conditional on the primary notebook correction and do not establish absolute sample purity.

## Suggested figure caption

**SEC-SAXS characterization of Tau-5\* WT and W397A/W433A.** Scattering profiles were averaged over zero-indexed frames 159–175 for WT and 169–188 for W397A/W433A, selected for high signal-to-noise and consistency of profile shape and Guinier fits. (A) I(q) normalized by the corresponding Guinier I(0); points show averaged data with propagated uncertainties and lines show BIFT fits. (B) Bayesian indirect Fourier transforms, normalized by their BIFT I(0), with conditional uncertainty shading. (C,D) Guinier fits using qmax Rg≤1, with excluded q points shown in gray and standardized fit residuals below. Guinier Rg values are 25.8±0.6 Å and 32.3±0.5 Å; BIFT Rg values are 27.2 Å and 34.6 Å for WT and mutant, respectively. Error estimates exclude background-choice and other systematic effects. These ensemble-averaged profiles do not by themselves establish oligomeric state.

## Reproduction audit added after the new context

These are candidate presentation windows under the declared method, not final oligomer assignments. Notebook 01 reproduces Natalie's original Guinier **24.90 / 30.78 Å** and isolates the effect of shortening the q fit range. `Rg_reproduction_and_qrange.png` shows that comparison. Compare current P(r) Rg **27.2 / 34.6 Å** to her GNOM **27.3 / 34.23 Å**, not to her Guinier values. AA Dmax **159 vs 139 Å** is unresolved; the P(r) method and tail uncertainty must remain visible in any presentation.

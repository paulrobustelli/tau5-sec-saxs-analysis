# Checkpoint — paused at the user's request, 21 September 2026

## Scientific objective
Determine whether AA's larger apparent SAXS size can be explained by oligomer contamination, and obtain the least-contaminated WT–AA comparison. The objective is **not** to identify or isolate clean oligomeric species. No conclusion of intrinsic AA expansion or absence of oligomer contamination is justified yet.

## Completed and reviewable
- Original Natalie repository cloned without changing its files.
- Three new notebooks executed: 01 reproduction/correction audit, 02 window/profile/EFA analysis, 03 report EFA sensitivity. 19 executed code cells, no cell errors; embedded figures and numerical outputs included.
- Natalie's saved-profile Guinier results reproduced: WT 24.904 Å and AA 30.785 Å. AA's fit-range convention accounts for most of the change to the earlier ~32 Å low-q estimate.
- Existing presentation plots use arithmetic mean windows WT 159–175 / AA 169–188, not EFA-isolated species. These plots predate the independent baseline assessment.
- EFA did not reliably remove or exclude oligomer contamination. One dominant clean-peak signal does not exclude a weak co-eluting contaminant.
- Source poster GNOM P(r) Rg 27.3 / 34.23 Å must be distinguished from Guinier estimates; our BIFT P(r) results are a different regularization method. AA Dmax remains method dependent.

## Baseline provenance of the existing figures
Natalie's window choices were retained, then recalculated: WT constant residual background 220–236; AA linear endpoint windows 0–10 and 439–449 with original HDF5 error weights. Shared background uncertainty was propagated. Alternatives were tested, but this was not an independent baseline search.

## New independent-baseline work: PROVISIONAL, NOT VALIDATED
The user subsequently requested an independent best-practice subtraction. Work has started in `analysis/independent_baseline.py`, `analysis/inspect_independent_baseline.py`, and `independent_baseline/`. It is an assessment of the **residual** baseline only: supplied HDF5 data are already facility-subtracted. The user will obtain unsubtracted sample/buffer frames next time for a full subtraction from raw data.

A raw q-band inventory shows elevated post-peak low-q baselines for both samples, and later large disturbances in WT. This motivates comparing baseline mechanisms and excludes simply assuming that full-run endpoint interpolation is optimal. It does not establish capillary fouling rather than residual sample scattering.

The first exploratory search selected WT pre 55–94/post 200–229 and AA pre 35–74/post 225–254. The broad search zones were chosen by visual inspection of the raw trace; candidate windows were screened using q-band slope and split-profile tests. These selections are **not approved buffer regions**: absence of analyte is unproven, and the WT post interval may still include tail signal. Do not automatically promote the longest-passing-window heuristic to a best-practice baseline.

Pre-only, post-only, interpolated-linear, and an integral-shaped alternative were computed with a common q=0.012–0.030 Å^-1 Guinier range. The integral implementation is custom, uses smoothing/positive deposition constraints, and is **not an exact RAW reproduction**. It needs synthetic validation, review of clipping bias, and mechanistic justification before any use as a preferred correction. The provisional results show that WT size is particularly baseline sensitive; the apparent AA–WT difference cannot be treated as settled.

No new final figures or executed notebook have been made from this independent-baseline work. All files in `independent_baseline/` are exploratory. Do not replace the established reproduction plots with these values yet.

## Resume here
1. Inspect the raw q-band traces and candidate anchors; reject tails, spikes and later peaks. Test pre/post stationarity and profile similarity with sensitivity to conservative input uncertainties.
2. Decide which baseline mechanisms are supported (constant residual, instrumental drift, fouling-like offset). Keep uncertainty across plausible mechanisms when the transition through the protein peak is unidentifiable. Do not choose a model because it produces a desired Rg or smaller AA–WT difference.
3. Validate any custom numerical correction using synthetic known baselines. Propagate baseline covariance and assess original facility-subtraction covariance limitations.
4. Recompute leading/central/trailing AA and WT profiles with **identical q ranges**, inspect normalized low-q shape and apparent Rg across elution and across baseline choices. This is the outstanding direct contamination check. Its preliminary variable-range numbers were inspected but not accepted as a result.
5. Put the independent correction, alternatives, audit and contamination check into a new executed notebook; retain notebooks 01–03 as prior-method references. Recompute I(q), Guinier, P(r), EFA and plots only after the baseline choice is defensible.
6. When raw sample and buffer frames arrive, redo the original solvent subtraction with acquisition metadata and suitable blank validation.

## References for the next pass
- RAW baseline documentation: https://bioxtas-raw.readthedocs.io/en/v2.1.2/tutorial/s2_baseline.html
- RAW SEC workflow and the distinction between buffer subtraction and baseline correction: https://bioxtas-raw.readthedocs.io/en/v2.1.0/api/ex_sec_saxs.html
- Brookes et al., integral baseline/capillary fouling method: https://journals.iucr.org/j/issues/2016/05/00/vg5038/index.html
- IUCr reporting guidance: https://journals.iucr.org/d/issues/2017/09/00/jc5010/

Working folder: `/Users/f0044gk/Desktop/SAXS`. Portable local repository: `/Users/f0044gk/Desktop/SAXS/repository`. Private remote: https://github.com/paulrobustelli/tau5-sec-saxs-analysis.

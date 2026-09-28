# Tau-5* chemical-shift and SAXS snapshot

Open [Tau5_current_results.ipynb](Tau5_current_results.ipynb). The notebook contains executed plots and reruns with NumPy and Matplotlib using `snapshot_data.json`.

Includes full-sequence/R2 DSSP helicity versus D2D, experimental/raw/fitted secondary shifts, separate residual bar plots for every measured nucleus in **both sign conventions**, and current SAXS predictions versus the earlier joint fit.

This is an interim 2026-09-28 snapshot: latest fits use shifts only, with no entropy penalty; SAXS is a prediction for those weights. Earlier pooled joint fits are labeled separately. Unequal candidate counts and small effective ensemble sizes mean these are not final ~500-member selected ensembles or a convergence claim. GP is included structurally; terminal shifts are excluded.

Weights reference local candidate IDs; coordinate pools and licensed predictor executables are not included. Source scripts preserve methods and require the original workspace to run. The EOM results elsewhere in the repository are unchanged.

## Pool versus selected absolute residuals

The notebook now includes per-nucleus **|Exp − Calc|** bar plots for the uniform pools versus the selected continuous weights, for WT/AA and both replicates. Standalone PNG/PDF plots and numerical CSVs are in [CS_absolute_errors](CS_absolute_errors/). These preserve the original snapshot and do not incorporate unfinished larger pools.

# Paper Figures outputs

Start with [Paper Figures.ipynb](../Paper%20Figures.ipynb). It is executed and explains all selection rules and limits. `RECOMPUTE=True` repeats calculations; by default it reads cached results and regenerates figures.

Primary reference processing: WT buffer220–236, average162–178; AA buffer24–70, average168–188, zero-based inclusive. These are documented reference choices, not a uniquely established optimal baseline. The early and linear alternatives remain part of the result.

WT/AA extended Rg:26.41/32.17Å; conventional25.48/32.25Å; P(r)26.91/35.46Å, conditional on primary processing. P(r) tails are notably baseline sensitive. WT early24–70 produces a long-distance tail and P(r) Rg40.69Å; post and linear give26.91/27.50Å. Do not hide this uncertainty by constraining Dmax.

WT has a reproducible second scattering pattern, but no supported shoulder Guinier Rg or established oligomer identity. AA does not show a second mode above the conditional rank-one noise controls. Primary AA comparisons use frame averaging, not forced EFA purification.

PNG/PDF figures01–08, all curves,100source-frame perturbations per sample,60rank-one noise simulations per sample/background,108support/q-range rotations, window-selection tables and BIFT sensitivity results are included. Normalized I(q)/Kratky bands propagate joint fitted normalization uncertainty. All noise intervals are conditional and exclude baseline-model and biological variability.

Scripts: `analysis/paper_figures.py` and `analysis/paper_figures_controls.py`. Uses the original input HDF5 profiles, subtracting one buffer or one interpolated background. No clipping or smoothing of observed curves. A108case rotation grid is exploratory, not a confidence interval.

# EOM: WT versus AA
Analysis of completed ATSAS-online jobs WT120EOM and AA120EOM, downloaded 2026-09-25.

Inputs: WT frames 162–178 with buffer 220–236; AA frames 168–188 with buffer 24–70 (zero-based, inclusive). Both use the 120-residue constructs. Each run has a 10,000-conformer RANCH pool; GAJOE uses 100 cycles, 1,000 generations, 50 ensembles, and exactly 50 non-repeated conformers per selected ensemble. No additive constant was fitted.

The returned fit covers **0.012–0.50 Å⁻¹**, although submitted data extend to 0.60 Å⁻¹. These plots show the returned range, without extrapolation.

The server's selected distribution and its reported average summarize selection across cycles; the 50 structures listed in the final log describe the best ensemble. Keep these distinct. Histogram bin coordinates give slightly different averages from the server's unbinned header values.

## Interpretation and limits
The best ensembles have mean Rg approximately **26.7 Å (WT)** and **32.0 Å (AA)**. The server selected-distribution means are **27.40 Å** and **31.68 Å**. AA therefore favors larger conformations under this common EOM modeling protocol.

The ± values above describe **conformational spread**, not uncertainty in the mean, not a confidence interval, and not uncertainty from buffer subtraction. The best-ensemble RMS sizes, sqrt(mean(Rg²)), are also reported: they are the more relevant ensemble quantity for the low-q limit when conformers have equal I(0). Model-coordinate Rg and experimental scattering Rg need not match exactly because hydration and scattering contrast enter the latter.

These are single-chain ensembles by construction. Fits are compatible with a monomer ensemble, but do not independently determine molecular weight or rule out oligomer contamination. RANCH pools and selected populations are model-dependent; an apparent multimodal histogram is not proof of discrete physical states. One random pool per sample was used; the 100 optimizer cycles do not quantify independent-pool or experimental uncertainty.

The reported χ² is retained verbatim. The independently computed mean squared standardized residual uses N in the denominator and is listed separately in summary.json; it is not assigned a formal number of degrees of freedom. Error bars preserve the supplied errors. Shared buffer errors can correlate data; the diagonal residual analysis does not model full covariance. AA has structured residuals: the fit overpredicts much of q ≈ 0.23–0.36 Å⁻¹ and underpredicts the high-q tail. Its lag-one residual correlation is about 0.55 (WT: approximately 0.00). This flags systematic mismatch that may reflect the conformer pool, scattering model, or background; it does not identify oligomers. Inspect residual structure rather than using χ² alone.

## Reproducibility
This notebook reads preserved server tables and logs, recomputes ensemble summaries and residuals, and exports the figures. `selected_conformers.zip` contains both selected 50-model ensembles. `data/` includes full 10,000-member size tables, fit/distribution outputs, sequences, and submitted profiles. Full downloaded archives (including every starting-pool structure) are retained locally under `outputs/ATSAS_models/eom/`; their hashes are recorded in provenance.json. No smoothing or refitting is applied here.

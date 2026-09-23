# DENSS exploratory envelopes

DENSS 1.8.7, standard SLOW mode, 64³ grid, oversampling 3, 20 independent random starts per sample, 4 workers; default positivity, connectivity and shrinkwrap; no imposed symmetry or atomic reference. Input is the final buffer-subtracted frame average used in Paper Figures, retaining measured errors and q=0.012–0.25 Å⁻¹. WT: frames162–178, buffer220–236; AA: frames168–188, buffer24–70 (zero-based inclusive). Dmax is fixed to the existing BIFT values: WT91.7258Å, AA166.7131Å. DENSS internally regularizes/interpolates the input using Sasrec. Original inputs, fitted profiles and reconstruction logs are included.

These are **effective scattering envelopes for disordered proteins**, not unique molecular structures, atomic models, or recovered conformational ensembles. In general the ensemble average of squared scattering amplitudes is not the squared amplitude of an ensemble-average density. Averaging here is over independent reconstructions, not protein conformations. Detailed lobes, cavities and apparent topology are not established. Dmax, baseline and q-range sensitivity are not tested by repeated random starts. This is a low-q exploratory reconstruction, not a full-q DENSS analysis; DENSS guidance recommends using the full reliable scattering range for final reconstruction. AA's long P(r) tail is sensitive to q-range, so its envelope is conditional on the primary P(r) choice.

Surfaces enclose90%of positive density. Both aligned-average and subsequently refined maps are displayed, with equal spatial axes and no extra display smoothing. This is an explicit visualization threshold, not a molecular boundary or confidence surface. Density uses DENSS's arbitrary default electron normalization; these maps do not estimate molecular weight. DENSS's FSC diagnostic compares reconstructions to a common reference; it is not an independent experimental resolution claim.

Post-averaging refinement was run from each aligned-average map against the same input, with seed20260923 and default DENSS refinement settings. It restores scattering agreement by modifying the density; the refined result is a single effective map, not an average over conformations.

The fit-check plot directly compares back-calculated average-map scattering with the measured data, fitting scale only, no offset. Individual fits are also shown. Averaging can worsen agreement, so the averaged map should not be called a validated single structure even if individual reconstructions fit well.

## Reproduction
Install requirements.txt in an isolated environment, then run `python run_reconstructions.py`. Then run `python render_denss.py` , `python render_denss.py --refined`, and `python denss_report.py` in the same directory to regenerate figures and the notebook. This repeats20random starts; exact seeds from this run are recorded in the archived logs. The complete-run archive includes all individual and aligned maps. Open WT_avg.mrc and AA_avg.mrc in ChimeraX/PyMOL or use the self-contained interactive HTML.

References: [DENSS documentation](https://github.com/tdgrant1/denss), [Grant2018](https://doi.org/10.1038/nmeth.4581).

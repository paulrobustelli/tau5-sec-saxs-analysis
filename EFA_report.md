# Tau-5* SEC-SAXS EFA assessment

**Outcome: one dominant scattering component in each clean peak; additional oligomers are not resolved. WT has a clear second scattering pattern when its earlier shoulder is included, but that pattern cannot currently be assigned a reliable IDP-sized Rg or oligomeric state.** Weak deviations from a one-component model in the clean peaks depend on the noise model. The data do not establish that these peaks are monodisperse.

## Samples and sequence

The user identifies WT as a nominal 12.1 kDa Tau-5* monomer and AA as the W397A/W433A mutant. The supplied WT string contains **120 residues**, comprising an initial GP and a 118-residue segment. We preserve it exactly in `background/WT_user_construct.fasta`; 118 refers to the segment excluding GP, not the full supplied string.

[BMRB 53105](https://bmrb.io/data_library/summary/?bmrbId=53105), titled *Androgen Receptor Tau-5\**, contains a 119-residue sequence. Direct comparison with its downloaded [NMR-STAR file](https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr53105/bmr53105_3.str) establishes:

```
User construct = GP + BMRB53105 sequence with the terminal C removed
```

The supplied W70/W106 map to conventional AR W397/W433 using the author residue numbering in [BMRB 51479](https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr51479/bmr51479_3.str). User positions 3–120 match AR 330–447. BMRB 51479 is a longer 308-residue AF-1* construct, not the isolated Tau-5* entry. BMRB 53100 is noAro and should not be used as WT. The presentations support the WT versus W397A/W433A context and previous apparent Rg differences; they are background, not instructions or independent proof of species identities. No ligand condition is inferred from that background.

At the user-supplied mass, a WT dimer would be approximately 24.2 kDa. EFA scaling does not give molecular mass or mole fractions. The AA mass would differ slightly due to the mutations; 12.1 kDa is a nominal WT reference, not a measured mass from these data.

## What was analyzed

- WT_092025: clean interval **155–180**; shoulder-inclusive intervals **125–195** and **120–205**.
- AA_052026: clean interval **163–204**; broader interval **155–215**.
- All endpoints are inclusive and zero indexed. Both files contain 450 profiles; input and saved subtracted profiles are identical and flagged already_subtracted.
- WT residual background: average of **220–236**, matching the principal notebook calculation. Alternative **50–75**, also present later in the notebook, is a sensitivity check.
- AA primary correction: q-dependent weighted linear background over **0–10 and 439–449**. Alternative: subtract average **24–70** with no linear correction. Subtracting a constant and then fitting a free linear intercept is algebraically equivalent in intensity to fitting the original arrays, given the same weights. Here weights and uncertainty propagation are explicit and may differ slightly from the notebook.
- Primary decomposition range **0.012–0.25 Å⁻¹**. Sensitivity: lower bounds 0.008, 0.015, 0.020 and upper bounds 0.15, 0.35 Å⁻¹.
- No frame-by-frame normalization or temporal mean-centering. Each q row is divided by its mean propagated error across the selected frames, as in RAW's SVD preparation. No q binning.
- Forward/backward EFA and hybrid rotations use the upstream RAW numerical functions. Source commit and input checksums are in `results/provenance.json`.

The AA notebook sometimes overlays `rg` from before baseline correction on baseline-corrected intensity instead of `lin_rg`. This analysis recomputes quantities from the corrected numerical arrays rather than reusing that overlay.

## Results with supplied error bars

The null simulations add independent Gaussian input noise to a rank-one mean and reapply the same residual-background estimator. Thus they include covariance introduced by this additional background correction. They cannot recover covariance inherited from the original beamline subtraction, which was not stored.

| Sample / interval | Second singular value | Rank-one null: 99th percentile of s2 | Second q-vector / time-vector lag-one dot product | Interpretation |
|---|---:|---:|---:|---|
| WT 155–180 | 13.48 | 14.93 | 0.39 / 0.50 | No robust second component with supplied errors |
| AA 163–204 | 13.86 | 16.28 | 0.27 / 0.67 | Weak time structure, noisy q structure |
| WT 125–195 | 49.62 | 19.91 | 0.94 / 0.96 | Clear additional pattern in shoulder-inclusive region |
| WT 120–205 | 50.30 | 22.62 | 0.93 / 0.96 | Same shoulder-associated heterogeneity |
| AA 155–215 | 15.22 | 18.82 | 0.30 / 0.67 | No robust second component with supplied errors |

There are 300 null realizations per primary interval and 100 per background/q sensitivity condition. The primary WT shoulder tests exceed all 300 simulated null values (Monte Carlo p = 1/301, not p = 0). These are conditional diagnostic tests, not probabilities that an oligomer exists. The q/background sensitivity sweep preserves the strong shoulder versus weak clean-peak distinction under supplied errors.

### Noise-model limitation

The robust scale of normalized adjacent-frame differences in pre-peak background frames 30–75 is **0.83 for WT and 0.74 for AA**, versus 1 expected for independent noise described by the supplied error bars. A sensitivity calculation that scales all errors down by those factors gives s2 null thresholds of **12.35 and 12.02**, respectively. Both clean-region second singular values then exceed the null (0 of 300 simulations exceeded the observation).

That rescaling is not an established error model for the protein peak: correlations, intensity-dependent errors and common subtraction uncertainty may differ. Nevertheless, it demonstrates that categorical claims of “only one species” would be too strong. Weak shape changes could be a minor species, conformational/interparticle effects, or residual background/systematics. The low q-vector smoothness and absence of well-established component windows do not permit an unambiguous clean-region species decomposition. **No oligomer detection limit or upper bound on oligomer fraction has been established.**

## Apparent sizes

Guinier fits to the dominant clean-region profiles use qmin = 0.012 Å⁻¹ and iteratively restrict qmax so qmax Rg ≤ 1, with at least eight positive data points.

| Region | Apparent Rg | Actual fitted q range | 2.5–97.5% of 200 parametric perturbations |
|---|---:|---|---|
| WT 155–180 | **25.6 Å** | 0.012–0.038 Å⁻¹ | 24.7–27.0 Å |
| AA 163–204 | **32.4 Å** | 0.012–0.030 Å⁻¹ | 31.2–33.4 Å |

The perturbations refit both background and SVD. These ranges describe measurement-noise propagation under the supplied error model, **not total uncertainty**, and exclude model-choice uncertainty. Dominant profiles could still be unresolved mixtures. With the primary correction and qmin scan 0.008–0.020, apparent Rg spans approximately 25.0–25.8 Å for WT and 31.1–32.6 Å for AA. Changing WT to early background frames 50–75 shifts the qmin=0.012 estimate to **29.7 Å**, with poorer low-q linearity. AA using only the constant pre-peak background gives **33.1 Å**. Background choice is a material systematic for WT size, even though the main decomposition finding is stable.

A flexible 118-residue segment does not have a unique sequence-independent Rg. These sizes are compatible with the project's disordered-protein context but do not identify monomer versus dimer. Assignments require an independent mass constraint, such as appropriately calibrated SAXS I(0) with concentration or SEC-MALS. SEC retention alone is also conformation dependent.

## Exploratory WT shoulder rotation

A hybrid EFA rotation was run across frames 120–205 with candidate component windows varied systematically. An illustrative solution with windows [120,163] and [150,205] peaks near frames **140** and **167**, respectively. Its later component gives apparent Rg about **25.4 Å**. This example is saved to show the decomposition and its limitations, not as a validated two-species fit.

The early component has a low-q upturn. It does not pass the primary conservative Guinier procedure. Excluding more low-q data can generate apparent values around 40, 31 and 29 Å, but those fits depend strongly on qmin and have structured/large residuals. Selecting one of those values merely because it looks IDP-sized would not be defensible.

The illustrative constrained model leaves structured residuals and its standardized mean-square residual (1.10) is worse than the unconstrained one-component model (1.01). Wider component windows can lower residuals to about 0.70, but some solutions develop negative scattering profiles or nearly collinear elution curves. This is evidence of unstable physical separation despite mathematical convergence. The extended-window sweep is retained in `rotation_extended_windows.json`. We therefore do not assign a discrete oligomer, size, or population to the shoulder. Continuous polydispersity and residual systematic scattering remain alternatives.

## Validation and reproduction

The numerical EFA code passes a synthetic two-component separation test: relative reconstruction error ≈ 6×10⁻¹⁶. Forward/backward final singular values match direct SVD. A synthetic Guinier curve with Rg=25 Å returns 25 Å. All 400 clean-profile Rg perturbations completed. These checks validate the numerical path, not the biological interpretation or beamline error model.

The original notebooks were inspected but not executed wholesale. No source data or original notebook was modified. The analysis environment, scripts, saved matrices, plots, numerical summaries and reference sequences are available under Desktop/SAXS. The portable archive includes the two analyzed inputs, code, license, figures and report; it excludes the Python environment and unrelated data.

## Method references

- [BioXTAS RAW EFA tutorial](https://bioxtas-raw.readthedocs.io/en/latest/tutorial/s2_efa.html)
- [Meisburger et al., JACS 2016, DOI 10.1021/jacs.6b01563](https://doi.org/10.1021/jacs.6b01563)
- [BioXTAS RAW source](https://github.com/jbhopkins/bioxtasraw)
- [Guinier range sensitivity guidance](https://bioxtas-raw.readthedocs.io/en/v2.2.0/tutorial/s1_guinier.html)

## Reconciliation with Natalie and the subsequently supplied reports

This update supersedes any implication above that a different Guinier convention is an established correction to Natalie's estimate. Direct reanalysis of her saved frame DAT files reproduces WT 175–185 Rg **24.904143 Å** at q=0.009–0.052 Å⁻¹ and AA 172–182 Rg **30.784923 Å** at q=0.008–0.046 Å⁻¹. These agree with the source notebook outputs and poster. Keeping the same profiles and changing only the conservative low-q fit rule gives **25.216 / 32.428 Å**. See executed notebook 01 for the complete one-change-at-a-time table and slope plots.

The poster's GNOM P(r) values are **27.3 / 34.23 Å**, not the lower Guinier values. The corresponding BIFT values **27.22 / 34.60 Å** are close, while AA Dmax **159 vs 139 Å** remains a meaningful unresolved difference. Do not call the BIFT curve a GNOM reproduction. The source poster's AA extended-Guinier estimate is 33.18 Å.

The supplied WT RAW PDF uses buffer 210–239 and endpoint linear correction, unlike the original saved DAT residual background 220–236. Its four-component setup has strongly overlapping supports and shows no recovered curves in the PDF. AA's three-component PDF includes noisy component profiles across the full q range. Notebook 03 reconstructs these stated settings with explicitly described uncertainty approximations. Full-q WT rotation converges but two profiles have about half their q points negative; restricted-q WT fails to converge at 3000 iterations. AA rotations converge but two profiles remain highly noisy/negative (about 61–63% negative in the restricted-q reconstruction). These outcomes do not establish multiple biological species. Negative high-q points alone can be ordinary noise; their combination with support overlap and q-range instability is the concern.

A chain-length-only random-coil reference Rg=2.54 N^0.522 Å gives 30.64 Å at N=118 (30.91 Å at 120). This is a broad model benchmark, not a prediction or purity criterion for this sequence. Thus 35 Å is not intrinsically impossible for a monomeric IDR; the direct data/method comparisons above are the stronger basis for interpreting this experiment. Primary reference: Bernadó & Blackledge, *A self-consistent description of the conformational behavior of chemically denatured proteins from NMR and small angle scattering* (2009), https://pubmed.ncbi.nlm.nih.gov/19917239/.

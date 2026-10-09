# Tau-5 bAIes/OpenMM: completed 10 ns validation runs

This directory preserves the completed CPU validation trajectories for the
120-residue GP–Tau-5* WT and W397A/W433A (`AA`) constructs. Each trajectory has
1,000 frames saved every 10 ps over 10 ns. These short trajectories validate
the implementation and restart path; they are not converged production bAIes
ensembles. The bAIes paper-scale sampling target is 2 microseconds.

## Descriptive results

The table below excludes the first 0.2 ns when reporting the descriptive
averages. The exclusion is a simple startup convention, not a convergence
claim.

| Construct | Rg (A), mean ± SD | DSSP alpha fraction | R2 AR395–405 alpha | R3 AR430–440 alpha |
|---|---:|---:|---:|---:|
| WT | 42.61 ± 7.68 | 0.00024 | 0.00222 | 0.00000 |
| AA | 40.19 ± 7.03 | 0.00042 | 0.00260 | 0.00028 |

The near-zero helicity is the result of these 10 ns trajectories, not evidence
that the experimental helices are absent. Ten nanoseconds is far too short to
assess equilibrium secondary-structure populations for these IDPs. The raw
size ordering in this validation sample is also not treated as a converged
prediction.

## Distogram sanity check

The original AlphaFold rank-1 pickle was audited from logits through the final
OpenMM atom pairs.

- Both pickles contain `120 × 120 × 64` logits and 63 finite bin edges.
- The WT and AA logits are exactly symmetric across residue-pair indices.
- The distogram dimension, AlphaFold PDB sequence, and OpenMM topology sequence
  all contain the same 120 residues.
- Construct positions 70 and 106 are W/W in WT and A/A in AA, corresponding to
  AR397 and AR433 after the two-residue GP extension.
- OpenMM atom serials are contiguous, so the preprocessing script's 1-based
  PDB serials map exactly to OpenMM indices after subtracting one.
- Reprocessing with the residue-specific cutoff matrix and sequence separation
  greater than three regenerates exactly 19 WT and 15 AA atom pairs.
- Every regenerated atom pair matches the stored `baies_params.dat` file.
  Gaussian means and widths reproduce within `5 × 10^-7 nm`, the expected
  six-decimal text-rounding error.
- Recomputing the softmax in float64 changes probability normalization error
  from at most `9.77 × 10^-4` to below `6 × 10^-16`; it does not change any
  selected pair. This small float16 rounding effect is documented but does not
  indicate an indexing or distogram-selection error.
- For finite-mode residue pairs, the distogram modes correlate with distances
  in their matching AlphaFold coordinates: `r=0.839` (WT) and `r=0.772` (AA),
  with RMSE 0.170 and 0.196 nm, respectively. This is an independent check that
  each pickle is paired with the intended structure.

The audit therefore supports that the saved distograms were read and mapped
correctly. Numerical evidence is in `distogram_audit_summary.csv`,
`distogram_restraint_audit.csv`, and `distogram_audit.json`.

## Files

- `Raw/WT` and `Raw/AA`: full XTC, topology, final portable XML state, binary
  checkpoint, log, fitted restraints, atom list, and AlphaFold prediction.
- `summary_10ns.csv`: ensemble-level Rg and helicity summaries.
- `timeseries_10ns.csv`: frame-level Rg and helicity.
- `per_residue_helicity_10ns.csv`: DSSP propensities for every residue.
- `Rg_and_helicity_10ns.*`: time series and full-sequence helicity profiles.
- `Rg_distribution_and_restraint_distances.*`: Rg distributions and restrained
  pair distances compared with fitted distogram means.
- `../code/analyze_completed_and_audit.py`: complete analysis and audit code.

Use `final_state.xml` to continue on another machine. The `.chk` files are less
portable because OpenMM checkpoints depend on compatible software, platform,
and hardware.

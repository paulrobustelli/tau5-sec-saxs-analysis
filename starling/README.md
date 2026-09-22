# Tau-5 STARLING predictions

Executed in the official STARLING Colab on a Tesla T4 using idptools-starling
2.0.2 and torch 2.11.0+cu128. Four independent 400-conformation ensembles, DDIM,
30 diffusion steps, 188 mM ionic strength, no experimental restraints. The salt
setting follows the existing Colab execution; it is not independently verified
as the SAXS buffer. WT supplied sequence has 120 residues; 118 removes initial GP.
AA replaces both W with A (W397A/W433A). Seeds, in execution order:
WT120 20260922; AA120 20260923; WT118 20260924; AA118 20260925.

| Construct | Mean Rg (Å) | RMS Rg (Å) | Conformational SD (Å) |
|---|---:|---:|---:|
| WT118 | 30.95 | 31.70 | 6.86 |
| AA118 | 32.47 | 33.13 | 6.60 |
| WT120 | 31.03 | 31.72 | 6.60 |
| AA120 | 31.90 | 32.57 | 6.57 |

The histograms are exact counts in common 2 Å bins extracted from the generated
per-conformation arrays, not Gaussian approximations. The local summary preserves
full-precision moments and counts; raw distance-map ensembles and individual Rg
arrays were saved in /content/tau5_starling in the Colab runtime. Browser download
could not be captured here; preserve those runtime files before disconnecting if
full conformations are needed. Histogram counts each sum to 400.

These are sampled single-chain predictions, not experimental oligomer fractions.
Histogram width measures conformational heterogeneity, not uncertainty in the
predicted mean. Finite-sampling standard errors of means are roughly 0.33 Å;
model error is separate. Small 118/120 differences should not be overinterpreted.
RMS Rg is closer to the ensemble averaging relevant to SAXS than the arithmetic
mean, but these are coarse-grained distance-map radii, not explicit back-calculated
X-ray scattering fits with solvent contrast. STARLING predicts modest AA expansion,
less than the current experimental difference; it does not establish agreement
with the absolute experimental WT size or validate the experimental baseline.

Official notebook: https://colab.research.google.com/github/idptools/idpcolab/blob/main/STARLING/STARLING_demo.ipynb
Method: https://doi.org/10.1038/s41586-026-10141-2

# Extended Guinier fits: Zheng & Best (2018)

Reference: https://doi.org/10.1016/j.jmb.2018.03.007, equations 5 and 6.

The model is ln(I/I0) = -(q Rg)^2/3 + 0.0479 (nu - 0.212) (q Rg)^4.
Rg = sqrt[gamma(gamma+1)/(2(gamma+2nu)(gamma+2nu+1))] b N^nu,
with gamma=1.1615, b=5.5 Angstrom, and N=number of residues minus one.
The free parameters are I0 and nu; Rg follows from nu. This is the paper's
coupled single-chain model, not an unconstrained polynomial fit.

We use the supplied 120-residue construct (119 peptide bonds), with a
118-residue sensitivity calculation. WT and AA have the same length.
Weighted least squares is performed in intensity space with the supplied
diagonal errors; unlike a log-intensity fit this can retain negative observations.
The objective implementation differs from a weighted log-space regression.
No smoothing, positivity clipping, constant-background fitting or Rg target is used.

For each curve, scan every measured upper q endpoint through 0.12 Angstrom^-1
and retain the largest endpoint with fitted qmax Rg <= 2. Primary qmin=0.012;
sensitivity qmin=0.010, 0.015 and 0.020. At least eight points and six >3 sigma
positive intensities are required. The numerical nu bounds are 0.2–0.8; boundary
solutions are rejected. These are numerical guards, not validation of all nu values.
All scans are retained, so fit-range dependence is visible. The selected maximum
range is not optimized for chi-square. Input baselines/decomposition are unchanged.

| Curve | Standard Guinier, q=0.012–0.030 | Extended Rg | Extended qmax |
|---|---:|---:|---:|
| WT average | 28.60 | 26.49 | 0.074 |
| WT main candidate | 27.42 | 25.75 | 0.076 |
| AA average | 32.15 | 31.73 | 0.062 |
| AA main candidate | 31.80 | 31.23 | 0.064 |
| WT earlier component | unstable | 40.87 (exploratory) | 0.048 |
| AA secondary component | unsupported | no supported fit | — |

Rg units are Angstrom; q units are inverse Angstrom. Broader fitting windows and
model assumptions change these values; this is not a replacement of the earlier
measurements. Although the curvature term raises Rg at a fixed fit interval,
using a broader interval can yield a smaller fitted Rg, as seen here.

WT shoulder: fitted Rg spans 37.56–43.26 Angstrom over qmin choices; primary
reduced chi-square is 2.44 and the lowest-q excess is not fully described.
The assumed 120-residue single-chain relation is not validated for an oligomer
or mixture. This result is a model diagnostic, not identification or reliable
mass/stoichiometry determination. AA's secondary curve remains unsupported.
IDP curvature by itself is not evidence of contamination.

The 100 stored full-decomposition perturbations are refitted for average and main
component at the fixed selected window. Their percentile intervals are conditional
on the existing baseline/support/model and do not incorporate model error or
fit-window selection uncertainty. Original facility covariance is unknown. Fit
weights are diagonal even though EFA can correlate q points. Shoulder errors
condition on fixed elution profiles; no full shoulder bootstrap is available.
Small fit standard errors must not be treated as total uncertainties.

Reproduce with `python extended_guinier.py` or execute the accompanying notebook.
The input NPZ files are copies of the prior decomposition audit, not new data.
Synthetic exact-recovery validation checks the formula, units and chain length.

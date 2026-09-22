# AA direct frame-window comparison

Windows fixed before fitting: leading 160–169, center 170–179, trailing 180–189, late tail 190–199 (zero-based inclusive). These are equal-size 10-frame arithmetic averages, without EFA or smoothing. They test the proposed leading-edge selection rather than optimizing windows for low Rg.

Same anchor frames as previous component analysis: pre 35–74 and post 225–254. Main correction interpolates anchor means; pre-only and post-only are sensitivity checks. Background sharing is propagated by explicit linear weights; paired comparisons preserve cross-window covariance. Unknown upstream correlations are not included.

| Window | Conventional Rg (Å) | Extended Rg (Å) |
|---|---:|---:|
| Leading 160–169 | 31.87 | 31.80 |
| Center 170–179 | 32.47 | 32.20 |
| Trailing 180–189 | 31.81 | 31.96 |
| Late tail 190–199 | 31.30 | 31.87 |

Conventional fits use q=0.012–0.030; extended Zheng–Best fits use a common q=0.012–0.055 Å^-1 and 120 residues (119 bonds). All extended primary fits have qmax Rg <2. Adaptive extended fits are also retained in JSON. Extended results differ from the earlier 31.7 Å average because both the window and upper fitting q differ.

Leading-minus-center extended Rg is about -0.40 Å; the 300-perturbation paired conditional 95% interval is -1.17 to +0.34 Å. Every tested baseline yields an interval including zero. Shape comparisons over q=0.012–0.15 allow one intensity scale and propagate shared-background covariance. Reduced difference chi-squares for leading/trailing/tail versus center are 0.66/0.78/1.01 for interpolation. There is no resolved shape difference under the supplied error model. This is not a proof of equality or purity.

The leading-edge average has only about 28% of the center scattering amplitude, so selecting it does not improve counting precision. There is no evidence here that discarding the center yields a cleaner or smaller population. Retain the central high-signal region as the primary measurement; use the leading edge as a consistency check. All DAT files contain q, I and propagated sigma.

Provenance correction: input and subtracted groups are identical for all 450 profiles in both source HDF5 files. S_ names and near-zero backgrounds suggest prior subtraction, but the histories are empty. Earlier statements that facility subtraction was definitively documented should be read as inference, not established provenance. No separate original sample/buffer export was identified in the current source repository; the generic RAW input-group label does not settle the upstream history.

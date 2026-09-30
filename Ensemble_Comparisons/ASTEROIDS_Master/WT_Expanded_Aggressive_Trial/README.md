# WT expanded aggressive-fit trial

All selected ensembles contain 500 distinct structures at equal weight. Shift scales are 0.05 ppm for all carbons, 0.25 ppm for N, and 0.05 ppm for HN/HA.

SAXS curves are multiplicatively normalized to the experimental Guinier-extrapolated I(0); no additive offset or free post-selection scale is used. Low q is q<0.04 A^-1, the unchanged middle is 0.04-0.10 A^-1, and the tail is q>=0.10 A^-1.

| Ensemble | C RMSD | CA RMSD | N RMSD | NH RMSD | low-q chi2/point | all-q chi2/point | Guinier Rg (A) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Raw pool | 0.496 | 0.376 | 0.787 | 0.217 | 26.137 | 72.190 | 32.90 |
| Current fit | 0.417 | 0.303 | 0.773 | 0.218 | 13.050 | 31.186 | 31.06 |
| 2x low / 0.5x tail | 0.403 | 0.285 | 0.789 | 0.218 | 3.607 | 24.241 | 28.94 |
| 4x low / 0.25x tail | 0.410 | 0.290 | 0.788 | 0.217 | 2.488 | 24.149 | 28.48 |

Experimental diagnostic Guinier Rg: 27.82 A.

The aggressive schedules improve C and CA only slightly beyond the current fit, do not improve N or NH, and leave the selected DSSP-H profile essentially unchanged. The 4x schedule gives the best low-q SAXS agreement, but sacrifices tail agreement. These results should not yet be propagated to every pool.

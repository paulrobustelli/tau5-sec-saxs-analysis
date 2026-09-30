# AA expanded aggressive-fit trial

All selected ensembles contain 500 distinct structures at equal weight. Shift scales are 0.05 ppm for all carbons, 0.25 ppm for N, and 0.05 ppm for HN/HA.

SAXS curves are multiplicatively normalized to the experimental Guinier-extrapolated I(0); no additive offset or free post-selection scale is used. Low q is q<0.04 A^-1, the unchanged middle is 0.04-0.10 A^-1, and the tail is q>=0.10 A^-1.

| Ensemble | C RMSD | CA RMSD | N RMSD | NH RMSD | low-q chi2/point | all-q chi2/point | Guinier Rg (A) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Raw pool | 0.489 | 0.392 | 1.354 | 0.165 | 1.823 | 9.182 | 32.30 |
| Current fit | 0.396 | 0.307 | 1.338 | 0.156 | 18.192 | 18.317 | 31.16 |
| 2x low / 0.5x tail | 0.383 | 0.295 | 1.319 | 0.154 | 1.835 | 7.992 | 32.30 |
| 4x low / 0.25x tail | 0.383 | 0.293 | 1.340 | 0.153 | 1.048 | 9.298 | 32.49 |

Experimental diagnostic Guinier Rg: 33.02 A.

The aggressive schedules improve C and CA only slightly beyond the current fit, do not materially improve N or NH, and leave the selected DSSP-H profile essentially unchanged. The 4x schedule gives the best low-q SAXS agreement, but sacrifices tail agreement. These results should not yet be propagated to every pool.

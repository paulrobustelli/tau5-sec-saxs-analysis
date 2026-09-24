# DAMMIF and EOM — EOM submitted

WT120EOM and AA120EOM submissions were confirmed by ATSAS online. Results are pending. DAMMIF is not yet submitted: GNOM preprocessing requires local ATSAS, and the official macOS browser download stalled. See submission_status.json.

Inputs are the final WT/AA frame-averaged curves, with original propagated errors. WT uses buffer frames 220–236 and averaging frames 162–178; AA uses buffer 24–70 and averaging frames 168–188. Frame numbering is zero-based, inclusive. The supplied WT construct has 120 residues including the initial GP; the AA sequence replaces the two tryptophans with alanines. See inputs/manifest.json.

Planned DAMMIF: GNOM P(r) from 0.012–0.25 Å^-1, Dmax sensitivity around the supported range; 20 independent slow-mode P1 reconstructions per sample, no imposed anisometry; alignment/clustering and direct fit validation of the representative/refined result. GNOM input will be generated with ATSAS; a BIFT curve will not be mislabeled as GNOM output.

Planned EOM: fully disordered monomer pools, no rigid domains, at least 10,000 conformations per sample; repeated selections from multiple pools/seeds. Use 50 conformations per selected ensemble without repeats where the server supports these short-chain settings. Compare the selected and starting Rg distributions, report mean Rg and sqrt(mean(Rg^2)) separately, plus fits/residuals. Assess full-q fitting against low-q and baseline sensitivity; assess effective pool diversity. Distribution widths are model-dependent, not experimental confidence intervals. This tests the monomer ensemble hypothesis, not oligomer stoichiometry.

Official services: https://www.embl-hamburg.de/biosaxs/atsas-online/
EOM guidance: https://www.embl-hamburg.de/biosaxs/eom.html

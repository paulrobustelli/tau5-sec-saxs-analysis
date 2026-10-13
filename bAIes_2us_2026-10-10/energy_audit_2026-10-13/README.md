# bAIes OpenMM energy and force audit — 13 October 2026

Actual reference-code tests support the core implementation. They do not establish an exact reproduction of the paper protocol, nor explain or validate the low-helicity Tau-5 ensembles.

## Reference implementations used

Built LAMMPS 2 Aug 2023 from its official source tarball with the authors’ fix_cmap.cpp patch. Built PLUMED v2.10 with the ISDB module including BAIES. Tests execute those binaries; these are not only comparisons against another Python transcription.

The complete coil-system test uses the authors’ PaaA2 tutorial topology, original coordinates, coefficients and CMAP file. Three configurations are compared: original, deterministic 0.01 Å coordinate perturbations, and 1% isotropic expansion. Atom order and coordinates are shared between engines. The OpenMM production runner is reused; minimization is disabled and the Reference platform is selected. Bond constraints are disabled for this comparison because LAMMPS’s reference system retains all bond potentials. Production used hydrogen-bond constraints, an intentional protocol difference.

The BAIES bias is tested with the actual WT/AA restraint inputs on five saved production conformations per construct. PLUMED driver supplies both energies and Cartesian bias forces. Both OpenMM Reference and M5 Max OpenCL single precision are evaluated.

## Results

| Component | Finding |
|---|---|
| Bonds and angles | Energies agree within 4e-11 kJ/mol. Combined-force RMS differences below 3e-10 kJ/mol/nm. |
| CMAP | Energy differences below 1.1e-8 kJ/mol; force RMS differences below 2e-7 kJ/mol/nm. This directly validates units, map indexing and orientation for this benchmark. |
| Bayesian bias, Reference | Largest energy difference 1.5e-10 kJ/mol; largest Cartesian force-component difference 2.9e-10 kJ/mol/nm. |
| Bayesian bias, GPU single precision | Largest energy difference 5.8e-6 kJ/mol; largest force-component difference 3.2e-5 kJ/mol/nm. |
| Torsions | OpenMM energies exceed LAMMPS by approximately 391.993 kJ/mol in all three configurations. The offset changes by about 0.0004 kJ/mol under perturbation. Force RMS differences are 0.0112–0.0116 kJ/mol/nm, versus reference torsion-force RMS 72–86 kJ/mol/nm: roughly 0.014–0.016%. These closely agreeing forces are consistent with mostly different energy-zero conventions and small coefficient differences; exact coefficient attribution is not yet established. |
| Nonbonded interactions | Original OpenMM energy lacks the reference cutoff shift. A test-only +epsilon correction reduces the energy difference from hundreds of kJ/mol to 0.008–0.010 kJ/mol. Force RMS differences are 0.00077–0.00088 kJ/mol/nm, roughly 0.0012% of the reference nonbonded-force RMS. Remaining small differences have not been assigned uniquely to rounding/parameters. |

The complete-system force RMS differs by approximately 9–10 parts per million relative to the reference RMS for these three conformations. That total is dominated by bonds; the isolated force comparisons above are more informative than the total alone.

## Interpretation and corrections to earlier discussion

The reference Lennard–Jones energy is shifted to zero at its pair cutoff. Our XML reports an unshifted, truncated energy. The resulting total-energy difference is conformation-dependent because the number of pairs inside the cutoff changes. Therefore raw total energies are not directly interchangeable between the two implementations.

However, adding the shift does not change the forces returned by these implementations inside the cutoff or outside it; the test confirms unchanged force differences. OpenMM’s ordinary derivative of the step/select expression does not supply an impulse at a crossing. Consequently, this energy-reporting difference is not evidence that it generated the expanded or low-helicity ensembles. It matters for interpreting or using recorded energies in Boltzmann weighting or energy comparisons. All conclusions should distinguish the reported discontinuous energy from the smooth piecewise forces that actually drove MD.

A near-constant torsion energy offset also does not itself change conformational forces or normalized ensemble weights. The small residual force discrepancies merit exact coefficient auditing if strict engine equivalence is required, but these tests do not reveal a large force error.

The large bond mismatch in the initial constrained comparison is not a port failure: OpenMM omits constrained bond potential terms while the LAMMPS input does not constrain them. Removing constraints restores bond agreement.

There is no demonstrated catastrophic CMAP, bias-strength, unit-conversion or GPU-precision error in these tests. The poor agreement with desired Tau-5 conformational properties may reflect sparse/unsuitable AF2 information and the random-coil prior, but that remains a hypothesis. This audit does not establish equilibrium populations, adequate prior support, correctness of all protein/topology cases, or experimental agreement. No full LAMMPS Tau-5 production run has been performed.

## Reproduction

Run compare_engines.py for AUDIT_CASE=initial, perturbed, expanded. Set SHIFT_LJ=1 for the test-only nonbonded shift variant. Then run isolate_force_terms.py. Run compare_bias.py; set GPU_COMPARE=1 for OpenCL single precision. Scripts assume the original workspace layout and reference binaries under work/. Raw inputs, LAMMPS logs/dumps, PLUMED logs and force dumps, numerical JSON/CSV results, and force arrays are retained beside this report.

Reference sources: https://github.com/COSBlab/bAIes-IDP ; https://github.com/plumed/plumed2/tree/v2.10 ; https://download.lammps.org/tars/lammps-2Aug2023.tar.gz . The original production inputs and trajectories were not modified during this audit.

# Tau-5* SEC-SAXS: reproduction, window selection and EFA

**Latest update (22 September):** actual EFA component curves and an executed audit are now in [04_EFA_component_Iq.ipynb](04_EFA_component_Iq.ipynb) and [COMPONENT_RESULTS.md](COMPONENT_RESULTS.md). WT has a distinguishable shoulder/main separation; AA purification is not established. The earlier [CHECKPOINT.md](CHECKPOINT.md) is retained as historical context.

Start with the **executed notebooks, in order**:

1. `01_reproduce_and_audit_Rg.ipynb` reproduces Natalie's WT **24.904 Å** and AA **30.785 Å** Guinier fits from her saved profiles, then shows q-range, error-propagation, baseline and frame-window changes separately.
2. `02_windows_profiles_and_EFA.ipynb` recomputes candidate averaging windows, Guinier and Bayesian P(r) fits, and clean-region EFA/noise checks. It includes the poster-style I(q), P(r), Guinier/residual and SEC plots.
3. `03_report_EFA_sensitivity.ipynb` reconstructs the settings in the supplied RAW reports and checks the proposed four-WT/three-AA component rotations against q-range changes. This is not an exact recreation of the unavailable RAW sessions.

All code cells have been run in fresh Python kernels. The original source notebooks remain unchanged in `BNL_SAXS/`; they are references, not the executed notebooks created here. Original RAW GUI/ATSAS dependencies are not required for the new notebooks.

## Main interpretation

The apparent AA Guinier disagreement largely comes from changing the fit range: on exactly the same saved average, 30.785 Å becomes 32.428 Å with the conservative low-q convention. Shared background error propagation and the selected-window change have much smaller effects. This comparison is a method audit, not grounds to dismiss Natalie's result.

Current candidate presentation windows are **WT 159–175** and **AA 169–188** (zero based, inclusive). Direct-average Guinier Rg values are **25.8 ± 0.6 Å** and **32.3 ± 0.5 Å**, conditional on the analysis choices. Nearby windows are nearly tied; these are not unique optimal or proven pure-species windows.

P(r) Rg estimates must be compared separately: current BIFT **27.2 / 34.6 Å** versus supplied poster GNOM **27.3 / 34.23 Å**, WT/AA. The AA Dmax difference (**~159 vs 139 Å**) remains unresolved and model dependent. GNOM has not been rerun here.

The clean traces have one dominant scattering pattern. Weak extra modes depend on the error model. The earlier WT shoulder adds a robust pattern, but the recovered profile is not sufficient to assign an oligomer. Constrained multi-component report reconstructions produce noisy/negative component profiles and q-range sensitivity. **No oligomer identity or abundance is established.**

`EFA_equations.md` walks through the equations. `EFA_report.md` contains detailed results and caveats. `poster_results/README.md` documents figure settings. `results/Rg_correction_audit.json` records the numerical change log.

## Reproduce locally

Use Python 3.12. Run commands from this repository root:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m ipykernel install --prefix .jupyter --name saxs --display-name 'SAXS Python'
.venv/bin/python analysis/execute_notebooks.py
```

`execute_notebooks.py` locates the project-local kernel automatically. To browse interactively, open the notebooks in Jupyter or VS Code and select this environment. `analysis/build_*notebook*.py` regenerates notebook sources and clears outputs; use only when intentionally rebuilding, then execute again.

The package contains all source-repository files, the supplied context documents, results, code and the minimal unchanged BioXTAS RAW modules needed for BIFT. Environments, caches and duplicate Git history are omitted. The full upstream RAW checkout remains in the Desktop working folder. The original two requested HDF5 series are already facility-subtracted; no second subtraction of the original facility buffer is performed. Residual background models are explicit. q is interpreted in Å⁻¹ from the source notebooks; the HDF5 unit field is blank.

## Provenance and licenses

Input repository: https://github.com/Natalie-Loui/BNL_SAXS at `4a6ac1c514d25a589e671a91efd1f38c4bcb316c`.

Numerical EFA/rotation and BIFT code: https://github.com/jbhopkins/bioxtasraw at `90f86f96cdc907d0bd41d3ba96d4936a9b1b0e2a`, under its included GPL license. Original data and context materials retain their authors' rights; no new license is assigned to them. File hashes and package versions are in `results/provenance.json`, `results/reproduction_input_hashes.json`, and `MANIFEST.sha256`.

Private GitHub repository: https://github.com/paulrobustelli/tau5-sec-saxs-analysis.

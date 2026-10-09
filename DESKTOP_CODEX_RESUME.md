# Desktop Codex: resume the Tau-5 bAIes/OpenMM work

Read this file first. Do not reconstruct the workflow from chat history.

## Authoritative files

Download these two archives from the repository root:

1. `Tau5_bAIes_10ns_results_2026-10-09.zip` — completed WT and AA 10 ns
   trajectories, final portable XML states, checkpoints, logs, analysis, and
   numerical distogram audit.
2. `Tau5_bAIes_OpenMM_restart_2026-10-09_v2.zip` — patched OpenMM runner,
   force field, CMAP, and restart scripts.

Also read:

- `BAIES_RESULTS_10NS.md`
- `BAIES_DISTOGRAM_AUDIT.csv`
- `BAIES_HANDOFF.md`

The completed 10 ns `final_state.xml` files in the results archive supersede
the 9.30–9.34 ns states in the older restart archive. Resume from the completed
states unless explicitly asked to reproduce the earlier snapshots.

## What has already been established

- Both WT and W397A/W433A (`AA`) validation runs completed 5,000,000 steps at
  2 fs per step: 10 ns and 1,000 saved frames each.
- The constructs contain GP followed by AR330–447, for 120 residues total.
- WT construct positions 70 and 106 are W/W; AA positions 70 and 106 are A/A.
  These correspond to AR397 and AR433.
- The AlphaFold distograms are `120 × 120 × 64`, sequence-aligned, exactly
  symmetric, and mapped to the correct CB atoms (CA for glycine).
- Reprocessing regenerates every stored restraint exactly: 19 for WT and 15
  for AA. Means and widths agree within text-rounding error.
- The 10 ns trajectories are not converged production ensembles. They have
  essentially zero DSSP alpha helicity and should be treated as implementation
  validation only.

## Create the environment

Use Python 3.11 or 3.12:

```bash
python3 -m venv baies_env
source baies_env/bin/activate
python -m pip install "openmm==8.6.1" numpy scipy pandas matplotlib mdtraj
```

Unzip both packages into one working directory. The following names assume:

```text
Tau5_bAIes_OpenMM_restart_2026-10-09/
Tau5_bAIes_10ns_results_2026-10-09/
```

## Mandatory short restart test

Run only 10,000 additional steps first. On an Apple Silicon Mac, try OpenCL;
fall back to CPU if OpenCL is unavailable. On an NVIDIA machine, use CUDA.

WT example:

```bash
python Tau5_bAIes_OpenMM_restart_2026-10-09/code/bAIes_openMM.py \
  -ff Tau5_bAIes_OpenMM_restart_2026-10-09/code/coil_202602.xml \
  -pdb Tau5_bAIes_10ns_results_2026-10-09/Raw/WT/protein.pdb \
  -xtc WT_restart_test.xtc \
  -cmap Tau5_bAIes_OpenMM_restart_2026-10-09/code/cmap_20240524.cmap \
  -baies Tau5_bAIes_10ns_results_2026-10-09/Raw/WT/baies_params.dat \
  -restart_state Tau5_bAIes_10ns_results_2026-10-09/Raw/WT/final_state.xml \
  -Nsteps 10000 -save_xtc 5000 -save_stdout 5000 \
  -checkpoint WT_restart_test.chk -state_out WT_restart_test_state.xml \
  -platform OpenCL
```

Repeat for AA by replacing every WT path/name with AA. Confirm that each test:

1. loads at step 5,000,000;
2. completes without NaNs or exceptions;
3. reports a temperature near 298 K;
4. writes two trajectory frames, a checkpoint, and a portable XML state.

Do not begin a long run until both short tests pass.

## Production continuation

The published bAIes sampling scale was 20,000 frames separated by 100 ps,
equivalent to 2 microseconds. At a 2 fs time step, 100 ps corresponds to 50,000
steps. To extend the existing 10 ns trajectory to 2 microseconds total requires
995,000,000 additional steps.

For a production continuation, use:

```text
-Nsteps 995000000
-save_xtc 50000
-save_stdout 50000
```

Choose `-platform CUDA` on NVIDIA hardware, `-platform OpenCL` on a compatible
Mac, or `-platform CPU` as fallback. Keep WT and AA in separate directories.
Write checkpoints and portable XML states regularly. Do not overwrite the
completed 10 ns files.

Before launching 995 million steps, estimate wall time from a 1–10 ns benchmark
on the chosen hardware and report the estimate to the user.

## Analysis requirements

Preserve the completed 10 ns bundle unchanged. Analyze continuations with the
matching `protein.pdb` and concatenate segments in chronological order. Report:

- Rg distributions and convergence by time block;
- DSSP alpha helicity and per-residue helicity;
- R2 (AR395–405) and R3 (AR430–440) helicity;
- restrained-pair distance distributions versus each Gaussian mean/width;
- Sa, contact maps, and comparisons with the ASTEROIDS and 2022 MD ensembles.

The existing audit shows that distogram indexing is correct. It also shows that
the 10 ns mean restrained distances remain roughly 0.25–0.27 nm above their
Gaussian means. Track whether that discrepancy decreases during longer runs.

## Prompt to give Desktop Codex

> Clone or open `paulrobustelli/tau5-sec-saxs-analysis`, read
> `DESKTOP_CODEX_RESUME.md` and `BAIES_RESULTS_10NS.md`, download both bAIes ZIP
> archives, and resume from the completed WT and AA 10 ns `final_state.xml`
> files. Run the mandatory 10,000-step restart tests first and report the exact
> results before launching production sampling.


# Tau-5 bAIes/OpenMM handoff

Last updated: 2026-10-09 (America/New_York)

## Objective

Continue and analyze bAIes/OpenMM simulations for the 120-residue GP–Tau-5* WT
and W397A/W433A (`AA`) constructs.

## Download

- Restart package: `Tau5_bAIes_OpenMM_restart_2026-10-09_v2.zip`
- Checksum: `Tau5_bAIes_OpenMM_restart_2026-10-09_v2.zip.sha256`
- Chat snapshot: https://chatgpt.com/s/cx_6ac91f2da44081919c31013c1f164902

The package contains portable OpenMM XML states, binary checkpoints, partial
XTC trajectories, topology PDBs, bAIes restraint parameters, the patched
OpenMM runner, force-field/CMAP files, restart scripts, and early Rg/helicity
analysis. The compact package omits the large AlphaFold distogram pickles;
their derived restraints are already present in `baies_params.dat` and the
pickles are not required to restart.

## Uploaded restart points

| Construct | Step | Time | Remaining to 10 ns |
|---|---:|---:|---:|
| WT | 4,670,000 | 9.34 ns | 330,000 steps |
| AA | 4,650,000 | 9.30 ns | 350,000 steps |

Use `restart_state.xml` when changing computers. The binary `.chk` files are
less portable and require a compatible OpenMM version, platform, and hardware.

## Environment and restart

```bash
python3 -m venv baies_env
source baies_env/bin/activate
python -m pip install "openmm==8.6.1" numpy
unzip Tau5_bAIes_OpenMM_restart_2026-10-09_v2.zip
cd Tau5_bAIes_OpenMM_restart_2026-10-09
bash restart_WT.sh /absolute/path/to/baies_env/bin/python
bash restart_AA.sh /absolute/path/to/baies_env/bin/python
```

The scripts load the portable XML states and write continuation trajectories,
logs, checkpoints, and final portable states without replacing the partial
trajectories.

## Current local status

After the snapshots were packaged, both original local validation simulations
completed 5,000,000 steps (10 ns) successfully. Those completed local outputs
were not included in the v2 restart archive; the uploaded portable states are
the restart points listed above.

## Suggested next work

1. Verify a short continuation from each XML state on the new machine.
2. Continue on a GPU for production sampling; the bAIes paper-scale protocol is
   20,000 frames at 100 ps intervals (2 microseconds total).
3. Concatenate or analyze partial and continuation XTC files in time order with
   the matching `protein.pdb` topology.
4. Recompute Rg, DSSP helicity, Sa, and contact maps after sufficient sampling.


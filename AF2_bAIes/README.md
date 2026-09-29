# AF2 inputs for bAIes

Repaired, sequence-specific notebooks for GP + AR330–447 (120 residues). Use Colab runtime **2025.07**, T4 GPU. The original distogram fork is https://github.com/zshengyu14/ColabFold_distmats . Fixes: legacy TPU-import Python scoping and overlapping TensorFlow distributions (single tensorflow-cpu 2.18.0). Save-all retains full prediction dictionaries, including distogram logits and bin edges.

- [WT notebook](Tau5star_WT_AF2_distograms.ipynb)
- [AA notebook](Tau5star_AA_AF2_distograms.ipynb)

WT completed five models. Each distogram has shape 120 × 120 × 64; all five PDBs and prediction dictionaries are retained in the local results ZIP. AA is being run separately. These are AF2 outputs, not bAIes simulation ensembles. The patched local PLUMED/LAMMPS environment is being prepared.

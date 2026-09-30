# Tau-5* master ensemble analysis

[Open the executed notebook](Tau5_ASTEROIDS_master.ipynb) · [Download notebook](Tau5_ASTEROIDS_master.ipynb?raw=true) · [Standalone HTML](Tau5_ASTEROIDS_master.html?raw=true) · [Download portable bundle](Tau5_ASTEROIDS_master_bundle.zip?raw=true)

All four historical starting pools and 12 equal-weight selections, plus the preferred WT and AA aggressive expanded-pool fits under two low-q schedules. Includes full-chain and matched AR391–446 structural comparisons with57,144 apo WT MD frames, SAXS before/after, low-q diagnostics, secondary shifts, signed/absolute residuals, DSSP-H populations, Sα–Rg, contact maps and helix-length/Sα histograms.

[Numerical ensemble summary](ensemble_summary.csv) · [SAXS error by q region](saxs_q_region_summary.csv)

Source result folders: [WT aggressive trial](WT_Expanded_Aggressive_Trial) · [AA aggressive trial](AA_Expanded_Aggressive_Trial)

The notebook is executed with embedded outputs. To rerun, extract the bundle and open the notebook beside analysis.py and data/. Dependencies: numpy,pandas,matplotlib,seaborn. Cached data permit regeneration of plots without trajectories. master_compute.py records the workspace-relative raw-data regeneration procedure. Units, definitions, limitations and sources are in the notebook. These are subset selections, not continuous reweighting; three optimizer seeds are not independent candidate pools.

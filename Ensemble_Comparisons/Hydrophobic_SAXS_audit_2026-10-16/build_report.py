from pathlib import Path
import json,hashlib,platform
import pandas as pd,numpy as np,mdtraj
O=Path(__file__).resolve().parent
summary=pd.read_csv(O/'ensemble_packing_before_after.csv');fits=pd.read_csv(O/'diagnostic_fit_results.csv');matched=pd.read_csv(O/'Rg_matched_WT26_30_equal_1A_bins.csv')
text='''Hydrophobic packing and the SAXS fitting advantage of REST

The calculations support a deficit of compact, tertiary-packed conformations
in the other pools, but do NOT establish that Rg<19 A states are necessary.
All results are conditional on these pools and the common SAXS forward model.

1. Sampling: REST has976/1866 frames below19 A; IDPForge WT/AA and bAIes WT/AA
have0/1000 each. ASTEROIDS priors have2/2000 each. The existing combined fit
reduces REST weight below19 A from52.3% to11.5%.

2. Rg alone does not describe the difference. Equal weighting in four1-A bins
from26 to30 A gives mean Rg27.93/28.02/28.09/28.02 A for REST/IDPForge/bAIes/
ASTEROIDS WT, but nonlocal buried HSA1134/359/51/300 A^2 and hydrophobic contact
pairs18.4/4.7/0.37/3.2. Thus REST retains substantially more nonlocal packing
even at similar overall size. These describe geometry, not verified native cores.

3. Within each non-REST WT and AA pool, high-burial groups give better
intermediate-q agreement than low-burial groups with matched2-A Rg-bin counts.
For WT: IDPForge mid chi2=189.7 vs110.8, bAIes69.8 vs49.6, ASTEROIDS52.5 vs18.7
(low versus high burial, uniform group averages with a full-q scale). Differences
remain large: choosing somewhat more buried existing frames is not equivalent
to supplying REST-like tertiary packing.

4. Direct pool-addition experiments: adding160 compact REST WT donors improves
mid-q chi2 from125.75 to0.82 for IDPForge,178.04 to1.47 for bAIes, and34.95 to0.67
for ASTEROIDS. Equal-size additions of REST Rg>=25 A donors yield1.80/3.75/3.53.
Compact additions get approximately55/57/42% fitted weight, respectively.
Two additional compact-donor selections reproduce improvement: IDPForge
0.77–0.85, bAIes1.28–1.52, ASTEROIDS0.64–0.67. These ranges reflect donor choice
only; not confidence intervals. Synthetic mixed pools are not equilibrium models.

5. Necessity tests qualify the compact-state interpretation. REST mid-q chi2
is0.58 with all frames,0.59 after removing every Rg<19 A frame,0.63 after removing
Rg<20 A,0.93 after removing Rg<22 A, and5.06 after removing Rg<24 A. Absolute Kish
support stays at least186.6 in every case. Loss of fitting flexibility/cardinality
also contributes as the remaining pool shrinks. No unique physical cutoff is
established. The broad REST pool, or its less-extremely-compact members, already
contains adequate SAXS curves.

6. More burial/aromatic coordination is not monotonically better within REST:
Rg-matched low/high-burial subset fits give mid chi2=1.39/4.77; low/high-aromatic
subsets give0.93/2.70. The useful feature is a diversity of internal distance
patterns and packing, not simply maximal collapse or maximal contact count.

Interpretation: REST contributes scattering profiles missing or severely
undersampled in the other pools. Distinct nonlocal packing at matched Rg and
successful donor additions support the tertiary-packing hypothesis. The tests
do not uniquely prove a two-population hydrophobic-core/extended equilibrium,
that Rg<19 A states are required, or that other generated structures are fake.
SAXS cannot identify a unique core topology. All donors change multiple structural
features simultaneously, including known geometry/cis features.

New fits are SAXS-only Borthakur Kish diagnostics, distinct from the existing
combined CS/SAXS fits. One profiled global scale; same frozen REST block SEM
across WT tests; original absolute Kish support retained. All24 fits passed
numerical checks and independent curve reconstruction. See METHODS.txt for
precise definitions, uncertainty limitations, and scaling conventions.
'''
(O/'FINDINGS.txt').write_text(text)
cells=[]
def md(s):cells.append({'cell_type':'markdown','metadata':{},'source':s})
def code(s):cells.append({'cell_type':'code','metadata':{},'source':s,'execution_count':None,'outputs':[]})
md('# Hydrophobic packing and SAXS: hypothesis tests\n\nAll completed common121 ensembles. WT and AA remain separate. New harmonic bAIes production is not included. New SAXS-only diagnostic fits are distinct from published combined CS/SAXS weights.')
code("from pathlib import Path\nimport pandas as pd,json\nfrom IPython.display import display,Image,Markdown\nroot=Path.cwd()\nassert (root/'FINDINGS.txt').exists()\nprint((root/'FINDINGS.txt').read_text())")
md('## Packing at comparable dimensions\nPoints are mean values within2-A Rg bins with at least5 frames; the table standardizes WT26–30 Å to equal mass in each1-A bin. Counts and full per-frame data are retained.')
code("display(pd.read_csv(root/'Rg_matched_WT26_30_equal_1A_bins.csv'))\ndisplay(Image(filename=str(root/'packing_at_matched_Rg.png')))")
code("display(Image(filename=str(root/'coordination_vs_buried_HSA.png')))\ndisplay(pd.read_csv(root/'ensemble_packing_before_after.csv'))\ndisplay(pd.read_csv(root/'tertiary_contacts_before_after.csv'))")
md('## Do added compact structures rescue the curves?\nSAXS-only fits; same error model and original absolute Kish support. Gray: original pool. Blue:160 compact REST donors. Orange:160 more expanded REST donors. Compact-donor replicates are diagnostic sensitivity, not independent validation.')
code("display(Image(filename=str(root/'WT_pool_addition_experiments.png')))\ndf=pd.read_csv(root/'diagnostic_fit_results.csv')\ndisplay(df[['label','N','Neff','chi2_all','chi2_mid','added_REST_weight','Rg_CA_A','buried_HSA_nonlocal_A2']])")
md('## Are very compact or highly buried conformations necessary?\nAll fits retain Kish Neff≥186.6. The below19-Å population is not required for the REST SAXS fit. Stronger removals reduce both structural coverage and fitting flexibility.')
code("display(Image(filename=str(root/'REST_removal_experiments.png')))")
md('## SAXS contributions in the REST diagnostic fit\nEvery component shares the SAME fitted scale; components sum to the total. Contributions do not establish necessity; that is tested separately by removals.')
code("display(Image(filename=str(root/'REST_fitted_SAXS_components.png')))")
md('## Conditional SAXS by conformation class\nEach group has uniform weights and ONE scale fitted over the full q range. Rg boundaries19/25/35 Å; other boundaries are within-pool tertiles. See table for numeric cutpoints. These panels compare shape after optimizing group amplitude, not different per-conformer scales.')
code("display(Image(filename=str(root/'conditional_SAXS.png')))\ndisplay(pd.read_csv(root/'conditional_group_scores.csv'))")
md('## Alternative view with no group rescaling\nThe same conditional raw curves now share the fixed scale of their full unweighted parent pool. This makes the normalization distinction explicit.')
code("display(Image(filename=str(root/'conditional_SAXS_fixed_parent_scale.png')))")
md('## Conditional SAXS by general tertiary-contact count\nMinimum heavy-atom distance<5 Å, sequence separation>6. This includes polar contacts and is separate from hydrophobic coordination.')
code("display(Image(filename=str(root/'tertiary_contact_SAXS.png')))\ndisplay(pd.read_csv(root/'matched_Rg_burial_group_scores.csv'))")
md('## Methods, validation and reproducibility\nSASA uses [MDTraj Shrake–Rupley](https://mdtraj.readthedocs.io/en/latest/api/generated/mdtraj.shrake_rupley.html). Recompute structures with compute.py and supplement.py; analyze.py generates groups; fit_experiments.py and additional_experiments.py guard completed fits against reruns; finalize.py and decompose.py generate comparison figures. Original pools and published weights are unchanged.')
code("print((root/'METHODS.txt').read_text())\ndisplay(pd.DataFrame(json.loads((root/'SASA_resolution_check.json').read_text()))[['ensemble','mean_abs_difference_A2','max_abs_difference_A2','spearman']])\ndisplay(pd.DataFrame(json.loads((root/'independent_fit_validation.json').read_text())))")
nb={'nbformat':4,'nbformat_minor':5,'metadata':{'kernelspec':{'name':'python3','display_name':'Python 3','language':'python'}},'cells':cells}
for i,c in enumerate(cells):c['id']=f'hsa-audit-{i}'
(O/'Hydrophobic_SAXS_hypothesis.ipynb').write_text(json.dumps(nb,indent=1))
files={}
for p in O.rglob('*'):
 if p.is_file() and p.suffix in ['.csv','.json','.py','.txt','.npz'] and 'chunk_' not in p.name and '__pycache__' not in str(p) and p.name!='provenance_manifest.json':files[str(p.relative_to(O))]=hashlib.sha256(p.read_bytes()).hexdigest()
(O/'provenance_manifest.json').write_text(json.dumps({'software':{'python':platform.python_version(),'numpy':np.__version__,'mdtraj':mdtraj.__version__},'files':files},indent=2))
print('notebook and report created',len(cells),'cells',len(fits),'fits')

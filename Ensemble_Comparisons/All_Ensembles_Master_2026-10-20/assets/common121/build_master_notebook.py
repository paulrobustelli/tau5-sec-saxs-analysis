"""Build a portable master notebook from validated comparison outputs, never rerun predictors."""
from pathlib import Path
import json,sys,hashlib
O=Path(__file__).resolve().parent
names=['REST_WT','bAIes_WT','bAIes_AA','IDPForge_WT','IDPForge_AA'];pending=[n for n in names if not (O/n/'fit_complete.json').exists()]
if pending and '--draft' not in sys.argv:raise SystemExit('Awaiting validated fits: '+', '.join(pending))
cells=[]
def md(s):cells.append({'cell_type':'markdown','metadata':{},'source':s.splitlines(True)})
def code(s):cells.append({'cell_type':'code','execution_count':None,'metadata':{},'outputs':[],'source':s.splitlines(True)})
md('# Tau5* ensemble comparison — 121 residues\n'+('**DRAFT: pending '+', '.join(pending)+'**\n' if pending else '')+'GP + AR330–448, including terminal Cys448. R2 main: AR390–404; R2 minor: AR404–410 (shared residue404); R3: AR430–440.')
md('## Methods and interpretation\nAll new SAXS predictions use Pepsi3.0, dro5, default Fraser radius1.61 Å, q0–0.6 Å⁻¹, 601 samples, n80, no smearing, and one fitted ensemble intensity scale with no offset. No individual-conformer scales. Chemical shifts use SPARTA+; explicitly authorized reuse of original120-residue predictions for the terminally extended bAIes/IDPForge pools. Common chemical-shift observations cover AR330–446; missing observations are excluded and counts are reported.\n\nBorthakur per-datatype calibration followed by global Kish regularization retains normalized Kish≥0.10. Sampling error uses fixed uniform-prior IID estimates for IDPForge and100ns independent blocks for bAIes/REST, with block-size scans. REST saved frame order is not physical kinetics. Numerical gradient/primal-dual checks validate the optimization; these fitted experimental data are not independent validation.')
md('## Construct and sampling provenance\nIDPForge and bAIes dynamics originally contained120 residues. Cys448 was appended afterward with validated peptide geometry and clash checks; the terminal Pro carbonyl O was rotated when needed, and Pro OXT removed. This produces121-residue forward-model inputs, not newly equilibrated121-residue simulations. Original pools and results remain preserved. bAIes uses the exploratory8Å/seqsep2 selection sensitivity variant; it is not an exact paper reproduction or a validated upstream correction.\n\nASTEROIDS references preserve the existing500-member selected ensembles and2000-member priors. SAXS is reevaluated from archived canonicaldro5 intensities on the shared experimental grid; historical selection optimized hydration separately. These reference curves are not a new selection or Kish fit.')
code("from pathlib import Path\nimport json\nimport pandas as pd\nfrom IPython.display import display, Image\nroot=Path.cwd()\nassert (root/'structural_comparison.csv').exists(), 'Run with this notebook directory as cwd.'")
md('## Structural and experimental summary')
code("display(pd.read_csv(root/'master_comparison_table.csv'))\ndisplay(pd.read_csv(root/'experimental_agreement.csv'))\ndisplay(pd.read_csv(root/'chemical_shift_comparison_coverage.csv'))")
md('## Rg distributions\nRg is geometric Cα Rg, not the SAXS forward-model Rg.')
code("display(Image(filename=str(root/'master_figures/Rg_before_after.png')))")
md('## Residue helicity profiles')
code("display(Image(filename=str(root/'master_figures/helicity_before_after.png')))")
for name in ['ASTEROIDS_WT','ASTEROIDS_AA']+names:
 md('## '+name.replace('_',' '))
 code("for suffix in ['SAXS_before_after','CS_before_after','SaRg_contacts']:\n    path=root/'master_figures'/('"+name+"_'+suffix+'.png')\n    if path.exists(): display(Image(filename=str(path)))")
 if name in names:code("path=root/'"+name+"/fit_validation.json'\nif path.exists(): display(json.loads(path.read_text()))\npath=root/'"+name+"/weight_support_audit.json'\nif path.exists(): display(json.loads(path.read_text()))")
md('## Geometry and support limitations\nThe supplied REST pool is compact; at normalized Kish≥0.10 its maximum possible mean Cα Rg is27.556 Å. This bound does not directly constrain a SAXS-derived Rg. REST reweighting can fit the supplied SAXS curve despite a lower Cα Rg; residual structure must still be inspected.\n\nThe polyproline audit found substantial incoming cis occupancy at REST Pro377 and a long trapped cis state at bAIes AA Pro372. IDPForge also contains peptide ω values outside cis/trans wells and positive Pro φ values. These are retained and reported, not silently removed. No cis/trans switching rates are inferred from REST. The student’s planned all-trans continuation is a separate future ensemble.')
code("path=root/'REST_WT/PPII_cis_before_after.csv'\nif path.exists(): display(pd.read_csv(path))")
notebook={'cells':cells,'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.10'}},'nbformat':4,'nbformat_minor':5}
for i,c in enumerate(cells):c['id']='tau5-'+str(i)
name='Tau5_121_master_DRAFT.ipynb' if pending else 'Tau5_121_master.ipynb';(O/name).write_text(json.dumps(notebook,indent=1));print(name)

from pathlib import Path
import shutil,json,datetime,hashlib
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;A=O/'assets';A.mkdir(exist_ok=True)
# Copy curated common-protocol publication containing all matrices, weights and audits.
src=R/'work/tau5-sec-saxs-analysis/Ensemble_Comparisons/Common121_2026-10-16'
shutil.copytree(src,A/'common121',dirs_exist_ok=True,ignore=shutil.ignore_patterns('*.ipynb','*.html','SHA256SUMS*'))
H=R/'work/asteroids_reference'
shutil.copytree(H/'data',A/'asteroids_history/data',dirs_exist_ok=True)
shutil.copy2(H/'analysis.py',A/'asteroids_history/analysis.py')
for f in ['ensemble_summary.csv','saxs_q_region_summary.csv']:shutil.copy2(H/f,A/'asteroids_history'/f)
groups={'idp200':'IDPForge_MaxEnt_2026-10-13','idp1000_superseded':'IDPForge_D2D_MaxEnt_2026-10-15','baies_original':'Tau5_bAIes_OpenMM_restart_2026-10-09/Analysis','baies_full':'bAIes_8A_seqsep2_2us_2026-10-14/final_region_analysis','diagnostics':'AF2_distogram_recompute_2026-10-14','geometry':'PPII_dihedrals_2026-10-15'}
for dest,source in groups.items():
 d=A/dest;d.mkdir(exist_ok=True)
 for p in (R/'outputs'/source).iterdir():
  if p.is_file() and p.suffix in ['.csv','.png','.json'] and p.stat().st_size<8_000_000:shutil.copy2(p,d/p.name)
D=A/'harmonic';D.mkdir(exist_ok=True)
for f in ['equilibration_validation.json','harmonic_validation.json']:shutil.copy2(R/'outputs/bAIes_harmonic_2026-10-16'/f,D/f)
status={}
for v in ['WT','AA']:
 p=R/'outputs/bAIes_harmonic_2026-10-16'/v;lines=(p/'production.log').read_text().splitlines();numeric=[x for x in lines if x and x[0].isdigit() and '\t' in x];last=numeric[-1] if numeric else None
 status[v]={'status':'completed' if 'Completed 1000000000 steps' in lines else 'in progress (log snapshot)','last_log_record':last,'production_ns':(int(last.split('\t')[0])-10000000)*.000002 if last else 0,'reweighted':False}
(D/'status.json').write_text(json.dumps(status,indent=2))
cells=[]
def md(s):cells.append(dict(cell_type='markdown',metadata={},source=s))
def code(s):cells.append(dict(cell_type='code',metadata={},source=s,outputs=[],execution_count=None))
md('# Tau5 — all-ensemble comparison master\nSnapshot: '+datetime.datetime.now(datetime.timezone.utc).isoformat()+'\n\nA single review of structural properties and experimental agreement before/after fitting. Current common121 results lead; historical protocols are explicitly separated. Missing fits are unavailable, never zero. New harmonic production is ongoing.\n\n## Contents\n1. Current common121 comparison: ASTEROIDS WT/AA, REST WT, IDPForge WT/AA, bAIes WT/AA\n2. Full bAIes production and restraint diagnostics\n3. Harmonic equilibrations and ongoing production\n4. Original bAIes and historical IDPForge fits\n5. All historical ASTEROIDS selections, aggressive fits, and 2022 apo WT R2R3 MD\n6. Provenance and limitations')
code("from pathlib import Path\nimport os,json, pandas as pd, numpy as np\nimport matplotlib.pyplot as plt\nfrom IPython.display import display, Image\nROOT=Path.cwd(); assets=ROOT/'assets'\nassert assets.exists()\ndef table(folder,name):\n    p=assets/folder/name\n    if p.exists(): display(pd.read_csv(p))\ndef figures(folder,names):\n    for n in names:\n        p=assets/folder/n\n        if p.exists(): display(Image(filename=str(p)))\n        else: print('Unavailable:',p.name)")
# Reuse the complete scientifically checked common121 notebook, with a portable root.
n=json.loads((R/'outputs/Ensemble_Comparison_121_2026-10-15/Tau5_121_master.ipynb').read_text())
for c in n['cells']:
 s=''.join(c['source']);s=s.replace('root=Path.cwd()',"root=assets/'common121'")
 (md if c['cell_type']=='markdown' else code)(s)
md('## Compact before/after comparison\nPaired raw metrics are shown per ensemble. ASTEROIDS after = historical subset selection reevaluated under the common SAXS protocol; other after = combined Kish MaxEnt. Chemical-shift nuclei and coverage differ by sequence; compare each paired ensemble first. These are fitted observations, not held-out validation.')
code("df=pd.read_csv(assets/'common121/master_comparison_table.csv')\nmetrics=['Rg_CA_A','R2_main_H_percent','R2_minor_H_percent','R3_H_percent','Sa_full','SAXS_chi2_mean']\nfig,axs=plt.subplots(2,3,figsize=(16,9),layout='constrained')\nfor ax,m in zip(axs.flat,metrics):\n    for i,(name,g) in enumerate(df.groupby('ensemble',sort=False)):\n        g=g.sort_values('stage',key=lambda x:x.map({'prior':0,'selected':1,'MaxEnt':1,'combined':1}).fillna(1))\n        ax.plot([i]*len(g),g[m],color='gray')\n        ax.scatter(i,g.iloc[0][m],marker='o',color='gray',label='Before' if i==0 else None)\n        ax.scatter(i,g.iloc[-1][m],marker='D',color='tab:blue',label='After' if i==0 else None)\n    ax.set_xticks(range(df.ensemble.nunique()),df.ensemble.unique(),rotation=65,ha='right');ax.set_title(m);ax.legend()\nplt.show()")
md('## Full Jeffreys bAIes production: 20,000 frames each\nOriginal 120-residue trajectories, exploratory 8 Å / seqsep2 selection including i,i+3. These full-trajectory summaries differ from the 1,000-frame frozen-prefix Cys-extended pools used for fitting. μ ± σ occupancy describes soft-restraint agreement; it is not a hard pass/fail test.')
code("table('baies_full','convergence_summary.csv')\nfigures('baies_full',['convergence.png','production_restraint_satisfaction.png'])\nfor v in ['WT','AA']: table('baies_full',v+'_restraint_agreement.csv')")
md('## Harmonic bAIes-N likelihood: equilibrations complete; production ongoing\nFixed fitted σ as SD; U = ½ kBT ((r−μ)/σ)². Same 8 Å / seqsep2 selection. 20 ns WT/AA checks are complete; 2 μs additional production is ongoing. Short equilibration helicities are not equilibrium estimates. No CS/SAXS reweighting exists for these new runs yet.')
code("display(json.loads((assets/'harmonic/status.json').read_text()))\ndisplay(pd.DataFrame(json.loads((assets/'harmonic/equilibration_validation.json').read_text())).T)")
md('## Original bAIes production and short selection diagnostics\nHistorical 120-residue original-restraint ensembles and matched short diagnostics. No before/after CS/SAXS fit is available for the original full trajectories. Historical region definitions in these tables are retained verbatim and must not be confused with current main/minor definitions.')
code("table('baies_original','summary.csv')\nfigures('baies_original',['observables.png','helicity.png','contacts.png','restraints.png'])\ntable('diagnostics','sensitivity_summary.csv')\ntable('diagnostics','helicity_corrected_R2_regions.csv')")
for folder,title in [('idp200','Historical 200-structure IDPForge fits'),('idp1000_superseded','Historical 1,000-structure delta2D IDPForge fits — superseded protocol')]:
 md('## '+title+'\nThese archived fits used different SAXS settings/scaling and 120-residue inputs. Preserve for provenance; do not rank their χ² against common121 results. The 1,000-structure common121 fits above supersede the earlier weights. Historical R2/R3 definitions remain as saved in these tables.')
 code(f"table('{folder}','ensemble_comparison.csv')\ntable('{folder}','experimental_fit_comparison.csv')\ntable('{folder}','structural_ensemble_comparison.csv')\nfigures('{folder}',['WT_fit.png','AA_fit.png'])")
md('## Polyproline geometry audit\nCis fractions and peptide geometry flags are retained. REST frame-order changes are not physical transition rates. These geometry observables were not explicit fit restraints.')
code("table('geometry','ensemble_summary.csv')")
md('## Historical ASTEROIDS master — all pools and selections\nThe following is the earlier master regenerated from its cached data: four base/expanded pools, three selections per pool, aggressive fits, and the 2022 apo WT R2R3 MD comparison. Historical 120-residue structural exports, hydration/normalization choices, and plotting definitions are retained and labeled. Statements such as “preferred” below refer to that historical snapshot, not the current common121 analysis. 2022 short-fragment apo WT MD is distinct from the supplied full-length REST WT.')
code("import types\nhistory=types.ModuleType('tau5_history')\ndef historical(source):\n    old=Path.cwd()\n    try:\n        os.chdir(assets/'asteroids_history')\n        exec(source,history.__dict__)\n    finally: os.chdir(old)\nhistorical((assets/'asteroids_history/analysis.py').read_text())\ntable('asteroids_history','ensemble_summary.csv')\ntable('asteroids_history','saxs_q_region_summary.csv')")
h=json.loads((H/'Tau5_ASTEROIDS_master.ipynb').read_text())
for c in h['cells']:
 s=''.join(c['source'])
 if c['cell_type']=='code':
  if 'from analysis import' in s:continue
  code('historical('+repr(s)+')')
 else:md(s)
md('## Source manifest and reproducibility\nAll copied input files are SHA256-recorded in source_manifest.json. This notebook reads validated caches and rerenders plots; it does not rerun molecular dynamics, predictors or fits. Full matrices, weights, regularization scans, validation and structural numerical tables are under assets/common121. Constructed ideal helices are not sampled ensembles and are excluded. Historical fits are included for completeness without asserting cross-protocol comparability.')
for i,c in enumerate(cells):c['id']=f'all-{i}'
nb=dict(cells=cells,nbformat=4,nbformat_minor=5,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.10'}})
(O/'Tau5_All_Ensembles_Master.ipynb').write_text(json.dumps(nb,indent=1))
manifest=[{'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(A.rglob('*')) if p.is_file()]
(O/'source_manifest.json').write_text(json.dumps(manifest,indent=2));print(len(cells),'cells;',len(manifest),'source files')

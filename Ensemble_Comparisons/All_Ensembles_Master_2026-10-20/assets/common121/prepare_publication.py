from pathlib import Path
import shutil,json,hashlib
O=Path(__file__).resolve().parent;R=O.parents[1]
D=R/'work/tau5-sec-saxs-analysis/Ensemble_Comparisons/Common121_2026-10-16';D.mkdir(parents=True,exist_ok=True)
for p in O.rglob('*'):
 if not p.is_file() or '__pycache__' in p.parts or p.suffix in ['.log','.pid'] or 'DRAFT' in p.name or p.name.endswith('_in_progress.csv') or 'watch_ready_fits.lock' in p.parts:continue
 if p.name in ['structural_cache.npz','prior_structural_cache.npz']:continue
 q=D/p.relative_to(O);q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
report='''Tau5* common121 ensemble comparison, 16 October 2026

All five new ensembles have completed native SAXS and numerically validated Borthakur per-datatype then global Kish fits (normalized Kish >=0.10). The executed master notebook and HTML contain all seven ensemble comparisons. master_comparison_table.csv is the numerical summary; observable_matrix.npz, weights.csv and regularization_scan.csv retain the fit inputs and outputs.

REST WT SAXS mean chi-square improves 14.9072 to 0.97417, with structured high-q residuals still visible. IDPForge WT and AA remain inconsistent with SAXS after fitting (70.2879 and 262.738). bAIes combined fits improve chemical shifts but worsen SAXS: WT 69.5002 to 110.3662; AA 136.0893 to 260.6174. A joint regularized objective need not minimize the raw SAXS-only chi-square. Fitted observations are not independent validation.

The construct is GP plus AR330-448 including Cys448. IDPForge and bAIes were generated/simulated with120 residues; terminal Cys was appended for SAXS with validated geometry, while their existing chemical shifts were reused by explicit instruction. REST and archived ASTEROIDS are native121-residue inputs. All current SAXS uses Pepsi3.0 dro5, default Fraser radius, one ensemble intensity scale and no offset. Earlier120-residue mixed-protocol fits are superseded and preserved separately.

R2 main is AR390-404; minor AR404-410 shares404; R3 AR430-440. Rg is geometric C-alpha Rg, not SAXS-derived Rg. bAIes is an exploratory8A seqsep2 sensitivity variant, not exact paper reproduction or established correction. Its1000-frame predictor snapshot precedes the final full20000-frame analysis and is labeled separately. REST is enhanced sampling: frame adjacency and cis/trans counts are not kinetic rates. Geometry flags and trapped cis states are retained, documented in the notebook and support audits. ASTEROIDS archived selections optimized hydration separately; canonical dro5 curves are reevaluations, not new selections or Kish fits.

Code retains original source paths and hashes for provenance. Raw structures, trajectories, native VM and environments remain in the source workspace and are excluded from this curated bundle. Structural summaries, contact matrices, helix profiles and forward-model matrices are included.
'''
(D/'REPORT.txt').write_text(report)
files=[p for p in D.rglob('*') if p.is_file() and p.name!='SHA256SUMS'];(D/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(D))+'\n' for p in sorted(files)))
print(len(files),'files',sum(p.stat().st_size for p in files),'bytes',D)

from pathlib import Path
import sys,json,hashlib,base64,zipfile,html
import numpy as np,pandas as pd
from scipy.stats import rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;H=R/'outputs/Hydrophobic_SAXS_audit_2026-10-16';sys.path.insert(0,str(H));from analyze import load,C,NAMES
checks=[];fits=[]
for f in sorted((O/'motif_ablations').glob('*/complete.json')):
 a=json.loads(f.read_text());d,w,raw,q,y,sy=load(a['base_pool']);ii=a['source_indices'];ww=pd.read_csv(f.parent/'weights.csv').weight.to_numpy();pred=a['scale']*np.einsum('ij,j->i',raw[:,ii],ww);curve=pd.read_csv(f.parent/'curve.csv').prediction.to_numpy();assert np.isfinite(pred).all() and np.allclose(pred,curve,rtol=1e-10,atol=1e-10);assert np.isclose(ww.sum(),1) and (ww>=0).all();neff=1/np.einsum('i,i->',ww,ww);assert neff>=200-1e-5;mid=(q>.08)&(q<.25);chi=np.mean(((pred[mid]-y[mid])/sy[mid])**2);assert np.isclose(chi,a['chi2_mid']);checks.append({'fit':a['label'],'max_curve_difference':float(abs(pred-curve).max()),'Neff_independent':neff});fits.append({'experiment':a['label'],'N':a['N'],'Neff':neff,'chi2_mid':chi,'chi2_all':a['chi2_all']})
base=json.loads((H/'fit_experiments/ASTEROIDS_WT_baseline/complete.json').read_text());fits.insert(0,dict(experiment='AST_WT_baseline',N=base['N'],Neff=base['Neff'],chi2_mid=base['chi2_mid'],chi2_all=base['chi2_all']));pd.DataFrame(fits).to_csv(O/'ablation_results.csv',index=False)
# All finite distances and exact pair accounting; independent direct-vs-binned geometric Debye check.
fig,axes=plt.subplots(7,1,figsize=(10,20),layout='constrained');cr=[]
for row,name in enumerate(NAMES):
 z=np.load(O/(name+'.npz'));g=np.load(O/(name+'_distance_scattering.npz'));assert np.isfinite(z['xyz']).all();assert np.all(z['hist'].sum((1,2))==121*120//2)
 d,w,raw,q,y,sy=load(name);bins=np.floor(d.Rg_CA_A.to_numpy()/2).astype(int);x=z['hist'][:,2].reshape(len(d),50,5).sum(2).astype(float);v=(raw/raw[np.argmin(abs(q-.03)),None].reshape(1,-1)).T
 # Within-bin mean centering of ranks for each distance bin / q column.
 x=rankdata(x,axis=0);v=rankdata(v,axis=0)
 for b in np.unique(bins):
  mask=bins==b;x[mask]-=x[mask].mean(0);v[mask]-=v[mask].mean(0)
 den=np.sqrt((x*x).sum(0)[:,None]*(v*v).sum(0)[None,:]);corr=np.divide(np.einsum('ni,nj->ij',x,v,optimize=False),den,out=np.zeros_like(den),where=den>0)
 assert np.isfinite(corr).all() and abs(corr).max()<=1.000001
 np.savez_compressed(O/(name+'_r_q_correlations.npz'),r_A=np.arange(2.5,250,5),q=q,correlations=corr)
 im=axes[row].pcolormesh(q,np.arange(2.5,250,5),corr,cmap='coolwarm',vmin=-1,vmax=1,shading='nearest');axes[row].set(ylim=(0,100),title=name,xlabel='q (Å⁻¹)',ylabel='Long-range CA distance (Å)');axes[row].axvline(.08,color='k',ls=':');axes[row].axvline(.25,color='k',ls=':')
 for i in [0,len(d)//2]:
  xyz=z['xyz'][i];dist=np.linalg.norm(xyz[:,None]-xyz[None,:],axis=2);exact=np.sinc(q[:,None,None]*dist[None]/np.pi).mean((1,2));approx=g['geometric_CA_pair_terms'][i].sum(0)+g['self_term'];err=float(abs(exact-approx).max());checks.append(dict(ensemble=name,index=i,geometric_binning_max_absolute_I_over_I0_error=err));assert err<.02
fig.colorbar(im,ax=axes,label='Rank correlation after within-2Å-Rg-bin centering');fig.suptitle('Distance-resolved correlations with native I(q)/I(0.03)\nSequence separation >20; geometric CA histograms; no curve fitting')
for ext in ['png','pdf']:fig.savefig(O/f'distance_q_correlations.{ext}',dpi=150)
(O/'validation.json').write_text(json.dumps(checks,indent=2))
text='''Mid-q structural motif audit — 19 October 2026

Scope: 9,866 conformers across REST WT, ASTEROIDS WT/AA priors, IDPForge WT/AA and previous Jeffreys bAIes WT/AA. Common 121-residue constructs and existing canonical native PepsiSAXS predictions. New harmonic bAIes has no native SAXS yet and is not included here.

FINDINGS
1. A selected ASTEROIDS WT conformer (source index 365) matches REST index 486 at Rg 30.757 vs 30.369 A. Their native mid-q shapes, normalized at the grid point nearest q=.03, have log-RMS difference .02243. Buried nonlocal hydrophobic areas are 329 vs 1014 A^2; hydrophobic contact pairs are 3 vs 18. ASTEROIDS has 8 close antiparallel segment-chord pairs and greater anisotropy (.771 vs .571). It is an elongated, folded-back chain arrangement, rather than a densely packed hydrophobic core. See the PDBs and distance maps for exact geometry.
2. For WT coil pools, within-2A-Rg-bin ranked correlations with favorable infinitesimal mid-q enrichment are strongest for long-range CA pairs at 10–30A: IDPForge .755, bAIes .663, ASTEROIDS .829. Burial correlations are .579, .357, .549. ASTEROIDS anisotropy is .693. AA shows the same broad direction. These scores describe a perturbation around each uniform parent with its full-q scale fixed; they are not individual conformer fit scores and REST starts on a different side of the residual.
3. Broad chain reversals are not uniquely necessary. Removing 714 low-burial broad-reversal ASTEROIDS structures changes diagnostic SAXS-only mid-q chi2 from 34.95 to 37.11; Rg-matched random removals give 44.20 and 44.74. Removing the six structures with burial >=800 A^2 gives 36.66; six Rg-matched random removals give 35.89 and 35.88. This pool has multiple ways to supply useful distances; the broad reversal definition is too inclusive to identify a uniquely required motif.
4. The structural distinction to pursue is intermediate-distance crowding versus atomically tight packing. The example establishes an actual low-packing substitute for a packed conformer's native mid-q shape, and it belongs to the historical selected ASTEROIDS ensemble. Other example pairs are supplied, including candidates not selected historically.

METHODS
CA pair histograms: 1 A bins, sequence separations 1–6, 7–20, >20, all 121 residues including GP. Histograms count every unordered CA pair exactly once. The equal-weight CA Debye transform includes the self term; it is a geometric diagnostic, not a reconstruction of PepsiSAXS. Direct unbinned checks are recorded in validation.json. Representative heavy-atom P(r) uses ZiZj electron-count weights, normalized by all heavy-atom pair weight; no hydrogen, solvent, or q-dependent form factors. It is explicitly a structural fingerprint.
Packing metrics reuse the previous audit: nonlocal sidechain C/S surface burial outside +/-6 residues; 1.4 A probe; 240 SASA sphere points; hydrophobic contact cutoff 5 A, sequence separation >6. No structures were rebuilt.
Broad reversal metric: 10-residue chords starting every 5 residues; sequence-start separation >=20; chord dot product <-.7 and midpoint distance <20 A. It counts overlapping segment pairs, not independent hairpins. Local heterogeneity is the coefficient of variation of 20-residue window Rg, sampled every 5 residues. Anisotropy is relative shape anisotropy from the CA gyration tensor.
Native matching: REST burial >800 and ASTEROIDS <400 A^2, AST Rg 24–32 A, paired Rg difference <1 A; minimize RMS log native I(q)/I(qnear .03) difference over .08<q<.25. Three examples use distinct REST indices. These thresholds are exploratory and recorded, not fitted to experimental residuals. Shape normalization is labeled separately from experimental scaling.
Removal fits: same validated Borthakur SAXS-only engine and fixed REST-derived SEM as the prior diagnostic comparisons; one optimized global scale, no offset, Kish absolute Neff >=200 for all ASTEROIDS conditions. Random-removal controls match removed counts within each 2 A Rg bin, seeds17/29. These diagnostic fits are separate from the historical ASTEROIDS selection and the combined CS/SAXS fits. All fit curve/weights checks, gradients and primal/dual checks passed. Occasional Apple BLAS warning messages were followed by independent finite einsum verification.

DATA/PROVENANCE
analyze.py reads manifests and verified coordinate/prediction caches from Ensemble_Comparison_121_2026-10-15 and Hydrophobic_SAXS_audit_2026-10-16. conformer_ids in each distance_scattering NPZ match metric rows and native curve columns. Representative PDB source paths and SHA256 hashes are recorded in representatives/provenance.json. No prior results or weights were modified.
'''
(O/'FINDINGS.txt').write_text(text)
# Standalone HTML contains figures; no dependencies on external links.
imgs=['matched_motifs','representative_distance_fingerprints','burial_distance_midq','feature_correlations','distance_q_correlations'];parts=['<html><meta charset="utf-8"><title>Mid-q motifs</title><style>body{max-width:1400px;margin:40px auto;font:17px Arial;line-height:1.5}img{width:100%}pre{white-space:pre-wrap}table{border-collapse:collapse}td,th{padding:8px;border:1px solid #ddd}</style><h1>Mid-q structural motifs</h1><pre>'+html.escape(text)+'</pre>',pd.DataFrame(fits).to_html(index=False)]
for img in imgs:parts.append('<h2>'+img.replace('_',' ')+'</h2><img src="data:image/png;base64,'+base64.b64encode((O/(img+'.png')).read_bytes()).decode()+'">')
parts.append('</html>');(O/'Midq_motifs.html').write_text('\n'.join(parts))
# Portable reader notebook, with relative local file references.
cells=[{'cell_type':'markdown','metadata':{},'source':['# Mid-q structural motifs\n',text],'id':'methods'}]
code="from pathlib import Path\nimport pandas as pd\nfrom IPython.display import display, Image\np=Path.cwd()\ndisplay(pd.read_csv(p/'ablation_results.csv'))\ndisplay(pd.read_csv(p/'matched_REST_ASTEROIDS.csv').head(10))\n"
cells.append(dict(cell_type='code',metadata={},execution_count=None,outputs=[],source=[code],id='tables'))
for i,img in enumerate(imgs):cells.append(dict(cell_type='code',metadata={},execution_count=None,outputs=[],source=[f"display(Image(filename=str(p/'{img}.png')))"],id='figure'+str(i)))
(O/'Midq_motifs.ipynb').write_text(json.dumps(dict(cells=cells,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}},nbformat=4,nbformat_minor=5),indent=2))
print(pd.DataFrame(fits).to_string(index=False));print('validation and report complete')

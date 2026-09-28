"""Direct nonnegative least-squares ensemble weights; no entropy or acceptance tolerance."""
import csv,json,sys
from pathlib import Path
import numpy as np
from scipy.linalg import eigvalsh

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/Tau5_joint_refinement'

def simplex(v):
    u=np.sort(v)[::-1];c=np.cumsum(u)-1
    k=np.flatnonzero(u-c/np.arange(1,len(v)+1)>0)[-1]
    return np.maximum(v-c[k]/(k+1),0)

def solve(a,b):
    # Center first to avoid cancellation between large absolute chemical shifts.
    r=a-b; q=r@r.T/r.shape[1]
    lips=float(eigvalsh(q,subset_by_index=[len(q)-1,len(q)-1])[0])
    w=np.full(len(a),1/len(a));y=w.copy();t=1.
    for iteration in range(50000):
        new=simplex(y-(q@y)/lips)
        nt=(1+np.sqrt(1+4*t*t))/2
        y=new+(t-1)/nt*(new-w);w=new;t=nt
        if iteration%100==0:
            g=q@w;gap=float(w@g-g.min())
            if gap<1e-9:break
    g=q@w;gap=float(w@g-g.min())
    return w,dict(iterations=iteration+1,duality_gap=gap,converged=gap<1e-8)

def main(sample):
    folder=BASE/sample;out=folder/'direct_shift_fit';out.mkdir(exist_ok=True)
    # Use all fully back-calculated base and enrichment candidates available now.
    paths=sorted(p for p in (folder/'backcalc').glob('*/conformer_*/result.json')
                 if p.parent.parent.name.startswith(('pool_','helix_lengths_'))
                 and p.parent.parent.name.split('_')[-1].isdigit())
    records=[json.loads(p.read_text()) for p in paths]
    restraints=json.loads((ROOT/f'outputs/ASTEROIDS_setup/{sample}/chemical_shift_restraints.json').read_text())
    keys=set.intersection(*(set(r['shifts']) for r in records))
    restraints=[r for r in restraints if 3<=r['model_residue']<=119 and f'{r["model_residue"]}:{r["atom"]}' in keys]
    obs=np.array([r['shift_ppm'] for r in restraints])
    atoms=np.array([r['atom'] for r in restraints])
    cs=np.array([[r['shifts'][f'{x["model_residue"]}:{x["atom"]}'] for x in restraints] for r in records])
    groups=np.array([int(p.parent.parent.name.split('_')[-1])%2 for p in paths])
    q,iy,sy=np.loadtxt(ROOT/f'outputs/ASTEROIDS_setup/{sample}/{sample}_EOM_fit_range.dat',unpack=True)
    curves=np.array([np.interp(q,r['q'],r['intensity']) for r in records])
    scales={'CA':1.,'CB':1.,'C':1.,'N':2.45,'NH':.49,'HA':.25}
    sigma=np.array([scales[a] for a in atoms])
    report={'method':'nonnegative least squares, sum weights=1; no entropy penalty or discrepancy acceptance rule',
            'relative_nucleus_scales_ppm':scales,'note':'N/H scaling balances nuclei; these are not fit acceptance thresholds. Carbon-only optimum is a diagnostic, not the final all-shift fit.', 'ensembles':[]}
    weight_rows=[];shift_rows=[]
    for e,g in [(1,1),(2,0)]:
        ids=np.flatnonzero(groups==g);a=cs[ids];raw=np.full(len(ids),1/len(ids))
        for label in ['raw','carbon_only','all_shifts']:
            mask=np.isin(atoms,['CA','CB','C']) if label=='carbon_only' else np.ones(len(atoms),bool)
            if label=='raw':w=raw;status={}
            else:w,status=solve(a[:,mask]/sigma[mask],obs[mask]/sigma[mask])
            residual=w@a-obs;pred=w@curves[ids]
            scale=float(np.dot(pred/sy,iy/sy)/np.dot(pred/sy,pred/sy))
            entry=dict(ensemble=e,label=label,n_candidates=len(ids),**status,
                RMSD_ppm={k:float(np.sqrt(np.mean(residual[atoms==k]**2))) for k in sorted(set(atoms))},
                R2_RMSD_ppm={k:float(np.sqrt(np.mean(residual[(atoms==k)&np.array([390<=r['author_residue']<=414 for r in restraints])]**2))) for k in sorted(set(atoms))},
                kish_effective_size=float(1/(w@w)),nonzero_weights=int((w>1e-8).sum()),
                saxs_chi2_mean=float(np.mean(((scale*pred-iy)/sy)**2)))
            report['ensembles'].append(entry)
            for i,ww in zip(ids,w):weight_rows.append([records[i]['pdb'],e,label,ww])
            for j,r in enumerate(restraints):
                key=f'{r["model_residue"]}:{r["atom"]}'
                refs=np.array([records[i]['random_coil_shifts'][key] for i in ids])
                if np.ptp(refs)>1e-6:raise ValueError('Nonconstant random coil baseline')
                shift_rows.append([e,label,r['author_residue'],r['atom'],obs[j]-refs[0],(w@a)[j]-refs[0]])
            print(sample,e,label,entry['RMSD_ppm'],status,flush=True)
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    for name,header,rows in [('weights.csv',['source_pdb','ensemble','fit','weight'],weight_rows),('secondary_shifts.csv',['ensemble','fit','AR_residue','atom','experimental','predicted'],shift_rows)]:
        with (out/name).open('w') as f:
            writer=csv.writer(f);writer.writerow(header);writer.writerows(rows)

if __name__=='__main__':main(sys.argv[1])

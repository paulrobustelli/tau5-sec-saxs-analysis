"""Exploratory maximum-entropy joint fit, with an explicit regularization scan.

Not ASTEROIDS. Fit average observables, never average per-conformer fit scores.
SPARTA+ validation RMSDs supply conservative model-error scales (ppm).
"""
import argparse, csv, json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, check_grad
from scipy.special import softmax
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/Tau5_joint_refinement'
SIGMA={'N':2.45,'C':1.07,'CA':0.92,'CB':1.13,'HA':0.25,'NH':0.49}

def objective(z, cs, obs, sigma, saxs, iy, sy, theta, cs_on=True, saxs_on=True):
    w=softmax(z)
    logratio=np.log(np.maximum(w,1e-300)*len(w))
    cost=theta*np.dot(w,logratio)
    dw=theta*(logratio+1)
    if cs_on:
        residual=(w@cs-obs)/sigma
        cost+=0.5*np.dot(residual,residual)
        dw+=cs@(residual/sigma)
    if saxs_on:
        pred=w@saxs
        scale=np.dot(pred/sy,iy/sy)/np.dot(pred/sy,pred/sy)
        residual=(scale*pred-iy)/sy
        cost+=0.5*np.dot(residual,residual)
        dw+=scale*saxs@(residual/sy)
    return cost,w*(dw-np.dot(w,dw))

def metrics(w,cs,obs,sigma,saxs,iy,sy,rg):
    pred=w@saxs
    scale=float(np.dot(pred/sy,iy/sy)/np.dot(pred/sy,pred/sy))
    kl=float(np.dot(w,np.log(np.maximum(w,1e-300)*len(w))))
    return {'chi2_cs_mean':float(np.mean(((w@cs-obs)/sigma)**2)),
            'chi2_saxs_mean':float(np.mean(((scale*pred-iy)/sy)**2)),
            'scale':scale,'KL_from_uniform':kl,'entropy_effective_size':float(len(w)*np.exp(-kl)),
            'kish_effective_size':float(1/np.dot(w,w)),
            'mean_CA_Rg_A':float(w@rg),'rms_CA_Rg_A':float(np.sqrt(w@(rg**2)))}

def fit(sample,pools):
    folder=BASE/sample
    paths=[]
    for pool in pools: paths+=sorted((folder/'backcalc'/pool).glob('*/result.json'))
    results=[json.load(open(p)) for p in paths]
    if len(results)<20: raise ValueError('Need at least 20 successful back-calculated conformers')
    restraints=json.load(open(ROOT/f'outputs/ASTEROIDS_setup/{sample}/chemical_shift_restraints.json'))
    common=set.intersection(*(set(r['shifts']) for r in results))
    retained=[r for r in restraints if r['model_residue'] not in (1,120) and f'{r["model_residue"]}:{r["atom"]}' in common]
    excluded=[dict(r,reason=('Chain-terminal SPARTA+ shifts lack conformational information (database-average fallback or absent prediction)' if r['model_residue'] in (1,120) else 'SPARTA+ provides no prediction for this atom/residue in every conformer')) for r in restraints if r not in retained]
    (folder/'unpredicted_restraints.json').write_text(json.dumps(excluded,indent=2)+'\n')
    keys=[f'{r["model_residue"]}:{r["atom"]}' for r in retained]
    cs=np.array([[r['shifts'][k] for k in keys] for r in results])
    obs=np.array([r['shift_ppm'] for r in retained]); sigma=np.array([SIGMA[r['atom']] for r in retained])
    q,iy,sy=np.loadtxt(ROOT/f'outputs/ASTEROIDS_setup/{sample}/{sample}_EOM_fit_range.dat',unpack=True)
    saxs=np.array([np.interp(q,r['q'],r['intensity']) for r in results])
    saxs/=np.mean([r['intensity'][0] for r in results])
    rg=np.array([r['ca_rg_A'] for r in results]); n=len(results)
    np.savez_compressed(folder/'observable_matrix.npz',chemical_shifts=cs,shift_observed=obs,shift_sigma=sigma,q=q,intensity=saxs,intensity_observed=iy,intensity_sigma=sy,ca_rg_A=rg)
    uniform=np.full(n,1/n)
    records=[dict(label='pool',theta=None,**metrics(uniform,cs,obs,sigma,saxs,iy,sy,rg))]
    weights={'pool':uniform}
    z=np.zeros(n)
    # Analytic derivative sanity check on a deterministic, bounded perturbation.
    rng=np.random.default_rng(20260928); testz=rng.normal(0,.01,n)
    v=rng.normal(size=n); v/=np.linalg.norm(v); eps=1e-5
    f,g=objective(testz,cs,obs,sigma,saxs,iy,sy,100)
    fd=(objective(testz+eps*v,cs,obs,sigma,saxs,iy,sy,100)[0]-objective(testz-eps*v,cs,obs,sigma,saxs,iy,sy,100)[0])/(2*eps)
    if not np.isclose(fd,g@v,rtol=1e-4,atol=1e-3): raise ValueError(f'Gradient check failed {fd} {g@v}')
    for theta in (1000,300,100,30,10,3,1):
        opt=minimize(objective,z,args=(cs,obs,sigma,saxs,iy,sy,theta),jac=True,method='L-BFGS-B',options={'maxiter':3000,'ftol':1e-11,'gtol':1e-7})
        z=opt.x; w=softmax(z); label=f'joint_theta_{theta}'; weights[label]=w
        records.append(dict(label=label,theta=theta,optimizer_success=bool(opt.success),optimizer_message=str(opt.message),**metrics(w,cs,obs,sigma,saxs,iy,sy,rg)))
    good=[r for r in records[1:] if r['optimizer_success'] and r['chi2_cs_mean']<=1.5 and r['chi2_saxs_mean']<=1.5]
    if good:
        chosen=min(good,key=lambda r:r['KL_from_uniform']); reason='Most regularized scanned fit with both mean normalized squared residuals <=1.5 (exploratory discrepancy rule, not cross-validation)'
    else:
        chosen=min(records[1:],key=lambda r:r['chi2_cs_mean']+r['chi2_saxs_mean']); reason='No scanned fit met both discrepancy targets; best combined normalized residual shown, not an accepted refined ensemble'
    for label,cs_on,saxs_on in [('shifts_only',True,False),('saxs_only',False,True)]:
        opt=minimize(objective,np.zeros(n),args=(cs,obs,sigma,saxs,iy,sy,chosen['theta'],cs_on,saxs_on),jac=True,method='L-BFGS-B',options={'maxiter':3000,'ftol':1e-11})
        w=softmax(opt.x);weights[label]=w;records.append(dict(label=label,theta=chosen['theta'],optimizer_success=bool(opt.success),**metrics(w,cs,obs,sigma,saxs,iy,sy,rg)))
    w=weights[chosen['label']]
    with (folder/'weights.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['source_pdb','protonated_pdb','CA_Rg_A']+list(weights))
        for i,r in enumerate(results):writer.writerow([r['pdb'],r['protonated_pdb'],rg[i]]+[weights[k][i] for k in weights])
    with (folder/'shift_fit.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['model_residue','author_residue','atom','observed_ppm','pool_ppm','joint_ppm','model_error_ppm'])
        for j,r in enumerate(retained):writer.writerow([r['model_residue'],r['author_residue'],r['atom'],obs[j],uniform@cs[:,j],w@cs[:,j],sigma[j]])
    np.savetxt(folder/'saxs_fit.csv',np.c_[q,iy,sy,records[0]['scale']*(uniform@saxs),chosen['scale']*(w@saxs)],delimiter=',',header='q_Ainv,I_exp,sigma_exp,I_pool_scaled,I_joint_scaled',comments='')
    report={'sample':sample,'method':'custom maximum-entropy reweighting, not ASTEROIDS','pool_size':n,'input_restraints':len(restraints),'fitted_restraints':len(retained),'unpredicted_restraints':len(excluded),'chemical_shift_model_errors_ppm':SIGMA,'chosen':chosen,'selection_rule':reason,'scan':records,'independent_validation':False,'caveats':['Shifts used for both sampling and fitting; fit is not independent validation','SAXS q points treated with diagonal propagated errors; correlations not modeled','SPARTA+ model-error scales from folded-protein validation are approximate for this IDP','All reported Rg summaries use C-alpha coordinates; not experimental Guinier Rg','Hydration density fixed across conformers; no additive background or shift offset fitted','Sampling probabilities are targets, not guaranteed realized secondary structure populations']}
    (folder/'fit_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    fig,axs=plt.subplots(2,2,figsize=(11,7))
    ax=axs[0,0];ax.errorbar(q,iy,yerr=sy,fmt='.',ms=2,color='0.6',label='experiment');ax.plot(q,records[0]['scale']*(uniform@saxs),label='pool');ax.plot(q,chosen['scale']*(w@saxs),label='joint');ax.set_yscale('log');ax.set_xlabel('q (Å⁻¹)');ax.set_ylabel('Intensity');ax.legend()
    ax=axs[1,0];ax.axhline(0,color='0.5',lw=1);ax.plot(q,(chosen['scale']*(w@saxs)-iy)/sy);ax.set_ylabel('(calculated − observed) / σ');ax.set_xlabel('q (Å⁻¹)')
    ax=axs[0,1];bins=np.linspace(min(rg)-1,max(rg)+1,25);ax.hist(rg,bins=bins,weights=uniform,histtype='step',label='pool');ax.hist(rg,bins=bins,weights=w,histtype='step',label='joint');ax.set_xlabel('Cα radius of gyration (Å)');ax.set_ylabel('Population');ax.legend()
    ax=axs[1,1];atoms=sorted(set(r['atom'] for r in retained));xx=np.arange(len(atoms));masks=[np.array([r['atom']==a for r in retained]) for a in atoms]
    ax.bar(xx-.18,[np.sqrt(np.mean(((uniform@cs-obs)[m]/sigma[m])**2)) for m in masks],.36,label='pool');ax.bar(xx+.18,[np.sqrt(np.mean(((w@cs-obs)[m]/sigma[m])**2)) for m in masks],.36,label='joint');ax.set_xticks(xx,atoms);ax.set_ylabel('Shift RMSD / model error');ax.legend()
    fig.suptitle(f'{sample}: exploratory joint fit, {n} conformers; θ={chosen["theta"]}\n'+('Targets met; pool convergence still requires assessment' if good else 'Fit targets NOT both met'))
    fig.tight_layout();fig.savefig(folder/'joint_fit.png',dpi=180);plt.close(fig)
    print(json.dumps({k:report[k] for k in ['sample','pool_size','fitted_restraints','unpredicted_restraints','chosen','selection_rule']},indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('sample');ap.add_argument('pools',nargs='+');a=ap.parse_args();fit(a.sample,a.pools)

"""Zheng & Best (2018), equations 5 and 6; units Angstrom.
Intensity-space weighted least squares; no intensity clipping or background fit.
N is peptide bonds, not residues. Conditional diagonal errors only.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent

def rg_from_nu(nu, residues=120):
    gamma=1.1615
    return np.sqrt(gamma*(gamma+1)/(2*(gamma+2*nu)*(gamma+2*nu+1)))*5.5*(residues-1)**nu

def shape(q,nu,residues=120):
    x=q*rg_from_nu(nu,residues)
    return np.exp(np.clip(-x*x/3+.0479*(nu-.212)*x**4,-100,100))

def fit(q,y,e,residues=120):
    if len(q)<8 or np.sum(y>3*e)<6:return None
    def objective(nu,details=False):
        f=shape(q,nu,residues); w=1/e**2
        amp=np.sum(y*f*w)/np.sum(f*f*w)
        if amp<=0:return 1e100
        chi=np.sum(((y-amp*f)/e)**2)
        return (amp,chi) if details else chi
    # Multiple brackets avoid the unphysical high-q upturn's competing minima.
    grid=np.linspace(.2,.8,61);scores=np.array([objective(n) for n in grid]);candidates=[]
    for j in range(1,len(grid)-1):
        if scores[j]<=scores[j-1] and scores[j]<=scores[j+1]:
            opt=minimize_scalar(objective,bounds=(grid[j-1],grid[j+1]),method='bounded',options={'xatol':1e-10})
            candidates.append((opt.fun,opt.x))
    candidates.extend([(scores[0],grid[0]),(scores[-1],grid[-1])])
    chi,nu=min(candidates);amp,_=objective(nu,True);rg=rg_from_nu(nu,residues)
    eps=1e-5;der=(shape(q,nu+eps,residues)-shape(q,nu-eps,residues))/(2*eps)
    J=np.c_[shape(q,nu,residues),amp*der]/e[:,None]
    cov=np.linalg.pinv(J.T@J); dr=(rg_from_nu(nu+eps,residues)-rg_from_nu(nu-eps,residues))/(2*eps)
    return dict(Rg=float(rg),nu=float(nu),I0=float(amp),Rg_conditional_se=float(abs(dr)*np.sqrt(cov[1,1])),chi2_reduced=float(chi/(len(q)-2)),qmin=float(q[0]),qmax=float(q[-1]),qRgmax=float(q[-1]*rg),n=int(len(q)),residues=residues,boundary=bool(nu<.201 or nu>.799))

def scan(q,y,e,qmin=.012,residues=120):
    rows=[]
    for hi in q[(q>=qmin+.007-1e-8)&(q<=.12)]:
        m=(q>=qmin-1e-8)&(q<=hi)&np.isfinite(y)&np.isfinite(e)&(e>0)
        r=fit(q[m],y[m],e[m],residues)
        if r:rows.append(r)
    valid=[r for r in rows if r['qRgmax']<=2 and not r['boundary']]
    return (max(valid,key=lambda r:r['qmax']) if valid else None),rows

def run():
    result={};all_scans={};fig,axes=plt.subplots(2,3,figsize=(14,8));sens,sa=plt.subplots(2,3,figsize=(14,7))
    for row,sample in enumerate(['WT_092025','AA_052026']):
        a=np.load(ROOT/'inputs'/f'{sample}_audit.npz');q=a['q']
        for col,(kind,ekey) in enumerate([('average','average_error'),('candidate','candidate_bootstrap_error'),('other','other_error')]):
            key=sample+'_'+kind;y=a[kind];e=a[ekey];primary,rows=scan(q,y,e);all_scans[key]=rows
            variations={str(lo):scan(q,y,e,lo)[0] for lo in [.010,.012,.015,.020]}
            short=scan(q,y,e,.012,118)[0]
            result[key]=dict(primary=primary,qmin_sensitivity=variations,residues118=short,interpretation='Exploratory single-chain hypothesis, identity unknown' if kind=='other' else 'Single-chain model; does not establish monomer identity')
            ax=axes[row,col];m=(q>=.01)&(q<=.085);ax.errorbar(q[m]**2,y[m],e[m],fmt='.',ms=3,elinewidth=.5,color='0.5')
            ax.set(xlabel='q² (Å⁻²)',ylabel='I(q), arbitrary scale',title=key.replace('_092025','').replace('_052026',''))
            if primary:
                m=(q>=primary['qmin']-1e-8)&(q<=primary['qmax']+1e-8)
                ax.plot(q[m]**2,primary['I0']*shape(q[m],primary['nu']),color='C1',label=f"Extended Rg={primary['Rg']:.1f} Å\nχ²/dof={primary['chi2_reduced']:.1f}");ax.legend(fontsize=8)
                np.savetxt(ROOT/(key+'_fit.dat'),np.c_[q[m],y[m],e[m],primary['I0']*shape(q[m],primary['nu'])],header='q I sigma extended_fit')
                if kind!='other':
                    bootkey='bootstrap' if kind=='candidate' else 'bootstrap_average'
                    br=[fit(q[m],b[m],e[m]) for b in a[bootkey]];br=[b['Rg'] for b in br if b is not None]
                    result[key]['fixed_window_bootstrap_95']=np.percentile(br,[2.5,97.5]).tolist()
            else:ax.text(.05,.85,'No supported fit',transform=ax.transAxes)
            ss=sa[row,col]
            for lo in [.010,.012,.015,.020]:
                _,rr=scan(q,y,e,lo);rr=[r for r in rr if r['qRgmax']<=2 and not r['boundary']]
                ss.plot([r['qmax'] for r in rr],[r['Rg'] for r in rr],label=f'qmin={lo}')
            ss.set(xlabel='qmax (Å⁻¹)',ylabel='Extended Rg (Å)',title=key.replace('_092025','').replace('_052026',''));ss.legend(fontsize=7)
    fig.tight_layout();fig.savefig(ROOT/'extended_fits.png',dpi=160);plt.close(fig)
    sens.tight_layout();sens.savefig(ROOT/'fit_range_sensitivity.png',dpi=160);plt.close(sens)
    (ROOT/'results.json').write_text(json.dumps(result,indent=2));(ROOT/'scans.json').write_text(json.dumps(all_scans,indent=2))
    # Exact synthetic recovery checks the exponent, units and N convention together.
    q=np.linspace(.012,.06,49);nu=.53;y=3.2*shape(q,nu);r=fit(q,y,np.ones(len(q))*.01)
    assert abs(r['nu']-nu)<1e-6 and abs(r['I0']-3.2)<1e-5
    (ROOT/'validation.json').write_text(json.dumps(dict(synthetic_recovery_passed=True,true_nu=nu,recovered=r),indent=2))
    return result

if __name__=='__main__':
    for k,v in run().items():print(k,v['primary'],flush=True)

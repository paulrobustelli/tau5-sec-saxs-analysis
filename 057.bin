from build_audit_notebooks import write,setup
write('03_report_EFA_sensitivity.ipynb',[
('m',r'''# Check the EFA settings in the new RAW reports
This notebook tests the *reported component hypotheses*, not assumed biological species counts. It repeats forward/backward EFA and constrained rotation over both the full recorded q range and q=0.012–0.25 Å⁻¹. The report's selected profiles are not supplied, so this is a reconstruction of its stated settings from the repository HDF5 series, **not an exact reproduction of the report's RAW session**.

WT: overall frames 130–200; four proposed supports 130–154, 130–194, 154–195, 155–195. AA: overall frames 117–250; three supports 160–222, 168–197, 168–222. All bounds are inclusive.

The PDF WT buffer range (210–239) differs from the notebook (220–236). Both reports list endpoint linear baseline windows 0–10 and 439–449. Here we reconstruct a weighted line after constant residual-buffer subtraction, and propagate its prediction variance. Since fitting a line includes an intercept, the constant subtraction cancels in the final intensities for fixed weights; its effect on weighting remains. Facility-subtraction covariance cannot be reconstructed from the supplied data.'''),('c',setup+"\nimport raw_efa_core as efa\n"),
('m',r'''## Reconstruct the stated background model
For q bin i, fit $y_i(t)=\alpha_i+\beta_i t$ to the endpoint windows. With $x_t=(1,t)^T$, the corrected mean is $y_i(t)-x_t^T\hat\beta_i$. The added baseline variance is $x_t^T\mathrm{Cov}(\hat\beta_i)x_t$, and the covariance between two corrected frames includes $x_t^T\mathrm{Cov}(\hat\beta_i)x_u$. SVD row scaling below is an approximation to noise weighting, not full whitening of that covariance.'''),
('c',"""def report_background(name):
    q,Y,E=load(name)
    a,b=(210,239) if name.startswith('WT') else (24,70)
    bg=Y[:,a:b+1].mean(1)
    bgv=(E[:,a:b+1]**2).sum(1)/(b-a+1)**2
    Y0=Y-bg[:,None]; E0=np.sqrt(E**2+bgv[:,None])
    idx=np.r_[0:11,439:450]; X=np.c_[np.ones(len(idx)),idx]
    T=np.c_[np.ones(Y.shape[1]),np.arange(Y.shape[1])]
    out=Y0.copy(); err=E.copy()
    for j in range(len(q)):
        w=1/E0[j,idx]**2
        cov=np.linalg.inv(X.T@(w[:,None]*X))
        beta=cov@(X.T@(w*Y0[j,idx]))
        out[j]-=T@beta
        err[j]=np.sqrt(E[j]**2+np.einsum('ij,jk,ik->i',T,cov,T))
    return q,out,err
specs=[('WT_092025',130,200,[[130,154],[130,194],[154,195],[155,195]]),
       ('AA_052026',117,250,[[160,222],[168,197],[168,222]])]
"""),
('m',r'''## SVD, evolving factors, and rotation
$$A_{it}=D_{it}/\bar\sigma_i,\qquad A=U\Sigma V^T.$$
Forward EFA computes the singular values of $A_{:,0:t}$; backward EFA uses $A_{:,t:T}$. Supports inferred from these curves constrain component elution profiles. A rotation of the retained subspace can give $A\approx S_w C^T$, but overlapping supports can make that rotation poorly determined. Smooth singular vectors, residual structure, positivity, and stability across ranges all matter.

The code below uses RAW's numerical EFA/rotation functions, included with their GPL license. It records convergence and negative scattering values, instead of equating convergence with physical validity.'''),
('c',"""results=[]
for name,a,b,supports in specs:
    q,Y,E=report_background(name)
    fig,axes=plt.subplots(2,3,figsize=(13,7))
    for row,(lo,hi) in enumerate([(float(q[0]),float(q[-1])),(.012,.25)]):
        m=(q>=lo)&(q<=hi); err=E[m,a:b+1]; w=err.mean(1)
        D=Y[m,a:b+1]/w[:,None]
        U,s,Vt=np.linalg.svd(D,full_matrices=False)
        fw=efa.runEFA(D); bw=efa.runEFA(D,False)[:,::-1]
        frames=np.arange(a,b+1)
        axes[row,0].semilogy(np.arange(1,9),s[:8],'o-')
        axes[row,0].set(title=f'q={lo:g}–{hi:g} Å⁻¹',xlabel='Mode',ylabel='Singular value')
        for k in range(min(5,len(s))):
            axes[row,1].semilogy(frames,fw[k],label=str(k+1))
            axes[row,2].semilogy(frames,bw[k],label=str(k+1))
        axes[row,1].set(title='Forward EFA',xlabel='Frame',ylabel='Singular value')
        axes[row,2].set(title='Backward EFA',xlabel='Frame',ylabel='Singular value')
        for ax in axes[row,1:]: ax.legend(fontsize=7,ncol=3)
        with np.errstate(all='ignore'):
            ok,conv,rot=efa.runRotation(D,Y[m,a:b+1],err,np.array(supports)-a,
                                        [True]*len(supports),Vt.T,niter=3000,tol=1e-9)
        r=dict(sample=name,q=[lo,hi],frames=[a,b],supports=supports,
               singular_values=s[:8].tolist(),converged=bool(ok),iterations=conv['iterations'])
        if ok:
            C=rot['M']*rot['C']; S=Y[m,a:b+1]@np.linalg.pinv(C.T)
            r.update(concentration_condition=float(np.linalg.cond(C)),
                     negative_profile_fraction=np.mean(S<0,axis=0).tolist(),
                     residual_ms=float(np.mean(((Y[m,a:b+1]-S@C.T)/err)**2)))
            np.savez_compressed(ROOT/'results'/f'{name}_report_rotation_{row}.npz',q=q[m],frames=frames,C=C,S=S,error=err)
        results.append(r)
        print(r)
    fig.suptitle(name+' — report hypotheses, not validated species')
    fig.tight_layout(); fig.savefig(ROOT/'results'/f'{name}_report_efa.png',dpi=170)
    display(fig);plt.close(fig)
(ROOT/'results/report_EFA_sensitivity.json').write_text(json.dumps(results,indent=2));
"""),
('m',r'''## Inspect the recovered curves
Negative values can be ordinary high-q noise, but large negative fractions or ill-conditioned, nearly collinear elution profiles undermine a unique species assignment. Full-q backgrounds may introduce modes unrelated to protein conformations. Restricting q also loses information, so it is a sensitivity check, not a guarantee of a correct answer.'''),
('c',"""for r in results:
    if not r['converged']: continue
    row=0 if r['q'][1]>.3 else 1
    d=np.load(ROOT/'results'/f\"{r['sample']}_report_rotation_{row}.npz\")
    fig,ax=plt.subplots(1,2,figsize=(10,3))
    for k in range(d['C'].shape[1]):
        ax[0].plot(d['frames'],d['C'][:,k],label=f'Component {k+1}')
        ax[1].plot(d['q'],d['S'][:,k],label=f'Component {k+1}')
    ax[0].set(xlabel='Frame',ylabel='Relative coefficient');ax[0].legend(fontsize=8)
    ax[1].set(xlabel='q (Å⁻¹)',ylabel='Recovered I (arbitrary scale)'); ax[1].set_yscale('symlog',linthresh=.1)
    fig.suptitle(f\"{r['sample']} q={r['q']}: exploratory constrained rotation\")
    fig.tight_layout();display(fig);plt.close(fig)
"""),
('m',r'''## What this can establish
The source reports contain proposed component counts. They do not provide independent molecular weights. A species assignment requires stable separation with plausible scattering profiles, defensible uncertainties and additional evidence. EFA cannot distinguish monomer and oligomer when their elution is identical or their scattering curves are indistinguishable in the measured range.

For a mixture, the low-q slope gives
$$R_{g,\mathrm{app}}^2=\frac{\sum_k I_k(0)R_{g,k}^2}{\sum_k I_k(0)}.$$
Consequently, an Rg value alone cannot determine oligomer fractions. Read the clean-region/noise sensitivity results in notebook 02 before interpreting the constrained rotations here.''')])

from pathlib import Path
import nbformat as nb
ROOT=Path(__file__).resolve().parents[1]
def write(name,cells):
 n=nb.v4.new_notebook(cells=[(nb.v4.new_markdown_cell(s) if kind=='m' else nb.v4.new_code_cell(s)) for kind,s in cells],metadata={'kernelspec':{'name':'python3','display_name':'Python 3','language':'python'}});nb.write(n,ROOT/name)
setup="""from pathlib import Path
import sys, json, hashlib
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, Markdown, Image
ROOT=Path.cwd()
sys.path.insert(0,str(ROOT/'analysis'))
from inspect_data import load, correct
from select_windows import average
from run_efa import guinier
%matplotlib inline
"""
write('01_reproduce_and_audit_Rg.ipynb',[
('m',r'''# Reproduce Natalie first, then audit each change
This notebook starts from Natalie's saved per-frame profiles and reproduces **24.90 Å WT / 30.78 Å AA**. It then changes one choice at a time. The source notebooks and HDF5 files are not modified. Frames and q indices are zero based, with both endpoints included.

A changed fit convention is a methodological choice, not evidence that Natalie's result is wrong. q-range sensitivity, background covariance, and a newly selected window must be distinguished. P(r) and Guinier estimates are compared separately.'''),('c',setup),
('m','## 1. Data provenance and exact averaging\nThe saved DAT files contain Natalie’s processed intensities and errors. An arithmetic mean has variance sum(sigma²)/n² when frames are treated as independent; we first reproduce that convention.'),
('c',"""saved={}; hashes={}
for name,a,b in [('WT_092025',175,185),('AA_052026',172,182)]:
    paths=[ROOT/'BNL_SAXS/dat_files'/name/f'profile_{j}.dat' for j in range(a,b+1)]
    d=np.array([np.loadtxt(p) for p in paths])
    saved[name]=(d[0,:,0],d[:,:,1].mean(0),np.sqrt((d[:,:,2]**2).sum(0))/len(d))
    for p in paths: hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    print(name, 'frames',a,b,'number',len(d))
print('All input file hashes retained in results/reproduction_input_hashes.json')
(ROOT/'results/reproduction_input_hashes.json').write_text(json.dumps(hashes,indent=2));
"""),
('m',r'''## 2. Show the fit equation and reproduce the reported values
$$y_i=\ln I_i=a+bq_i^2,\quad \sigma_{y_i}=\sigma_{I_i}/I_i,$$
$$\hat\beta=(X^TWX)^{-1}X^TWy,\quad R_g=\sqrt{-3b},\quad I(0)=e^a.$$
This is RAW's error-weighted linear Guinier fit. The uncertainty below is the conditional covariance error; RAW can additionally report the larger fit-range variation error. WT uses q indices 3–45; AA uses 3–41, verified against saved notebook outputs.'''),
('c',"""def fixed_fit(q,I,E,lo,hi):
    m=(q>=lo-1e-12)&(q<=hi+1e-12)&(I>0)&(E>0)
    X=np.column_stack([np.ones(m.sum()),q[m]**2]); y=np.log(I[m]); sigma=E[m]/I[m]
    cov=np.linalg.inv(X.T@((1/sigma**2)[:,None]*X))
    beta=cov@(X.T@(y/sigma**2)); rg=np.sqrt(-3*beta[1])
    return dict(Rg=float(rg),fit_se=float(3*np.sqrt(cov[1,1])/(2*rg)),
                qmin=float(q[m][0]),qmax=float(q[m][-1]),qRgmax=float(q[m][-1]*rg),
                chi2=float(np.sum(((y-X@beta)/sigma)**2)/(len(y)-2)),I0=float(np.exp(beta[0])))
reproduced={}
for name,end,target in [('WT_092025',45,24.90414304),('AA_052026',41,30.78492264)]:
    q,I,E=saved[name]; f=fixed_fit(q,I,E,q[3],q[end]); reproduced[name]=f
    assert abs(f['Rg']-target)<1e-6
    print(name, f)
"""),
('m',r'''## 3. One change at a time
Stage A reproduces Natalie. B changes only the q fit convention (start 0.012 Å⁻¹, initial upper bound 0.04 Å⁻¹, iteratively reduce until qRg≤1). This conservative iteration can retain a slightly shorter interval than the longest admissible prefix. C retains exactly B's q points and intensities but propagates shared background covariance. D updates the residual-background estimator while keeping those same q points. E changes the averaging window only. F reapplies the conservative q-range rule to the selected window.

The WT corrected intensities reproduce the DAT values to their saved precision. AA differs slightly because the new weighted baseline fit uses original HDF5 errors rather than errors after residual-buffer subtraction. These are numerical changes with explicit assumptions, not an independently established correction to the instrument data.'''),
('c',"""rows=[]
for name,a,b,newa,newb in [('WT_092025',175,185,159,175),('AA_052026',172,182,169,188)]:
    q,I,E=saved[name]; B=guinier(q,I,E); lo,hi=B['qmin'],B['qmax']
    q2,I2,E2=average(name,a,b); assert np.allclose(q,q2)
    q3,I3,E3=average(name,newa,newb)
    stages=[('A saved profiles, original fit',reproduced[name]),
            ('B change q range only',B),
            ('C shared background error only',fixed_fit(q,I,E2,lo,hi)),
            ('D updated baseline estimator',fixed_fit(q,I2,E2,lo,hi)),
            ('E selected frames, same q points',fixed_fit(q3,I3,E3,lo,hi)),
            ('F selected frames, qRg rule',guinier(q3,I3,E3))]
    rows.extend(dict(sample=name,stage=label,**f) for label,f in stages)
    print(name, 'maximum profile change at original window:',np.max(abs(I-I2)))
header='| Sample | Stage | Rg (Å) | conditional SE | qmin–qmax (Å⁻¹) | qmax Rg |\\n|---|---|---:|---:|---|---:|\\n'
display(Markdown(header+'\\n'.join(f\"|{r['sample']}|{r['stage']}|{r['Rg']:.3f}|{r['fit_se']:.3f}|{r['qmin']:.3f}–{r['qmax']:.3f}|{r['qRgmax']:.3f}|\" for r in rows)))
(ROOT/'results/Rg_correction_audit.json').write_text(json.dumps(rows,indent=2));
"""),
('m','## 4. Observe q-range sensitivity directly\nThe data and errors are unchanged across each curve below. Curvature outside the lowest-q regime can change the fitted slope, particularly for the expanded AA sample. Shorter fits also have larger statistical uncertainty.'),
('c',"""fig,axes=plt.subplots(2,2,figsize=(11,8))
for col,name in enumerate(saved):
    q,I,E=saved[name]; ref=reproduced[name]; narrow=guinier(q,I,E)
    m=(q>=.008)&(q<=.06)&(I>0)
    ax=axes[0,col]; ax.errorbar(q[m]**2,np.log(I[m]),yerr=E[m]/I[m],fmt='.',alpha=.7,color='black')
    for f,label in [(ref,'Natalie reproduction'),(narrow,'Conservative qRg≤1')]:
        x=np.array([f['qmin'],f['qmax']])**2
        ax.plot(x,np.log(f['I0'])-f['Rg']**2*x/3,label=f\"{label}: {f['Rg']:.2f} Å\")
    ax.set(title=name,xlabel='q² (Å⁻²)',ylabel='ln I'); ax.legend(fontsize=8)
    ends=q[(q>=.026)&(q<=.055)]; fits=[fixed_fit(q,I,E,q[3],end) for end in ends]
    axes[1,col].errorbar(ends,[f['Rg'] for f in fits],yerr=[f['fit_se'] for f in fits],fmt='o-',markersize=3)
    axes[1,col].set(xlabel='Upper fit q (Å⁻¹)',ylabel='Guinier Rg (Å)',title='Same average, changing upper q only')
fig.tight_layout();fig.savefig(ROOT/'poster_results/Rg_reproduction_and_qrange.png',dpi=200)
display(fig);plt.close(fig)
"""),
('m',r'''## 5. P(r) is a separate estimate
Natalie's supplied poster reports GNOM Rg **27.3 Å WT / 34.23 Å AA**, with Dmax **95 / 139 Å**. Our current BIFT estimates are about **27.22 / 34.60 Å** with Dmax **94 / 159 Å**. The AA Guinier discrepancy is largely reproduced by changing the fitted q range on the *same data*. The P(r) Rg estimates are already close, but the AA Dmax difference remains material. GNOM and BIFT regularization/Dmax selection differ; neither should be silently substituted for the other. We have not rerun GNOM or established a unique AA Dmax.

A 35 Å Rg is not by itself evidence of oligomerization in an IDR. The size of this sequence must be assessed with the full scattering curve, concentration behavior, uncertainty and ideally an independent mass measurement. No window was selected by closeness to a desired Rg.

## New report context
WT_raw_report.pdf records buffer frames 210–239 plus a linear baseline using 0–10 and 439–449; this differs from the source DAT correction (220–236). Its four proposed component intervals overlap strongly; the report does not display recovered component curves. AA_raw_report.pdf shows a three-component rotation over 117–250, including q out to 3.19 Å⁻¹; two recovered curves are very noisy. These are analysis settings/results, not proof of four WT or three AA biological species. Notebook 03 checks their sensitivity.''')])
# Fix all literal LaTeX strings in the earlier notebook generator.
p=ROOT/'analysis/build_notebook.py';s=p.read_text().replace("md('''", "md(r'''")
s=s.replace("This executed notebook accompanies", "This notebook accompanies")
s=s.replace("'name':'saxs'", "'name':'python3'")
s=s.replace("ROOT/'EFA_analysis.ipynb'", "ROOT/'02_windows_profiles_and_EFA.ipynb'")
p.write_text(s)

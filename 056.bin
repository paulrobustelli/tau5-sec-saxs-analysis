from pathlib import Path
import nbformat as nb
ROOT=Path(__file__).resolve().parents[1]
c=[]
def md(s):c.append(nb.v4.new_markdown_cell(s))
def code(s):c.append(nb.v4.new_code_cell(s))
md(r'''# Tau-5* SEC-SAXS: clean-window selection, standard plots and EFA

This notebook accompanies the files in `poster_results/`. It uses the requested WT_092025 and AA_052026 HDF5 series. The original source notebooks remain unchanged. The selected windows are **WT 159–175** and **AA 169–188**, inclusive and zero indexed. See `EFA_report.md` for biological interpretation and `EFA_equations.md` for the full equation walkthrough.

No oligomeric identity is assigned from Rg alone. Window scores are a presentation-quality criterion conditional on the supplied errors and chosen residual-background correction.''')
code('''from pathlib import Path
import sys, json, numpy as np
ROOT = Path.cwd()
sys.path.insert(0, str(ROOT / 'analysis'))
from run_efa import analyze, guinier
from select_windows import average, candidates
from IPython.display import display, Image, Markdown
assert (ROOT / 'BNL_SAXS/data/WT_092025.hdf5').exists()
''')
md(r'''## Verify the WT construct

Your exact supplied sequence is the experimental reference. It has 120 residues, including the N-terminal GP; the segment after GP has 118. BMRB53105 is a 119-residue Tau-5* construct. The user sequence adds GP and removes its terminal C. User W70/W106 correspond to AR W397/W433. Expected WT monomer mass is the user-supplied 12.1 kDa.''')
code('''sequence_info = json.loads((ROOT/'background/sequence_comparison.json').read_text())
sequence_info
''')
md(r'''## Recompute the selected averaging windows

Candidates contain the local SEC maximum, have at least seven frames, and lie inside WT 150–195 or AA 155–215. Passing criteria: ≤5% Rg span over four low-q cutoffs, Guinier reduced weighted residual ≤1.5, and amplitude-scaled first/second-half profile residual mean square ≤1.5. Among passing candidates we maximize the mean-profile signal/error norm. Shared residual-background uncertainty is retained when averaging.

The optimization does not use a target Rg. Neighboring windows are nearly tied.''')
code('''selection = []
for name, lo, hi in [('WT_092025',150,195), ('AA_052026',155,215)]:
    rows, best = candidates(name, lo, hi)
    selection.append(best)
    print(name, 'evaluated:',len(rows), 'selected:',best['start'],best['end'],
          'Rg:',round(best['Rg'],3), 'A', 'qmin span:',round(100*best['qmin_Rg_relative_span'],2),'%')
''')
code('''display(Image(filename=str(ROOT/'poster_results/SEC_selected_windows.png')))
''')
md(r'''## Arithmetic mean profiles and Guinier fits

The Guinier relation is

$$\ln I(q)=\ln I(0)-\frac{R_g^2}{3}q^2,$$

so a fitted slope m gives $R_g=\sqrt{-3m}$. The primary fit starts at q=0.012 Å⁻¹ and enforces qmax Rg≤1. Points excluded from the fitted range are gray in the plots. Fit errors are conditional and exclude background-choice systematics.''')
code('''guinier_results = {}
for s in selection:
    q, intensity, error = average(s['sample'],s['start'],s['end'])
    fit = guinier(q,intensity,error)
    guinier_results[s['sample']] = fit
    print(s['sample'], fit)
''')
md(r'''## Bayesian indirect Fourier transform for P(r)

The forward relation and real-space size are

$$I(q)=4\pi\int_0^{D_{\max}}P(r)\frac{\sin(qr)}{qr}\,dr,$$

$$R_g^2=\frac{\int r^2P(r)\,dr}{2\int P(r)\,dr}.$$

This cell reruns RAW's actual BIFT implementation. BIFT selects smoothing and Dmax using its Bayesian evidence procedure. The plotted shading is its conditional parameter-sampling uncertainty, not a full systematic uncertainty. q and window sensitivity results are saved separately.''')
code('''from run_bift import fit_bift
bift_results = {}
for s in selection:
    p = fit_bift(s['sample'],s['start'],s['end'])
    bift_results[s['sample']] = p
    a = np.load(ROOT/'poster_results'/f"{s['sample']}_bift_primary.npz")
    rg_check = np.sqrt(np.trapezoid(a['p']*a['r']**2,a['r'])/(2*np.trapezoid(a['p'],a['r'])))
    assert np.isclose(rg_check,p['rg'])
    assert np.isclose(4*np.pi*np.trapezoid(a['p'],a['r']),p['i0'])
    print('Verified P(r) integrals:',s['sample'],rg_check)
''')
code('''display(Image(filename=str(ROOT/'poster_results/poster_saxs_panel.png')))
display(Image(filename=str(ROOT/'poster_results/fit_and_window_diagnostics.png')))
''')
md(r'''## EFA: independent shape patterns across elution

We model the background-corrected matrix as $D=SC^T+E$, where S contains component scattering curves and C their elution amplitudes. Divide each q row by its mean error to obtain A, then calculate $A=U\Sigma V^T$. EFA repeats SVD on progressively longer leading and trailing frame blocks. Component time windows constrain a rotation toward nonnegative elution profiles; mathematical convergence alone does not prove physical separation.

The following recomputes the clean-region and shoulder-inclusive SVD tests. It does not force a second species into the clean regions.''')
code('''efa_results = []
for name,a,b in [('WT_092025',155,180),('AA_052026',163,204),('WT_092025',125,195)]:
    r,_ = analyze(name,a,b,simulate=100)
    efa_results.append(r)
    print(name,(a,b),'s2 =',round(r['sv'][1],2),
          '99% rank-one null =',round(r['null_s2_99'],2),
          'q/time mode-2 smoothness =',round(r['q_ac'][1],2),round(r['t_ac'][1],2))
''')
code('''display(Image(filename=str(ROOT/'results/WT_092025_120-205_default_0.012-0.25_diagnostics.png')))
''')
md(r'''## Interpretation and uncertainty

The clean peaks have one dominant scattering pattern. With the supplied errors, additional components are not robustly resolved. Scaling errors to smaller pre-peak empirical fluctuations makes weak clean-region second modes detectable, illustrating error-model sensitivity rather than establishing oligomers. WT has a clear additional pattern when the earlier shoulder is included, but its recovered low-q shape and component-window dependence prevent a reliable oligomer/Rg assignment.

The selected-window plot values differ slightly from the broader clean-region EFA profile values because these are direct arithmetic averages over optimized windows, with shared-background uncertainty retained. Guinier and P(r) Rg values are intentionally reported separately. Nearby q/window sensitivity is small for BIFT, but background choice remains a material systematic, particularly for WT.

Read `poster_results/README.md` for the figure caption and exact parameters. The full EFA derivation is in `EFA_equations.md`.''')
code('''print((ROOT/'results/noise_scale_sensitivity.json').read_text())
print((ROOT/'results/provenance.json').read_text())
''')
n=nb.v4.new_notebook(cells=c,metadata={'kernelspec':{'name':'python3','display_name':'SAXS Python','language':'python'},'language_info':{'name':'python','version':'3.12'}})
nb.write(n,ROOT/'02_windows_profiles_and_EFA.ipynb')

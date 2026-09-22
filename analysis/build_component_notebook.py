from pathlib import Path
import nbformat as nb
R=Path(__file__).resolve().parents[1];c=[]
def md(s):c.append(nb.v4.new_markdown_cell(s))
def code(s):c.append(nb.v4.new_code_cell(s))
md(r'''# Actual EFA component I(q): WT and AA
This notebook generates the **component-level scattering curves**, rather than presenting a selected-frame average as an EFA result. Start with the two `*_EFA_candidate_main.dat` files and `component_Iq_comparison.png`.

The question is whether contaminating larger material could inflate AA's apparent size. We do not assume that an extra fitted component is a clean oligomer, or that the main fitted component has been proven monomeric. A two-component calculation is tested as a hypothesis; it is not a selected biological species count.

The existing notebooks 01–03 remain prior-method/reproduction references. This notebook uses a independently assessed residual-background hypothesis and does **not** redo the unavailable original facility solvent subtraction.''')
code('''from pathlib import Path
import sys, json, numpy as np
from IPython.display import display, Image, Markdown
ROOT=Path.cwd();sys.path.insert(0,str(ROOT/'analysis'))
from component_recovery import *
from audit_component_recovery import audit, validate
from component_anchor_sensitivity import anchors_audit
from component_bift_controls import controls
from plot_component_recovery import make_plots
''')
md(r'''## Baseline assumptions and uncertainty
The working baseline is a q-dependent linear interpolation between residual-background anchor means. WT anchors are **55–94 and 210–234**; AA anchors **35–74 and 225–254**. These are independently assessed intervals, not a claim of known analyte-free blanks. WT's post region is moved beyond the EFA fit interval to avoid fitting its tail directly. AA 230–259 was rejected because its mid/high-q drift failed the stationarity screen. All bounds are zero indexed and inclusive.

We compare pre-only and post-only corrections as stress tests, and shift anchors by ±5 frames, rejecting intervals that fail the same stationarity screen. The interpolated model is a working hypothesis, **not a uniquely established optimal baseline**. The provisional custom integral correction is not used here.

Write the corrected matrix as $D=YA$, where each column of $A$ selects a frame and subtracts the appropriate weighted anchor averages. For a component extraction vector $l_k$, $s_k=YAl_k$, so its fixed-elution variance is
$$\mathrm{Var}(s_{ik}\mid C)=\sum_j \sigma_{ij}^2(Al_k)_j^2.$$
This includes the shared residual-background covariance, including overlapping frame/anchor effects. Unknown correlations introduced by the facility's original subtraction remain unavailable.''')
code('''checks, anchor_results=anchors_audit()
assert all(r['passes'] for r in checks)
for r in checks: print(r)
print('Synthetic recovery / covariance checks:', validate())
(ROOT/'component_results/validation.json').write_text(json.dumps(validate(),indent=2));
''')
md(r'''## Forward/backward EFA and actual curve recovery
The row-scaled matrix is $D_w=U\Sigma V^T$. Forward and backward SVD diagnose evolving signal rank; the support intervals then constrain a nonnegative elution rotation. For the resulting elution matrix $C$, the scattering curves are
$$S=D(C^T)^+.$$
We retain the original signed scattering values; **no smoothing or positivity clipping is applied to the recovered I(q)**. Elution coefficients are constrained nonnegative. Their normalization is arbitrary and is not a mole fraction.

Fiducial WT: frames 120–205, supports 120–175 and 140–205. Fiducial AA: frames 145–215, supports 145–185 and 160–215. These are explicit trial supports, not unique intervals established by EFA. The exported main candidate is the component with the second support, scaled to its contribution to the selected frame average. Other supports and q ranges are tested below.''')
code('''import itertools
rows=[]
for name,spec in SPECS.items():
    for kind,end,start,qmax in itertools.product(['linear','pre','post'],spec['ends'],spec['starts'],[.15,.25]):
        r,_=evaluate(name,kind,end,start,qmax);rows.append(r)
    print(name, '72 baseline/support/q-range combinations evaluated')
(OUT/'rotation_grid.json').write_text(json.dumps(rows,indent=2));
''')
md(r'''## Refit the decomposition under noise; do not treat C as exact
Each sample gets 100 parametric perturbations of all 450 original frame profiles, followed by the same residual-background transform and a fresh EFA rotation. The DAT uncertainty column is the standard deviation of these recovered curves. Bootstrap intervals remain conditional on this baseline, component count and supports.

A separate 200-replicate rank-one null simulation tests whether the second singular value exceeds what the supplied errors and the new shared baseline could produce. It tests an additional signal pattern, **not oligomer identity**. Error calibration and unmodelled background structure can affect this test.''')
code('''results={name:audit(name,bootstrap=100,null_reps=200) for name in SPECS}
for name,r in results.items():
    print(name, 'average / candidate Rg:',r['average_guinier']['Rg'],r['candidate_guinier']['Rg'])
    print('Candidate bootstrap 95%:',r['bootstrap_Rg_95'])
    print('Paired candidate-minus-average Rg 95%:',r['bootstrap_delta_Rg_vs_average_95'])
    print('s2 / null 99%:',r['singular_values'][1],r['null_s2_99'])
''')
md(r'''## Separate P(r) curve changes from regularization/error changes
BIFT estimates P(r) for the candidate and the frame average. Decomposition changes both I(q) and its uncertainties, which can change regularization and Dmax even without a corresponding structural change. Therefore we also fit the **unchanged average using candidate errors**, and the candidate using average errors. The latter is a diagnostic control, not its proper uncertainty model.

BIFT here uses diagonal uncertainties and cannot account for the full q covariance induced by decomposition. Its reported errors do not describe all model ambiguity. In particular, a shortened AA P(r) tail cannot automatically be interpreted as removal of an oligomer.''')
code('''controls()
make_plots()
print((OUT/'bift_weight_controls.json').read_text())
display(Image(filename=str(OUT/'component_Iq_comparison.png')))
''')
code('''display(Image(filename=str(OUT/'component_Guinier.png')))
display(Image(filename=str(OUT/'decomposition_diagnostics.png')))
''')
md(r'''## Stability and the contamination question
The plotted sensitivity set includes converged rotations with condition number <30, no main-curve q bins below −3 conditional errors, a valid common-range Guinier slope, and residual mean square no worse than 1.1 times rank one. The candidate main peak must lie within 5 frames of the observed main peak; for WT, the other component must peak before frame 155 to represent the observed shoulder. These screens do **not** certify biological purity. AA's extra fitted curve remains largely noise/negative even when the main curve passes.

All Guinier comparisons below use exactly q=0.012–0.030 Å⁻¹. Differences from earlier reported values reflect both the changed baseline and the fit convention; compare within this notebook, not across methods without the audit.''')
code('''display(Image(filename=str(OUT/'component_model_sensitivity.png')))
print((OUT/'summary.json').read_text())
print('Leading / central / trailing direct-average check:')
for r in json.loads((OUT/'elution_window_check.json').read_text()):
    print(r['sample'],r['frames'],round(r['Rg'],2),'±',round(r['se'],2),'Å; qmaxRg',round(r['qRgmax'],3))
''')
md(r'''## What has and has not been obtained
- **WT:** an actual main-peak component curve with a distinguishable earlier pattern. Its low-q shape changes relative to the average, but its absolute size depends materially on the residual-background model. Neither component is assigned a molecular identity from this alone.
- **AA:** an actual candidate curve from the forced two-component hypothesis. Its second singular value is not resolved above the supplied-error null; the candidate-minus-average Guinier shift includes zero in the conditional bootstrap. Its decomposition is not established as successful removal of oligomer contamination.
- Across AA's peak, direct-average fixed-range Rg remains roughly 31–32 Å, with increasing uncertainty in the tail. This does not prove absence of co-eluting contaminants.
- **No validated purified monomer I(q) is claimed.** The requested component curves, removed-pattern curves, residuals, uncertainties and alternatives are now explicit. Additional raw buffer/sample data and/or independent mass evidence are needed to settle the remaining interpretation.

The component DAT files contain q, intensity and conditional bootstrap sigma. Their arbitrary intensity scale is the component's contribution to the selected mean, not absolute concentration. A separate rank-one curve is exported and explicitly labelled **not purified**. The paired-bootstrap comparison and baseline/support spread are more informative than a single smooth-looking curve.''')
n=nb.v4.new_notebook(cells=c,metadata={'kernelspec':{'name':'python3','display_name':'Python 3','language':'python'}});nb.write(n,R/'04_EFA_component_Iq.ipynb')

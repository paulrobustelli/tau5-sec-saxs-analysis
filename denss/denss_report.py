from pathlib import Path
import json,re,shutil,zipfile,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--runs',type=Path,default=Path.cwd());parser.add_argument('--output',type=Path,default=Path.cwd());args=parser.parse_args();RUN=args.runs.resolve();OUT=args.output.resolve();stats=json.loads((OUT/'summary.json').read_text());cols={'WT':'navy','AA':'orangered'}
fig,axs=plt.subplots(3,2,figsize=(11,9),gridspec_kw={'height_ratios':[3,1,1]},sharex='col')
for c,s in enumerate(cols):
 a=np.loadtxt(OUT/(s+'_average_fit.mrc2sas.fit'));raw=np.loadtxt(OUT/(s+'.dat'));sf=np.median(a[:,1]/np.interp(a[:,0],raw[:,0],raw[:,1]));a[:,1:]/=sf;q,y,e,p=a.T
 indiv=sorted((RUN/(s+'_denss')).glob(s+'_denss_*_map.fit'));chis=[]
 for path in indiv:
  b=np.loadtxt(path);chi=float(np.mean(((b[:,1]-b[:,3])/b[:,2])**2));chis.append(chi);axs[0,c].plot(b[:,0],b[:,3],color=cols[s],lw=.5,alpha=.2)
 chi=float(np.mean(((y-p)/e)**2));axs[0,c].errorbar(q,y,e,fmt='.',ms=3,color='black',elinewidth=.4,label='Measured');axs[0,c].plot(q,p,color=cols[s],label=f'Aligned average (mean χ²={chi:.2f})');axs[0,c].set(yscale='log',ylabel='I(q)',title=s);axs[0,c].legend(fontsize=8)
 axs[1,c].plot(q,(y-p)/e,'.',color=cols[s]);axs[1,c].axhline(0,color='gray');axs[1,c].set(xlabel='q (Å⁻¹)',ylabel='Residual / σ')
 stats[s]['individual_measured_mean_chi2']=chis;stats[s]['average_measured_mean_chi2']=chi
 b=np.loadtxt(OUT/(s+'_refined_fit.mrc2sas.fit'));sf=np.median(b[:,1]/np.interp(b[:,0],raw[:,0],raw[:,1]));b[:,1:]/=sf;rc=float(np.mean(((b[:,1]-b[:,3])/b[:,2])**2));stats[s]['refined_measured_mean_chi2']=rc
 axs[0,c].plot(b[:,0],b[:,3],'--',color='limegreen',lw=1,label=f'Refined (mean χ²={rc:.2f})');axs[0,c].legend(fontsize=8)
 axs[2,c].plot(b[:,0],(b[:,1]-b[:,3])/b[:,2],'.',ms=2,color='green');axs[2,c].axhline(0,color='gray');axs[2,c].set(xlabel='q (Å⁻¹)',ylabel='Refined / σ');axs[1,c].set(ylabel='Average / σ',xlabel='')
 assert len(chis)==20
fig.suptitle('Back-calculated scattering: individual reconstructions (faint) and aligned average');fig.tight_layout();fig.savefig(OUT/'DENSS_fit_checks.png',dpi=180);fig.savefig(OUT/'DENSS_fit_checks.pdf');plt.close(fig)
(OUT/'summary.json').write_text(json.dumps(stats,indent=2))
methods='''# DENSS exploratory envelopes

DENSS 1.8.7, standard SLOW mode, 64³ grid, oversampling 3, 20 independent random starts per sample, 4 workers; default positivity, connectivity and shrinkwrap; no imposed symmetry or atomic reference. Input is the final buffer-subtracted frame average used in Paper Figures, retaining measured errors and q=0.012–0.25 Å⁻¹. WT: frames162–178, buffer220–236; AA: frames168–188, buffer24–70 (zero-based inclusive). Dmax is fixed to the existing BIFT values: WT91.7258Å, AA166.7131Å. DENSS internally regularizes/interpolates the input using Sasrec. Original inputs, fitted profiles and reconstruction logs are included.

These are **effective scattering envelopes for disordered proteins**, not unique molecular structures, atomic models, or recovered conformational ensembles. In general the ensemble average of squared scattering amplitudes is not the squared amplitude of an ensemble-average density. Averaging here is over independent reconstructions, not protein conformations. Detailed lobes, cavities and apparent topology are not established. Dmax, baseline and q-range sensitivity are not tested by repeated random starts. This is a low-q exploratory reconstruction, not a full-q DENSS analysis; DENSS guidance recommends using the full reliable scattering range for final reconstruction. AA's long P(r) tail is sensitive to q-range, so its envelope is conditional on the primary P(r) choice.

Surfaces enclose90%of positive density. Both aligned-average and subsequently refined maps are displayed, with equal spatial axes and no extra display smoothing. This is an explicit visualization threshold, not a molecular boundary or confidence surface. Density uses DENSS's arbitrary default electron normalization; these maps do not estimate molecular weight. DENSS's FSC diagnostic compares reconstructions to a common reference; it is not an independent experimental resolution claim.

Post-averaging refinement was run from each aligned-average map against the same input, with seed20260923 and default DENSS refinement settings. It restores scattering agreement by modifying the density; the refined result is a single effective map, not an average over conformations.

The fit-check plot directly compares back-calculated average-map scattering with the measured data, fitting scale only, no offset. Individual fits are also shown. Averaging can worsen agreement, so the averaged map should not be called a validated single structure even if individual reconstructions fit well.

## Reproduction
Install requirements.txt in an isolated environment, then run `python run_reconstructions.py`. Then run `python render_denss.py` , `python render_denss.py --refined`, and `python denss_report.py` in the same directory to regenerate figures and the notebook. This repeats20random starts; exact seeds from this run are recorded in the archived logs. The complete-run archive includes all individual and aligned maps. Open WT_avg.mrc and AA_avg.mrc in ChimeraX/PyMOL or use the self-contained interactive HTML.

References: [DENSS documentation](https://github.com/tdgrant1/denss), [Grant2018](https://doi.org/10.1038/nmeth.4581).
'''
(OUT/'README.md').write_text(methods)
n=nbformat.v4.new_notebook();n.metadata.kernelspec={'display_name':'Python 3','language':'python','name':'python3'}
n.cells=[nbformat.v4.new_markdown_cell(methods),nbformat.v4.new_code_cell("from pathlib import Path\nimport json\nfrom IPython.display import display, Image, Markdown\nROOT=Path.cwd()\nsummary=json.loads((ROOT/'summary.json').read_text())\ndisplay(Image(filename=str(ROOT/'DENSS_refined_envelopes.png')))\ndisplay(Image(filename=str(ROOT/'DENSS_envelopes.png')))"),nbformat.v4.new_markdown_cell('## Agreement with measured scattering\nFaint curves show all20individual reconstructions; the bold curve is the aligned average. Separate residual panels show the average (sample color) and refined map (green). The reported mean χ² is the mean squared normalized residual, not adjusted for the many effective reconstruction parameters.'),nbformat.v4.new_code_cell("display(Image(filename=str(ROOT/'DENSS_fit_checks.png')))\nfor s,r in summary.items():\n    x=r['individual_measured_mean_chi2']\n    display(Markdown(f\"**{s}**: individual mean χ² range {min(x):.2f}–{max(x):.2f}; aligned-average mean χ² {r['average_measured_mean_chi2']:.2f}; refined mean χ² {r['refined_measured_mean_chi2']:.2f}.\"))"),nbformat.v4.new_markdown_cell('## Interactive maps\nOpen [refined maps](DENSS_refined_interactive.html) or [aligned averages](DENSS_interactive.html) to rotate the two envelopes. Refined maps: [WT](WT_refined.mrc), [AA](AA_refined.mrc). Aligned averages: [WT](WT_avg.mrc), [AA](AA_avg.mrc).'),nbformat.v4.new_code_cell("for s,r in summary.items():\n    print(s)\n    print(r['alignment_log'])")]
nbformat.write(n,OUT/'DENSS Envelopes.ipynb')
with zipfile.ZipFile(OUT/'DENSS_complete_runs.zip','w',zipfile.ZIP_DEFLATED) as z:
 for s in cols:
  for f in (RUN/(s+'_denss')).glob('*'):
   if f.is_file() and 'current' not in f.name:z.write(f,f.relative_to(RUN))

print({s:{'individual_chi2_range':[min(r['individual_measured_mean_chi2']),max(r['individual_measured_mean_chi2'])],'average_chi2':r['average_measured_mean_chi2']} for s,r in stats.items()})

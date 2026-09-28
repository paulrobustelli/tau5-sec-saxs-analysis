"""Absolute errors of ensemble-average shifts, before and after selection."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/Tau5_joint_refinement'
OUT=BASE/'CS_absolute_errors';OUT.mkdir(exist_ok=True)
data=json.loads((ROOT/'outputs/publish_snapshot/Chemical_shift_SAXS/snapshot_data.json').read_text())
atoms={'CA':'Cα','C':'C′','CB':'Cβ','N':'N','NH':'HN'}
exports=[];summary=[]
for atom,name in atoms.items():
    fig,axes=plt.subplots(2,2,figsize=(16,7),sharex=True,sharey=True)
    for row,s in enumerate(('WT','AA')):
        for col,e in enumerate(('1','2')):
            ax=axes[row,col]
            values={}
            for stage in ('raw','all_shifts'):
                rr=[r for r in data['shifts'][s] if r['atom']==atom and r['ensemble']==e and r['fit']==stage]
                values[stage]={int(r['AR_residue']):float(r['experimental'])-float(r['predicted']) for r in rr}
            assert values['raw'].keys()==values['all_shifts'].keys()
            x=np.array(sorted(values['raw']))
            if len(x):
                for stage,offset,color,label in [('raw',-.2,'#aaaaaa','Pool'),('all_shifts',.2,'#2466a8','Selected weights')]:
                    residual=np.array([values[stage][i] for i in x]);absolute=abs(residual)
                    rms=float(np.sqrt(np.mean(residual**2)));mae=float(np.mean(absolute))
                    ax.bar(x+offset,absolute,width=.38,color=color,label=f'{label}: RMSD {rms:.3f} ppm')
                    summary.append([s,e,atom,stage,len(x),rms,mae])
                    exports.extend([s,e,atom,stage,int(i),float(v)] for i,v in zip(x,absolute))
                ax.legend(fontsize=8,loc='upper right')
            else:ax.text(.5,.5,'No measured shifts',transform=ax.transAxes,ha='center')
            ax.axvspan(391,414,color='.94',zorder=0)
            ax.set(title=f'{s} — ensemble {e}',xlim=(329,447),ylabel='|Experimental − calculated| (ppm)')
            ax.set_ylim(bottom=0);ax.spines[['top','right']].set_visible(False)
            if row==1:ax.set_xlabel('AR residue')
    fig.suptitle(f'{name}: unweighted pool vs selected ensemble weights\nAbsolute error of the ensemble-average shift; direct chemical-shift fits, no SAXS restraint. R2 shaded.',fontsize=12)
    fig.tight_layout(rect=[0,0,1,.92])
    for ext in ('png','pdf'):fig.savefig(OUT/f'{atom}_absolute_errors.{ext}',dpi=160)
    plt.close(fig)
for name,header,rows in [('absolute_errors.csv',['protein','ensemble','atom','stage','AR_residue','absolute_error_ppm'],exports),('summary.csv',['protein','ensemble','atom','stage','n_shifts','RMSD_ppm','MAE_ppm'],summary)]:
    with (OUT/name).open('w') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)
print(OUT)

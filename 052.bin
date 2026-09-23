from pathlib import Path
import json,re,shutil,subprocess,os,sys
import numpy as np
import denss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--refined',action='store_true');parser.add_argument('--runs',type=Path,default=Path.cwd());parser.add_argument('--output',type=Path,default=Path.cwd());args=parser.parse_args();RUN=args.runs.resolve();OUT=args.output.resolve();OUT.mkdir(exist_ok=True,parents=True)
COL={'WT':'navy','AA':'orangered'}
fig=plt.figure(figsize=(11,9));interactive=make_subplots(rows=1,cols=2,specs=[[{'type':'scene'},{'type':'scene'}]],subplot_titles=['WT','AA']);stats={}
for c,s in enumerate(['WT','AA']):
 d=RUN/(s+'_denss');prefix=d/(s+'_denss');mapfile=(RUN/(s+'_refined.mrc')) if args.refined else Path(str(prefix)+'_avg.mrc');rho,side=denss.read_mrc(str(mapfile));positive=np.maximum(rho,0);ordered=np.sort(positive.ravel())[::-1];level=ordered[np.searchsorted(np.cumsum(ordered),.90*ordered.sum())];dx=float(side)/rho.shape[0]
 vertices,faces,_,_=marching_cubes(positive,level=level,spacing=(dx,dx,dx));vertices-=float(side)/2
 # Same spatial scale for both samples and both views. No display smoothing.
 for row,az in enumerate([-60,30]):
  ax=fig.add_subplot(2,2,2*row+c+1,projection='3d');mesh=Poly3DCollection(vertices[faces],alpha=.85,facecolor=COL[s],edgecolor='none',shade=True);ax.add_collection3d(mesh);ax.set(xlim=(-95,95),ylim=(-95,95),zlim=(-95,95),xlabel='x (Å)',ylabel='y (Å)',zlabel='z (Å)',title=f'{s} — '+('view 1' if row==0 else 'view 2'));ax.view_init(elev=20,azim=az);ax.set_box_aspect((1,1,1))
 interactive.add_trace(go.Mesh3d(x=vertices[:,0],y=vertices[:,1],z=vertices[:,2],i=faces[:,0],j=faces[:,1],k=faces[:,2],color=COL[s],opacity=.85,name=s),row=1,col=c+1)
 for suffix in ['_avg.mrc','_reference.mrc','_fsc.dat','_allfscs.dat','_map.fit','_final.log']:
  shutil.copy2(str(prefix)+suffix,OUT/(s+suffix))
 shutil.copy2(RUN/(s+'.dat'),OUT/(s+'.dat'));shutil.copy2(RUN/(s+'_denss_fitdata.fit'),OUT/(s+'_fitdata.fit'))
 # Evaluate the averaged map against the measured curve; scale only, no offset.
 command=[sys.executable,'-m','denss.scripts.denss_mrc2sas','-f',str(mapfile),'-d',str(OUT/(s+'.dat')),'--plot_off','-o',str(OUT/(s+('_refined_fit' if args.refined else '_average_fit')))]
 result=subprocess.run(command,capture_output=True,text=True,check=True);(OUT/(s+('_refined_fit_stdout.txt' if args.refined else '_average_fit_stdout.txt'))).write_text(result.stdout+result.stderr)
 log=(Path(str(prefix)+'_final.log')).read_text();stats[s]={'grid':list(rho.shape),'box_A':float(side),'voxel_A':dx,'surface_density_threshold':float(level),'enclosed_positive_density_fraction':.9,'dmax_input_A':91.7257817 if s=='WT' else 166.7131223,'map_density_Rg_A':float(denss.rho2rg(rho,side=side,dx=dx)),'alignment_log':log}
interactive.update_layout(title='DENSS effective envelopes: 90% positive density; equal spatial scale',height=650)
for sn in ['scene','scene2']:interactive.update_layout(**{sn:dict(xaxis=dict(title='x (Å)',range=[-95,95]),yaxis=dict(title='y (Å)',range=[-95,95]),zaxis=dict(title='z (Å)',range=[-95,95]),aspectmode='cube')})
interactive.write_html(OUT/('DENSS_refined_interactive.html' if args.refined else 'DENSS_interactive.html'),include_plotlyjs=True)
fig.suptitle(('DENSS refined envelopes — initialized from aligned averages' if args.refined else 'DENSS aligned averages — 20 random starts per sample')+'\n90% positive-density surfaces; equal spatial scale; no imposed symmetry',fontsize=12);fig.tight_layout();fig.savefig(OUT/('DENSS_refined_envelopes.png' if args.refined else 'DENSS_envelopes.png'),dpi=200);fig.savefig(OUT/('DENSS_refined_envelopes.pdf' if args.refined else 'DENSS_envelopes.pdf'));plt.close(fig)
(OUT/('summary_refined.json' if args.refined else 'summary.json')).write_text(json.dumps(stats,indent=2));print('rendered')

if args.refined:
 for s in ['WT','AA']:
  for f in RUN.glob(s+'_refined*'):
   if 'current' not in f.name:shutil.copy2(f,OUT/f.name)

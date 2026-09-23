"""Authentic PyMOL rendering of DENSS maps; inspired by tdgrant.com/tips/.
Run with a working PyMOL Python installation. No density smoothing or refitting.
"""
from pathlib import Path
import json,struct,numpy as np
import pymol2
OUT=Path(__file__).resolve().parent;SRC=OUT.parent
thresholds={}
for sample in ['WT','AA']:
 b=(SRC/(sample+'_refined.mrc')).read_bytes();nx,ny,nz,mode=struct.unpack_from('<4i',b);extra=struct.unpack_from('<i',b,92)[0];assert mode==2
 a=np.frombuffer(b,dtype='<f4',count=nx*ny*nz,offset=1024+extra);a=np.sort(np.maximum(a,0))[::-1];thresholds[sample]={str(f):float(a[np.searchsorted(np.cumsum(a),f*a.sum())]) for f in [.90,.60,.30]}
with pymol2.PyMOL() as p:
 c=p.cmd
 c.set('normalize_ccp4_maps',0);c.set('orthoscopic',1);c.set('ray_opaque_background',1);c.bg_color('white');c.set('antialias',2);c.set('ray_trace_mode',0);c.set('ray_shadows',0);c.set('ambient',.35);c.set('specular',.35);c.set('transparency_mode',1);c.set('two_sided_lighting',1)
 c.set_color('wt_navy',[.05,.08,.42]);c.set_color('aa_orange',[.9,.23,.04])
 for sample in ['WT','AA']:
  c.load(str(SRC/(sample+'_refined.mrc')),sample+'_density');c.map_double(sample+'_density');c.disable(sample+'_density')
  for fraction,color,trans in [(.90,'marine',.90),(.60,'cyan',.78),(.30,'red',.0)]:
   name=sample+'_'+str(int(fraction*100));c.isosurface(name,sample+'_density',thresholds[sample][str(fraction)]);c.color(color,name);c.set('transparency',trans,name)
  c.group(sample+'_shells',sample+'_90 '+sample+'_60 '+sample+'_30')
  # DENSS tutorial's continuous volume workflow, with thresholds adapted to map scale.
  t=thresholds[sample];ramp=[t['0.9'],0,0,1,0.015,t['0.6'],0,1,1,.09,t['0.3'],1,.5,0,.3,float(t['0.3']*1.5),1,0,0,.45]
  c.volume_ramp_new(sample+'_ramp',ramp);c.volume(sample+'_volume',sample+'_density',sample+'_ramp');c.disable(sample+'_volume')
 # Render each at the exact same camera scale; molecular coordinates stay unchanged.
 for sample in ['WT','AA']:
  c.disable('all');c.enable(sample+'_shells');c.enable(sample+'_90');c.enable(sample+'_60');c.enable(sample+'_30');c.reset();c.turn('x',-15);c.turn('y',25)
  # Set camera position and clipping equally for both maps.
  view=list(c.get_view());view[9:12]=[0,0,-450];view[12:15]=[0,0,0];view[15:18]=[250,650,20];c.set_view(view)
  for i in range(2):
   if i:c.turn('y',90)
   c.png(str(OUT/(sample+'_density_view'+str(i+1)+'.png')),width=1400,height=1400,dpi=300,ray=1)
  c.scene(sample+'_density','store')
 # Side-by-side objects in a saved native PyMOL session; translated for display only.
 c.disable('all')
 for sample,shift in [('WT',-100),('AA',100)]:
  for obj in [sample+'_density',sample+'_90',sample+'_60',sample+'_30',sample+'_volume']:c.translate([shift,0,0],object=obj,camera=0)
  c.enable(sample+'_shells');c.enable(sample+'_90');c.enable(sample+'_60');c.enable(sample+'_30')
 c.reset();view=list(c.get_view());view[9:12]=[0,0,-850];view[12:15]=[0,0,0];view[15:18]=[600,1100,20];c.set_view(view);c.scene('raytraced_shells','store');c.save(str(OUT/'DENSS_density_shells.pze'))
 c.disable('WT_shells');c.disable('AA_shells');c.enable('WT_volume');c.enable('AA_volume');c.scene('DENSS_volume_ramps','store');c.save(str(OUT/'DENSS_volume_ramps.pze'))
 print('PyMOL version:',c.get_version())
(OUT/'thresholds.json').write_text(json.dumps(thresholds,indent=2))
print('Rendered in PyMOL; shell fractions are enclosed positive density, not uncertainty.')

"""Run inside the PyMOL GUI. Actual PyMOL continuous volume rendering, not ray-traced shells."""
from pymol import cmd
from pathlib import Path
import numpy as np,json,struct
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'rainbow'
cmd.reinitialize();cmd.set('normalize_ccp4_maps',0);cmd.bg_color('white');cmd.set('orthoscopic',1);cmd.set('antialias',2);cmd.set('ray_opaque_background',1);cmd.set('ray_shadows',0);cmd.set('ambient',.4);cmd.set('specular',.3)
settings={}
for s in ['WT','AA']:
 fn=ROOT/(s+'_refined.mrc');data=fn.read_bytes();nx,ny,nz,mode=struct.unpack_from('<4i',data);skip=struct.unpack_from('<i',data,92)[0];a=np.frombuffer(data,dtype='<f4',offset=1024+skip,count=nx*ny*nz);ordered=np.sort(np.maximum(a,0))[::-1];thr=float(ordered[np.searchsorted(np.cumsum(ordered),.9*ordered.sum())]);peak=float(a.max())
 cmd.load(str(fn),s+'_map');cmd.map_double(s+'_map')
 cmd.isosurface(s+'_blue',s+'_map',thr);cmd.color('blue',s+'_blue');cmd.set('transparency',.55,s+'_blue')
 ramp=[]
 for frac,col,opacity in [(0.015,'blue',0),(.04,'blue',.025),(.14,'cyan',.025),(.30,'green',.025),(.50,'yellow',.025),(.70,'orange',.025),(.9,'red',.035),(1.1,'red',.04)]:ramp.extend([peak*frac,col,opacity])
 cmd.volume_ramp_new(s+'_rainbow_ramp',ramp);cmd.volume(s+'_rainbow',s+'_map',s+'_rainbow_ramp')
 # Native PyMOL pseudoatoms as explicit density beads, not an atomic protein model.
 # Sample the actual map's 90% contour on a common 4 A grid.
 field=cmd.get_volume_field(s+'_map');extent=np.array(cmd.get_extent(s+'_map'));spacing=(extent[1]-extent[0])/(np.array(field.shape)-1)
 # Trilinear interpolation using numpy only, on equal spacing for both samples.
 grid=np.stack(np.meshgrid(*[np.arange(extent[0,j],extent[1,j],4.) for j in range(3)],indexing='ij'),-1).reshape(-1,3)
 u=(grid-extent[0])/spacing;ix=np.floor(u).astype(int);t=u-ix;valid=np.all((ix>=0)&(ix<np.array(field.shape)-1),axis=1);grid=grid[valid];ix=ix[valid];t=t[valid];vals=np.zeros(len(ix))
 for d0 in [0,1]:
  for d1 in [0,1]:
   for d2 in [0,1]:
    d=np.array([d0,d1,d2]);vals+=field[ix[:,0]+d0,ix[:,1]+d1,ix[:,2]+d2]*np.prod(np.where(d,t,1-t),axis=1)
 coords=grid[vals>=thr]
 lines=[]
 for i,(x,y,z) in enumerate(coords,1):lines.append(f'HETATM{i:5d}  D   DEN A{i%9999:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C')
 cmd.read_pdbstr('\n'.join(lines)+'\nEND',s+'_density_beads');cmd.alter(s+'_density_beads','vdw=2.0');cmd.hide('everything',s+'_density_beads');cmd.show('spheres',s+'_density_beads');cmd.color('gray50',s+'_density_beads');cmd.set('sphere_quality',2)
 settings[s]={'threshold_90':thr,'peak':peak,'ramp':ramp,'bead_spacing_A':4,'bead_radius_A':2,'bead_count':len(coords)}
cmd.disable('all')
def show_denss(sample='WT',style='rainbow',view=1):
 cmd.disable('all');cmd.enable(sample+'_'+style);cmd.reset();cmd.turn('x',-15);cmd.turn('y',25+(90 if int(view)==2 else 0));v=list(cmd.get_view());v[9:12]=[0,0,-450];v[12:15]=[0,0,0];v[15:18]=[250,650,20];cmd.set_view(v);cmd.set('orthoscopic',1);cmd.viewport(1000,1000)
cmd.extend('show_denss',show_denss)
show_denss('WT','rainbow',1)
(OUT/'render_settings.json').write_text(json.dumps(settings,indent=2));print('RAINBOW READY')

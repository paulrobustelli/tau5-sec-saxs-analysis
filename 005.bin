"""CPU orthographic surface rendering. Interpolated normals; no density/mesh smoothing."""
from pathlib import Path
import json,numpy as np
import denss
from skimage.measure import marching_cubes
from PIL import Image,ImageDraw,ImageFont
from matplotlib import font_manager
import plotly.graph_objects as go
from plotly.subplots import make_subplots
OUT=Path(__file__).resolve().parent;SRC=OUT.parent;W,H=2400,1800
colors={'WT':(.06,.10,.42),'AA':(.85,.24,.055)};metadata={};meshes={}
def render(v,f,n,color,angle):
 width,height=1200,900
 theta=np.deg2rad(-60+angle);eye=np.array([np.cos(theta),np.sin(theta),.32]);eye/=np.linalg.norm(eye);right=np.cross([0,0,1],eye);right/=np.linalg.norm(right);up=np.cross(eye,right);rot=np.array([right,up,eye]).T
 p=v@rot;norm=n@rot;screen=np.column_stack((width/2+p[:,0]*height/170,height/2-p[:,1]*height/170,p[:,2]));zbuf=np.full((height,width),-np.inf);rgb=np.ones((height,width,3))
 lights=[(np.array([-.45,.65,1.]),.65),(np.array([.7,.1,.7]),.23)]
 lights=[(x/np.linalg.norm(x),w) for x,w in lights]
 for face in f:
  t=screen[face];lo=np.maximum(np.floor(t[:,:2].min(0)).astype(int),0);hi=np.minimum(np.ceil(t[:,:2].max(0)).astype(int),[width-1,height-1]);
  if np.any(hi<lo):continue
  yy,xx=np.mgrid[lo[1]:hi[1]+1,lo[0]:hi[0]+1];xx=xx+.5;yy=yy+.5
  den=(t[1,1]-t[2,1])*(t[0,0]-t[2,0])+(t[2,0]-t[1,0])*(t[0,1]-t[2,1])
  if abs(den)<1e-10:continue
  a=((t[1,1]-t[2,1])*(xx-t[2,0])+(t[2,0]-t[1,0])*(yy-t[2,1]))/den;b=((t[2,1]-t[0,1])*(xx-t[2,0])+(t[0,0]-t[2,0])*(yy-t[2,1]))/den;c=1-a-b;z=a*t[0,2]+b*t[1,2]+c*t[2,2];region=zbuf[lo[1]:hi[1]+1,lo[0]:hi[0]+1];mask=(a>=0)&(b>=0)&(c>=0)&(z>region)
  if not mask.any():continue
  N=a[...,None]*norm[face[0]]+b[...,None]*norm[face[1]]+c[...,None]*norm[face[2]];N/=np.maximum(np.linalg.norm(N,axis=-1,keepdims=True),1e-9)
  lum=np.full(a.shape,.22);spec=np.zeros(a.shape)
  for light,weight in lights:
   lum+=weight*np.maximum(N@light,0);half=light+[0,0,1];half/=np.linalg.norm(half);spec+=.20*weight*np.maximum(N@half,0)**32
  shade=np.clip(np.array(color)*lum[...,None]+spec[...,None],0,1);region[mask]=z[mask];rr=rgb[lo[1]:hi[1]+1,lo[0]:hi[0]+1];rr[mask]=shade[mask]
 return Image.fromarray(np.uint8(np.clip(rgb,0,1)*255))
im=Image.new('RGB',(W,H),'white')
for col,sample in enumerate(colors):
 rho,side=denss.read_mrc(str(SRC/(sample+'_refined.mrc')));dx=float(side)/rho.shape[0];ordered=np.sort(np.maximum(rho,0).ravel())[::-1];level=ordered[np.searchsorted(np.cumsum(ordered),.9*ordered.sum())]
 v,f,n,_=marching_cubes(rho,level=level,spacing=(dx,dx,dx));v-=float(side)/2
 if np.mean(np.sum(n*(v-v.mean(0)),axis=1))<0:n=-n
 meshes[sample]=(v,f,n)
 metadata[sample]={'threshold':float(level),'density_fraction':.9,'map':sample+'_refined.mrc','smoothing':'lighting normals only; unchanged density and mesh vertices','voxel_A':dx}
 for row,angle in enumerate([0,90]):im.paste(render(v,f,n,colors[sample],angle),(col*1200,row*900))
draw=ImageDraw.Draw(im);fontpath=font_manager.findfont('DejaVu Sans');boldpath=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold'))
def font(size,bold=False):return ImageFont.truetype(boldpath if bold else fontpath,size)
for col,s in enumerate(colors):draw.text((col*1200+600,35),s,font=font(56,True),fill=tuple(int(c*255) for c in colors[s]),anchor='mt')
for row,label in enumerate(['View 1','90° rotation']):draw.text((70,row*900+135),label,font=font(30),fill='#555555')
bar=int(round(50*900/170))
for col in range(2):
 x=col*1200+820;y=H-125;draw.line((x,y,x+bar,y),fill='#333333',width=7);draw.text((x+bar/2,y+16),'50 Å',font=font(30),fill='#333333',anchor='mt')
draw.text((W/2,H-27),'DENSS refined effective envelopes • 90% positive-density contour • identical spatial scale',font=font(25),fill='#555555',anchor='mm')
im.save(OUT/'WT_AA_DENSS_presentation.png',dpi=(240,240));im.save(OUT/'WT_AA_DENSS_presentation.pdf',resolution=240)
(OUT/'render_settings.json').write_text(json.dumps(metadata,indent=2))
p=make_subplots(rows=1,cols=2,specs=[[{'type':'scene'},{'type':'scene'}]],subplot_titles=['WT','AA'])
for col,s in enumerate(colors):
 v,f,_=meshes[s];rgb='rgb('+','.join(str(int(x*255)) for x in colors[s])+')';p.add_trace(go.Mesh3d(x=v[:,0],y=v[:,1],z=v[:,2],i=f[:,0],j=f[:,1],k=f[:,2],color=rgb,flatshading=False,lighting=dict(ambient=.35,diffuse=.8,specular=.35,roughness=.35,fresnel=.1),lightposition=dict(x=200,y=-300,z=450),hoverinfo='skip',name=s),row=1,col=col+1)
scene=dict(xaxis=dict(visible=False,range=[-85,85]),yaxis=dict(visible=False,range=[-85,85]),zaxis=dict(visible=False,range=[-85,85]),aspectmode='cube',bgcolor='white',camera=dict(projection=dict(type='orthographic'),eye=dict(x=1.4,y=-2,z=.7)))
p.update_layout(scene=scene,scene2=scene,paper_bgcolor='white',height=750,margin=dict(l=10,r=10,t=65,b=50),title='DENSS effective envelopes · equal scale · drag to rotate',showlegend=False)
p.write_html(OUT/'DENSS_presentation_interactive.html',include_plotlyjs=True)
print('Rendered unchanged density maps at matching scale')

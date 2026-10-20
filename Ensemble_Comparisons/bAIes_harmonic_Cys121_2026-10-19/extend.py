"""Append empirically sampled terminal Cys geometry; preserve all original core coordinates."""
from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd,mdtraj as md
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
assert not (O/'complete.json').exists()
ref=md.load(str(R/'outputs/WT_MD_300K_2026-10-15/300K_skip5.xtc'),top=str(R/'outputs/WT_MD_300K_2026-10-15/strip.pdb'))
def atom(top,r,n):return next(a.index for a in top.residue(r).atoms if a.name==n)
phiidx,phi=md.compute_phi(ref);omidx,omega=md.compute_omega(ref)
k=next(k for k,a in enumerate(phiidx) if ref.top.atom(a[1]).residue.index==120);ok=np.where((phi[:,k]<-.35)&(phi[:,k]>-3.1))[0]
omk=next(k for k,a in enumerate(omidx) if ref.top.atom(a[-1]).residue.index==120);ok=ok[np.abs(omega[ok,omk])>2.8]
choices=ok[np.linspace(0,len(ok)-1,min(128,len(ok)),dtype=int)]
anchors=[atom(ref.top,119,n) for n in ['CA','C','O']];ci=[a.index for a in ref.top.residue(120).atoms if a.element.symbol!='H'];cnames=[ref.top.atom(i).name for i in ci]
records=[];manifests=[]
for ensemble,source in [('Harmonic',R/'outputs/bAIes_harmonic_backcalc_2026-10-19/manifest.json')]:
 source_records=json.loads(source.read_text())
 for z,m in enumerate(source_records):
  t=md.load(m['file']);keep=[a.index for a in t.top.atoms if a.element.symbol!='H' and not (a.residue.index==119 and a.name=='OXT')];t=t.atom_slice(keep);assert t.n_residues==120
  target=t.xyz[0].astype(float);ai=[atom(t.top,119,n) for n in ['CA','C','O']];Y=target[ai];best=None
  axis=target[ai[1]]-target[ai[0]];axis/=np.linalg.norm(axis);origin=target[ai[1]];oi=atom(t.top,119,'O')
  for angle in [0.,60.,-60.,120.,-120.,180.,30.,-30.,90.,-90.,150.,-150.]:
   ang=np.deg2rad(angle);K=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]]);rotpsi=np.eye(3)*np.cos(ang)+(1-np.cos(ang))*np.outer(axis,axis)+np.sin(ang)*K
   trialO=(target[oi]-origin)@rotpsi.T+origin
   for ti in np.roll(choices,z%len(choices)):
    X=ref.xyz[ti,anchors].astype(float);u,ss,vt=np.linalg.svd((X-X.mean(0)).T@(Y-Y.mean(0)));rot=u@np.diag([1,1,np.linalg.det(u@vt)])@vt;c=(ref.xyz[ti,ci]-X.mean(0))@rot+Y.mean(0);c=(c-origin)@rotpsi.T+origin
    bond=np.linalg.norm(c[cnames.index('N')]-target[atom(t.top,119,'C')])
    if not .125<bond<.145:continue
    td=target.copy();td[oi]=trialO;d=np.linalg.norm(c[:,None]-td[None,:],axis=2)
    for name,old in [('N','C'),('N','CA'),('CA','C')]:d[cnames.index(name),atom(t.top,119,old)]=np.inf
    od=np.linalg.norm(trialO-target,axis=1);od[[oi,ai[0],ai[1]]]=np.inf
    mind=float(min(d.min(),od.min()));score=np.sum(np.maximum(.20-d,0)**2)+np.sum(np.maximum(.20-od,0)**2)
    if best is None or score<best[0]:best=(score,c,int(ti),mind,float(bond),trialO,float(angle))
    if score<1e-12:break
   if best is not None and best[0]<1e-12:break
  assert best is not None,(ensemble,m['id']);score,c,ti,mind,bond,newO,angle=best
  assert mind>=.16,(ensemble,m['id'],mind)
  target[oi]=newO
  top=t.top.copy();res=top.add_residue('CYS',list(top.chains)[0],resSeq=121);new={}
  for name,ix in zip(cnames,ci):new[name]=top.add_atom(name,ref.top.atom(ix).element,res)
  for a,b in [('N','CA'),('CA','CB'),('CB','SG'),('CA','C'),('C','O'),('C','OXT')]:top.add_bond(new[a],new[b])
  top.add_bond(top.atom(atom(t.top,119,'C')),new['N']);xyz=np.concatenate([target,c]);out=md.Trajectory(xyz[None].astype(np.float32),top)
  phidx,ph=md.compute_phi(out);phval=float(ph[0,-1]*180/np.pi);assert -180<phval<-20
  ident=f'HARMONIC121_{ensemble}_{m["variant"]}_{z+1:04d}';dest=O/ensemble/m['variant'];dest.mkdir(parents=True,exist_ok=True);p=dest/(ident+'.pdb');out.save_pdb(str(p));data=p.read_bytes();j=R/'work/maxent_vm_share/jobs'/('COMMON_SAXS_'+ident);j.mkdir(exist_ok=True);old=j/'input.pdb'
  if old.exists():assert old.read_bytes()==data
  else:old.write_bytes(data)
  (j/'ready').touch()
  rec={**m,'id':ident,'source_id':m['id'],'source_file':m['file'],'source_job':m['job'],'ensemble':ensemble,'file':str(p),'job':str(j),'sha256':hashlib.sha256(data).hexdigest(),'predictor_input_sha256':hashlib.sha256(data).hexdigest()};manifests.append(rec)
  records.append({'id':ident,'template_frame':ti,'Cys_phi_deg':phval,'peptide_bond_nm':bond,'min_nonbonded_extension_core_nm':mind,'clash_score_020nm':float(score),'terminal_Pro_psi_rotation_deg':angle,'core_coordinate_max_change_nm':float(np.max(abs(out.xyz[0,:len(target)]-target)))})
  if (z+1)%100==0:print(ensemble,z+1,flush=True)
(O/'manifest.json').write_text(json.dumps(manifests,indent=2));pd.DataFrame(records).to_csv(O/'extension_validation.csv',index=False)
(O/'complete.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'N':len(manifests),'method':'Terminal Cys heavy atoms sampled from supplied121resWT MD; proper rigid alignment to terminalPro CA/C/O, native peptide geometry, negative Cys phi and trans peptide; choose clash-free candidate. Core coordinates retained except terminalPro carbonyl O rotated when needed to choose a clash-free terminalPro psi; terminalPro OXT removed. No new equilibrium sampling or core relaxation.','reference_sha256':hashlib.sha256((R/'outputs/WT_MD_300K_2026-10-15/300K_skip5.xtc').read_bytes()).hexdigest()},indent=2))

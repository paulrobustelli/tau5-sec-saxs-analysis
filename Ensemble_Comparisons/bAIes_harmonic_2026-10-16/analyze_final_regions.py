"""Full production region/Rg analysis; refuses incomplete trajectories."""
from pathlib import Path
import datetime,hashlib,json
import mdtraj as md,numpy as np,pandas as pd
O=Path(__file__).resolve().parent;R=O.parents[1]
regions={'R1':(350,360),'R2_main':(390,404),'R2_minor':(404,410),'R2_historical_core':(394,402),'R3':(430,440)}
seqs=json.loads((R/'outputs/IDPForge_2026-10-13/sequences.json').read_text())
A=O/'final_region_analysis';assert not (A/'complete.json').exists()
for v in ['WT','AA']:
 assert (O/v/'production_final.xml').exists(),v+' not complete'
 with md.formats.XTCTrajectoryFile(str(O/v/'production.xtc')) as f:assert len(f)==20000
A.mkdir(exist_ok=True);summaries=[]
for v in ['WT','AA']:
 topfile=O/v/'protein.pdb';top=md.load(str(topfile)).top;assert top.to_fasta()==[seqs[v]]
 ca=top.select('name CA');rows=[];hs=[];count=0
 for t in md.iterload(str(O/v/'production.xtc'),top=top,chunk=200):
  assert np.isfinite(t.xyz).all() and np.allclose(t.time,np.arange(count+1,count+len(t)+1)*100)
  h=md.compute_dssp(t,simplified=False)=='H';hs.append(h[:,2:]);x=t.xyz[:,ca].astype(float);rg=np.sqrt(np.mean(np.sum((x-x.mean(1,keepdims=True))**2,axis=2),axis=1))*10
  for k in range(len(t)):
   row={'frame_index':count+k,'production_ns':float(t.time[k]/1000),'absolute_ns':float(t.time[k]/1000+20),'Rg_CA_A':float(rg[k])}
   for name,(lo,hi) in regions.items():row[name+'_DSSP_H']=float(h[k,lo-328:hi-327].mean())
   rows.append(row)
  count+=len(t)
 assert count==20000
 df=pd.DataFrame(rows);df.to_csv(A/f'{v}_timeseries.csv',index=False)
 pd.DataFrame({'AR_residue':np.arange(330,448),'DSSP_H':np.concatenate(hs).mean(0)}).to_csv(A/f'{v}_residue_helicity.csv',index=False)
 blocks=df.groupby(df.frame_index//2000).mean(numeric_only=True);blocks.to_csv(A/f'{v}_200ns_blocks.csv')
 summaries.append({'variant':v,'frames':count,'means':df.drop(columns=['frame_index','production_ns','absolute_ns']).mean().to_dict(),'topology_sha256':hashlib.sha256(topfile.read_bytes()).hexdigest(),'note':'Ten chronological 200 ns blocks are descriptive convergence diagnostics, not independent replicate uncertainty.'})
(A/'complete.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'regions':regions,'summary':summaries},indent=2))
print(json.dumps(summaries),flush=True)

from pathlib import Path
import numpy as np,mdtraj as md,openmm as mm,json,subprocess,shutil,os
from openmm import unit
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;PL=R/'work/plumed-reference/src/lib/plumed';rows=[]
for v in ['WT','AA']:
 r=R/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09'/v;dest=O/v;dest.mkdir(exist_ok=True)
 t=md.load(str(r/'production_2us_gpu.xtc'),top=str(r/'protein.pdb'),stride=4000)
 for n in ['baies_params.dat','atom_list.ndx']:shutil.copy2(r/n,dest/n)
 with (dest/'frames.xyz').open('w') as f:
  for frame in t.xyz:
   f.write(f'{len(frame)}\n0 0 0\n');np.savetxt(f,np.c_[np.zeros(len(frame)),frame],fmt=['%d','%.12g','%.12g','%.12g'])
 kt=298.15*8.314462618e-3
 (dest/'plumed.dat').write_text(f'batoms: GROUP NDX_FILE=atom_list.ndx NDX_GROUP=batoms\nbaies: BAIES ATOMS=batoms DATA_FILE=baies_params.dat PRIOR=JEFFREYS TEMP={kt:.15g}\nb: BIASVALUE ARG=baies.ene\nPRINT ARG=baies.ene FILE=COLVAR STRIDE=1 FMT=%0.12f\n')
 subprocess.run([str(PL),'driver','--ixyz','frames.xyz','--length-units','nm','--plumed','plumed.dat','--dump-forces','forces.xyz','--dump-forces-fmt','%0.12f'],cwd=dest,check=True,stdout=(dest/'plumed.log').open('w'),stderr=subprocess.STDOUT)
 pars=np.atleast_2d(np.loadtxt(r/'baies_params.dat'));sys=mm.System()
 for _ in range(t.n_atoms):sys.addParticle(1)
 force=mm.CustomBondForce('-kT*log(erf((r-mu)/(sigma*sqrt(2)))/(r-mu))');force.addGlobalParameter('kT',kt);force.addPerBondParameter('mu');force.addPerBondParameter('sigma')
 for x in pars:force.addBond(int(x[1]-1),int(x[2]-1),[x[3],x[4]])
 sys.addForce(force);ctx=mm.Context(sys,mm.VerletIntegrator(.001),mm.Platform.getPlatformByName('OpenCL' if os.environ.get('GPU_COMPARE') else 'Reference'),{'Precision':'single'} if os.environ.get('GPU_COMPARE') else {})
 es=[];fs=[]
 for xyz in t.xyz:
  ctx.setPositions(xyz);s=ctx.getState(getEnergy=True,getForces=True);es.append(s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole));fs.append(s.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
 pe=np.loadtxt(dest/'COLVAR')[:,1];fl=(dest/'forces.xyz').read_text().splitlines();pf=[]
 for i in range(len(t)):
  block=fl[i*(t.n_atoms+2)+2:(i+1)*(t.n_atoms+2)];pf.append(np.loadtxt(block,usecols=(1,2,3)))
 pf=np.array(pf);fs=np.array(fs)
 rows.append({'variant':v,'frames':len(t),'max_energy_difference_kJ_mol':float(max(abs(np.array(es)-pe))),'max_force_difference_kJ_mol_nm':float(np.max(abs(fs-pf))),'force_RMS_difference':float(np.sqrt(np.mean((fs-pf)**2)))})
 del ctx
np.savez(O/'bias_force_comparison_last.npz',openmm=fs,plumed=pf)
(O/('bias_comparison_gpu.json' if os.environ.get('GPU_COMPARE') else 'bias_comparison.json')).write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))

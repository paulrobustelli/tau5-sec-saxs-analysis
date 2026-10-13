"""Reference PaaA2 LAMMPS versus the exact OpenMM production runner."""
from pathlib import Path
import sys,importlib.util,json,subprocess,shutil,os
import numpy as np
import openmm as mm
from openmm import unit
from openmm.app import Simulation
ROOT=Path(__file__).resolve().parents[2]; OUT=Path(__file__).resolve().parent
case=os.environ.get('AUDIT_CASE','initial'); OUT=OUT/'cases'/(case+('_shifted' if os.environ.get('SHIFT_LJ') else ''));OUT.mkdir(parents=True,exist_ok=True)
UP=ROOT/'work/baies-upstream/tutorial/bAIes/4-simulation'; RUNNER=ROOT/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/code/bAIes_openMM.py'
for name in ['idp_nvt.data','cmap_20240524.cmap']:shutil.copy2(UP/name,OUT/name)
if case!='initial':
 lines=(OUT/'idp_nvt.data').read_text().splitlines();idx=lines.index('Atoms')+2;xyz=np.array([[float(x) for x in line.split()[4:7]] for line in lines[idx:idx+1165]])
 if case=='perturbed':xyz+=np.random.default_rng(73).normal(0,.01,xyz.shape)
 elif case=='expanded':xyz=(xyz-xyz.mean(0))*1.01+xyz.mean(0)
 for k,x in enumerate(xyz):lines[idx+k]=' '.join(lines[idx+k].split()[:4])+ ' '+ ' '.join(f'{z:.12g}' for z in x)
 (OUT/'idp_nvt.data').write_text('\n'.join(lines)+'\n')
source=(UP/'idp_nvt.in').read_text().split('thermo_style')[0]
(OUT/'reference.in').write_text(source+'\nthermo_style custom step pe ebond eangle edihed evdwl f_drycmap\nthermo 1\nthermo_modify format float %.15g\nrun 0\nwrite_dump all custom reference_forces.dat id x y z fx fy fz modify sort id format float %.15g\n')
subprocess.run([str(ROOT/'work/lammps-build/lmp'),'-in','reference.in','-log','reference.log'],cwd=OUT,check=True,stdout=(OUT/'lammps_stdout.log').open('w'))
spec=importlib.util.spec_from_file_location('runner',RUNNER);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
instances=[]
m.HBonds=None
class Capture(Simulation):
 def __init__(self,*a,**k):super().__init__(*a,**k);instances.append(self)
 def minimizeEnergy(self,**k):pass
m.Simulation=Capture
sys.argv=['runner','-ff',str(RUNNER.parent/'coil_202602.xml'),'-cmap',str(RUNNER.parent/'cmap_20240524.cmap'),'-pdb',str(ROOT/'work/baies-upstream/tutorial/bAIes/2-preprocessing/idp.pdb'),'-xtc',str(OUT/'unused.xtc'),'-platform','Reference','-Nsteps','0']
m.main();sim=instances[0];system=sim.system
# LAMMPS data/PDB share atom ordering and coordinates; use dump to verify.
lines=(OUT/'reference_forces.dat').read_text().splitlines();idx=next(i for i,l in enumerate(lines) if l.startswith('ITEM: ATOMS'))
dump=np.loadtxt(lines[idx+1:]);positions=dump[:,1:4]/10;refF=dump[:,4:]*41.84
assert len(positions)==system.getNumParticles()
for i,f in enumerate(system.getForces()):
 if isinstance(f,mm.CustomNonbondedForce) and os.environ.get("SHIFT_LJ"):
  f.setEnergyFunction(f.getEnergyFunction().replace("4*epsilon*((sigma/r)^12-(sigma/r)^6)","4*epsilon*((sigma/r)^12-(sigma/r)^6)+epsilon"))
 f.setForceGroup(i)
context=mm.Context(system,mm.VerletIntegrator(.001),mm.Platform.getPlatformByName('Reference'));context.setPositions(positions)
energies={};forces={}
for i,f in enumerate(system.getForces()):
 state=context.getState(getEnergy=True,getForces=True,groups=1<<i);label=f.__class__.__name__;energies[label]=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);forces[label]=state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
st=context.getState(getEnergy=True,getForces=True);of=st.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer);np.savez(OUT/'force_comparison.npz',reference_forces=refF,openmm_forces=of,positions_nm=positions,**forces)
log=(OUT/'reference.log').read_text().splitlines();i=next(i for i,l in enumerate(log) if l.startswith('Step ') or l.strip().startswith('Step '));values=np.array([float(x) for x in log[i+1].split()]);refE=dict(zip(['total','bond','angle','dihedral','nonbonded','cmap'],values[1:]*4.184))
result={'reference_energy_kJ_mol':refE,'openmm_energy_kJ_mol':energies,'openmm_total_kJ_mol':st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),'force_difference_RMS_kJ_mol_nm':float(np.sqrt(np.mean((of-refF)**2))),'force_difference_max_kJ_mol_nm':float(np.max(abs(of-refF))),'reference_force_RMS_kJ_mol_nm':float(np.sqrt(np.mean(refF**2)))}
(OUT/'engine_comparison.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

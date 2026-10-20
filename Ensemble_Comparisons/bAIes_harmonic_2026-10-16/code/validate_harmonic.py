from pathlib import Path
import ast,json,hashlib
import numpy as np
from openmm import System,CustomBondForce,VerletIntegrator,Context,Platform
from openmm.unit import nanometer,kilojoule_per_mole
O=Path(__file__).resolve().parents[1]
src=O/'code/bAIes_openMM.py'
tree=ast.parse(src.read_text());expr=[n.args[0].value for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='CustomBondForce'][0]
kt=298.15*8.314462618e-3;mu=.6;sigma=.1
system=System();system.addParticle(12);system.addParticle(12)
f=CustomBondForce(expr);f.addGlobalParameter('kT',kt);f.addPerBondParameter('mu');f.addPerBondParameter('sigma');f.addBond(0,1,[mu,sigma]);system.addForce(f)
integrator=VerletIntegrator(.002);ctx=Context(system,integrator,Platform.getPlatformByName('Reference'));rows=[]
for z in [-3,-1,0,1,3]:
 r=mu+z*sigma;ctx.setPositions([[0,0,0],[r,0,0]]*nanometer);state=ctx.getState(getEnergy=True,getForces=True)
 e=state.getPotentialEnergy().value_in_unit(kilojoule_per_mole);force=state.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole/nanometer)[1,0]
 assert np.isclose(e,.5*kt*z*z,atol=1e-10)
 assert np.isclose(force,-kt*z/sigma,atol=1e-10)
 rows.append(dict(z=z,energy_kBT=e/kt,force_kJ_mol_nm=float(force)))
(O/'harmonic_validation.json').write_text(json.dumps(dict(expression=expr,temperature_K=298.15,passed=True,runner_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),checks=rows),indent=2))
print('PASS: actual runner expression energy and force at 0, ±1, ±3 sigma')

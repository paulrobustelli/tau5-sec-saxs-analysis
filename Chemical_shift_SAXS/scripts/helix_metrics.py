"""Explicit DSSP and Ramachandran definitions for ensemble convergence."""
from pathlib import Path
import json
import numpy as np
import mdtraj as md

ROOT=Path(__file__).resolve().parents[1]
DEFINITIONS={
    'DSSP_H':'DSSP alpha helix H only',
    'DSSP_HGI':'DSSP H (alpha), G (3_10), or I (pi)',
    'Rama_alphaR':'-180 <= phi <= 0 and -120 <= psi <= 50 degrees; broad alphaR basin, not necessarily a hydrogen-bonded helix',
    'Rama_alpha_core':'-100 <= phi <= -30 and -80 <= psi <= -5 degrees; narrower helix-core sensitivity definition',
}

def get_assignments(results):
    for r in results:
        if 'phi_deg' in r and 'dssp' in r: continue
        p=ROOT/r['protonated_pdb'];t=md.load(str(p))
        r['dssp']=md.compute_dssp(t,simplified=False)[0].tolist()
        for name,fn,center in [('phi_deg',md.compute_phi,2),('psi_deg',md.compute_psi,1)]:
            indices,angles=fn(t);values=[None]*120
            for atoms,angle in zip(indices,np.rad2deg(angles[0])):values[t.topology.atom(int(atoms[center])).residue.index]=float(angle)
            r[name]=values
        (p.parent/'result.json').write_text(json.dumps(r)+'\n')
    dssp=np.array([r['dssp'] for r in results])
    phi=np.array([r['phi_deg'] for r in results],dtype=float)
    psi=np.array([r['psi_deg'] for r in results],dtype=float)
    valid=np.isfinite(phi)&np.isfinite(psi)
    broad=(phi>=-180)&(phi<=0)&(psi>=-120)&(psi<=50)
    core=(phi>=-100)&(phi<=-30)&(psi>=-80)&(psi<=-5)
    # Missing terminal dihedrals remain NaN, not a false zero occupancy.
    return {'DSSP_H':(dssp=='H').astype(float),'DSSP_HGI':np.isin(dssp,['H','G','I']).astype(float),
            'Rama_alphaR':np.where(valid,broad.astype(float),np.nan),'Rama_alpha_core':np.where(valid,core.astype(float),np.nan)}

def profile(weights,matrix):
    valid=np.isfinite(matrix)
    denominator=weights@valid
    return np.divide(weights@np.nan_to_num(matrix),denominator,out=np.full(matrix.shape[1],np.nan),where=denominator>0)

"""Measure actual dihedral run lengths and DSSP for helix-enriched candidates."""
import argparse
import json
from pathlib import Path
import mdtraj as md
import numpy as np


def verify(folder):
    folder = Path(folder)
    manifest = json.loads((folder/'generation_manifest.json').read_text())
    report = []
    for filename, target in manifest['conformers'].items():
        pdb = folder/filename
        if not pdb.exists():
            report.append(dict(file=filename, status='missing'))
            continue
        t = md.load(str(pdb))
        if t.n_residues != 120:
            raise ValueError(f'{pdb}: incomplete sequence')
        phi = np.full(120, np.nan); psi = phi.copy()
        for values, method, center in [(phi,md.compute_phi,2),(psi,md.compute_psi,1)]:
            ids, angles = method(t)
            for atoms, value in zip(ids, np.rad2deg(angles[0])):
                values[t.topology.atom(int(atoms[center])).residue.index] = value
        core = (phi>=-100)&(phi<=-30)&(psi>=-80)&(psi<=-5)
        start, end = target['start_model']-1, target['end_model']-1
        exact = bool(core[start:end+1].all() and not core[start-1] and not core[end+1])
        dssp = md.compute_dssp(t, simplified=False)[0]
        runs = []
        i = 0
        while i < 120:
            if dssp[i] != 'H':
                i += 1; continue
            j = i+1
            while j < 120 and dssp[j]=='H': j += 1
            if i <= end and j > start:
                runs.append(dict(start_AR=i+328, end_AR=j+327, length=j-i))
            i = j
        report.append(dict(file=filename, target=target, exact_alpha_core_run=exact,
            phi_target_deg=phi[start:end+1].tolist(), psi_target_deg=psi[start:end+1].tolist(),
            DSSP_target=''.join(dssp[start:end+1]), overlapping_DSSP_H_runs=runs))
    out = dict(requested=len(report), available=sum(x.get('status')!='missing' for x in report),
        exact_alpha_core_pass=sum(x.get('exact_alpha_core_run',False) for x in report), results=report,
        note='Requested lengths describe torsion runs. DSSP H requires hydrogen-bond geometry and can differ.')
    (folder/'helix_length_validation.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='results'}))
    if out['exact_alpha_core_pass'] != out['requested']:
        raise RuntimeError('Missing or nonmatching helix targets; inspect validation report')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folder');verify(p.parse_args().folder)

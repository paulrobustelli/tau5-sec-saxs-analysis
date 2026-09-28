"""Explicit, labeled canonical-helix enrichment around D2D >= 10% sites.

Extends IDPConfGen in this process only; upstream installation is unchanged.
Lengths refer to consecutive alpha-core phi/psi residues, not DSSP H counts.
"""
import argparse
import json
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'outputs/Tau5_joint_refinement'
ACTIVE = None


def make_plan(sample, copies):
    data = np.loadtxt(BASE/sample/'d2d_populations.csv', delimiter=',', skiprows=1)
    sequence = json.loads((ROOT/'outputs/ASTEROIDS_setup/manifest.json').read_text())['samples'][sample]['sequence']
    anchors = [int(r[0]) for r in data if r[2] >= .10]
    windows = {}
    for anchor in anchors:
        for length in (4, 6, 8, 10, 12):
            start = max(3, min(anchor-(length-1)//2, 119-length+1))
            windows.setdefault((start, start+length-1), []).append(anchor+327)
    return [dict(start_model=start, end_model=end, start_AR=start+327,
                 end_AR=end+327, length=end-start+1, anchors_AR=sites, copy=copy,
                 proline_positions_AR=[i+327 for i in range(start,end+1) if sequence[i-1]=='P'])
            for copy in range(1, copies+1)
            for (start, end), sites in sorted(windows.items())]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('sample', choices=['WT', 'AA'])
    p.add_argument('--replicate', type=int, choices=[1, 2], required=True)
    p.add_argument('--copies', type=int, default=2)
    p.add_argument('--pilot', action='store_true')
    args = p.parse_args()
    plan = make_plan(args.sample, args.copies)
    if args.pilot:
        plan = [next(x for x in plan if x['length']==n) for n in (4,6,8,10,12)]
    suffix = 'pilot' if args.pilot else str(args.replicate)
    dest = BASE/args.sample/f'helix_lengths_{suffix}'
    dest.mkdir(parents=True, exist_ok=True)
    if list(dest.glob('conformer_*.pdb')):
        raise RuntimeError('Refusing to overwrite existing enrichment structures')
    manifest = dict(sample=args.sample, replicate=args.replicate, cutoff=.10,
        interpretation='Supplementary proposal distribution; not equilibrium populations',
        requested_length_definition='consecutive alpha-core phi/psi; DSSP measured separately',
        torsions='omega=180; phi=-63 +/- 5 and psi=-42 +/- 5 degrees, Gaussian clipped at +/- 12',
        boundaries='immediate flanking psi set to +140 degrees to delimit alpha-core run',
        outside_window='original D2D-informed IDPConfGen fragment sampling',
        conformers={f'conformer_{i}.pdb':v for i,v in enumerate(plan,1)})
    (dest/'generation_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    mp.set_start_method('fork')
    from idpconfgen import cli_build as build
    from idpconfgen.cli import maincli
    original_angles = build.get_adjacent_angles
    original_generator = build.conformer_generator

    def get_angles(*a, **kw):
        original = original_angles(*a, **kw)
        def sample_angles(atom_index):
            template, raw = original(atom_index)
            angles = np.array(raw, copy=True).reshape(-1,3)
            first = build.calc_residue_num_from_index(atom_index)+1
            if ACTIVE is not None:
                for j, row in enumerate(angles):
                    pos = first+j
                    if ACTIVE['start_model'] <= pos <= ACTIVE['end_model']:
                        jitter = np.clip(np.random.normal(0,5,2), -12,12)
                        row[:] = np.deg2rad([180, -63+jitter[0], -42+jitter[1]])
                    elif pos in (ACTIVE['start_model']-1, ACTIVE['end_model']+1):
                        row[2] = np.deg2rad(140)
            return template, angles.ravel()
        return sample_angles

    def generator(**kw):
        global ACTIVE
        original = original_generator(**kw)
        yield next(original)  # atom labels
        for window in plan:
            ACTIVE = window
            yield next(original)
        ACTIVE = None

    build.get_adjacent_angles = get_angles
    build.conformer_generator = generator
    seed = 20269000 + (1000 if args.sample=='AA' else 0) + 100*args.replicate
    sys.argv = ['idpconfgen','build','-db',str(ROOT/'work/software/idpconfgen_database_2024.json'),
        '-seq',str(ROOT/f'outputs/ASTEROIDS_setup/{args.sample}/{args.sample}_120aa.fasta'),
        '-nc',str(len(plan)), '-csss',str(BASE/args.sample/'csss.json'),
        '--dloop-off','-et','pairs','-of',str(dest),'-rs',str(seed),'-n','1']
    maincli()
    if len(list(dest.glob('conformer_*.pdb'))) != len(plan):
        raise RuntimeError('Incomplete enrichment batch; inspect builder logs')


if __name__ == '__main__':
    main()

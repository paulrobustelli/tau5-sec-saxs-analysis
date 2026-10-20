#!/usr/bin/env python3
"""Count cis peptide bonds and plot time series and histogram.

Install: python -m pip install mdtraj numpy pandas matplotlib
Example:
  python plot_cis_bonds.py --trajectory trajectory.xtc --topology protein.pdb \
      --output cis_plots --title 'WT REST 300 K' --ar-start 330 --tag-residues 2

MDTraj reads trajectory times in ps; plots use ns. --ar-start refers to the
first residue AFTER --tag-residues. Omit --ar-start to use topology numbering.
Counts include ALL peptide bonds, including non-proline and tag bonds.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--trajectory', required=True)
    p.add_argument('--topology', required=True)
    p.add_argument('--output', default='cis_plots')
    p.add_argument('--title', default='Cis peptide bonds')
    p.add_argument('--threshold', type=float, default=30.0, help='Cis if |omega| < threshold degrees')
    p.add_argument('--ar-start', type=int, help='First biological residue number after tags; single-chain only')
    p.add_argument('--tag-residues', type=int, default=0)
    p.add_argument('--chunk', type=int, default=1000)
    args = p.parse_args()
    if not 0 < args.threshold < 180:
        p.error('--threshold must lie between 0 and 180 degrees')
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    times, states = [], []
    bond_indices = None
    for trajectory in md.iterload(args.trajectory, top=args.topology, chunk=args.chunk):
        indices, omega = md.compute_omega(trajectory, periodic=False)
        if bond_indices is None:
            bond_indices, topology = indices, trajectory.topology
            if len(indices) == 0:
                raise ValueError('No peptide omega dihedrals found; check topology and bonds.')
            if args.ar_start is not None and topology.n_chains != 1:
                raise ValueError('--ar-start mapping requires a single-chain topology.')
        else:
            assert np.array_equal(indices, bond_indices)
        degrees = np.rad2deg(omega)
        if not np.isfinite(degrees).all():
            raise ValueError('Nonfinite omega dihedrals found.')
        states.append(np.abs(degrees) < args.threshold)
        times.append(trajectory.time / 1000.0)
    if not states:
        raise ValueError('Trajectory contains no frames.')
    cis, time_ns = np.concatenate(states), np.concatenate(times)
    if not np.isfinite(time_ns).all() or np.any(np.diff(time_ns) <= 0):
        raise ValueError('Trajectory times must be finite and strictly increasing.')
    following = [topology.atom(int(ix[2])).residue for ix in bond_indices]
    labels, numbers = [], []
    for residue in following:
        if args.ar_start is None:
            number = residue.resSeq
            label = f'Chain {residue.chain.index}: {residue.name}{number}'
        elif residue.index < args.tag_residues:
            number = None
            label = f'Tag {residue.name} ({residue.index + 1})'
        else:
            number = args.ar_start + residue.index - args.tag_residues
            label = f'{residue.name}{number}'
        labels.append(label)
        numbers.append(number)
    proline = np.array([r.name == 'PRO' for r in following])
    polyproline = np.array([n is not None and 372 <= n <= 379 for n in numbers]) if args.ar_start is not None else np.zeros(len(labels), dtype=bool)
    total = cis.sum(axis=1)
    counts = pd.DataFrame({'time_ns': time_ns, 'all_peptide_bonds_cis': total,
                           'X_Pro_bonds_cis': cis[:, proline].sum(axis=1)})
    if polyproline.any():
        counts['AR372_379_incoming_bonds_cis'] = cis[:, polyproline].sum(axis=1)
    counts.to_csv(out / 'cis_counts_vs_time.csv', index=False)
    pd.DataFrame(cis.astype(int), columns=labels).assign(time_ns=time_ns).to_csv(out / 'cis_states_by_bond.csv', index=False)
    pd.DataFrame({'bond_label': labels, 'following_topology_residue_index': [r.index for r in following],
                  'cis_fraction': cis.mean(axis=0)}).to_csv(out / 'cis_bond_summary.csv', index=False)

    # Time series and bond-by-bond raster. Bonds are labeled by the following residue.
    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True, layout='constrained')
    axs[0].step(time_ns, total, where='mid', color='#333333', lw=0.8, label=f'All {cis.shape[1]} peptide bonds')
    if polyproline.any():
        axs[0].step(time_ns, cis[:, polyproline].sum(axis=1), where='mid', color='#0072b2', lw=0.65, label='Bonds entering AR372–379')
    axs[0].set(title=args.title, ylabel='Cis bonds per frame', yticks=np.arange(total.max() + 1))
    axs[0].legend(); axs[0].grid(axis='y', alpha=0.2)
    active = np.flatnonzero(cis.any(axis=0))
    for row, column in enumerate(active):
        axs[1].scatter(time_ns[cis[:, column]], np.full(cis[:, column].sum(), row), s=3, marker='|', color='#0072b2')
    axs[1].set(yticks=np.arange(len(active)), yticklabels=[labels[i] for i in active],
               xlabel='Saved trajectory time (ns)', ylabel='Residue following cis bond')
    if not len(active):
        axs[1].text(0.5, 0.5, 'No cis bonds detected', transform=axs[1].transAxes, ha='center')
    fig.supxlabel(f'Cis: |ω| < {args.threshold:g}°. Unweighted frames; REST frame changes are not physical kinetics.', fontsize=10)
    for extension in ('png', 'pdf'):
        fig.savefig(out / f'cis_vs_time.{extension}', dpi=200)
    plt.close(fig)

    # Integer-count histogram normalized by number of frames.
    k = np.arange(total.max() + 1)
    frequency = np.bincount(total, minlength=len(k))
    percent = frequency / len(total) * 100
    pd.DataFrame({'cis_bonds_per_frame': k, 'frames': frequency, 'percent': percent}).to_csv(out / 'cis_count_histogram.csv', index=False)
    fig, ax = plt.subplots(figsize=(7.5, 4.8), layout='constrained')
    ax.bar(k, percent, width=0.78, color='#0072b2')
    for count, pct, n in zip(k, percent, frequency):
        ax.text(count, pct + 0.7, f'{pct:.1f}%\n({n})', ha='center', fontsize=10)
    ax.axvline(total.mean(), color='#d55e00', ls='--', label=f'Mean = {total.mean():.2f}')
    ax.set(xticks=k, xlabel='Number of cis peptide bonds per frame', ylabel='Frames (%)',
           title=args.title, ylim=(0, percent.max() + 8))
    ax.legend(frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.supxlabel(f'All {cis.shape[1]} peptide bonds; cis: |ω| < {args.threshold:g}°. {len(total):,} unweighted frames.', fontsize=10)
    for extension in ('png', 'pdf'):
        fig.savefig(out / f'cis_count_histogram.{extension}', dpi=200)
    print(f'{len(total)} frames; {cis.shape[1]} peptide bonds; mean cis count {total.mean():.6f}; outputs: {out}')


if __name__ == '__main__':
    main()

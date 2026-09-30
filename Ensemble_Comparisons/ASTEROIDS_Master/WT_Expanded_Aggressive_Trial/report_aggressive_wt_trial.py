"""Plot and tabulate the WT-expanded aggressive-fit trial."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TRIAL = ROOT / "outputs/Tau5_ASTEROIDS_trials/WT_expanded_aggressive_CS_lowq"
BASELINE = ROOT / "outputs/Tau5_ASTEROIDS/WT/replicate_1/expanded/seed_20260929"
SETUP = ROOT / "outputs/ASTEROIDS_setup/WT"
COLORS = {"Raw pool": "#777777", "Current fit": "#2876b8", "2x low / 0.5x tail": "#df762f", "4x low / 0.25x tail": "#24946d"}


def guinier(q, intensity, sigma=None, cutoff=0.03):
    mask = (q <= cutoff) & (intensity > 0)
    w = np.ones(mask.sum()) if sigma is None else intensity[mask] / sigma[mask]
    slope, intercept = np.polyfit(q[mask] ** 2, np.log(intensity[mask]), 1, w=w)
    return float(np.exp(intercept)), float(np.sqrt(max(0, -3 * slope)))


def normalize_i0(q, intensity, target_i0):
    return intensity * target_i0 / guinier(q, intensity)[0]


def main():
    trial = json.loads((TRIAL / "comparison.json").read_text())
    current_summary = json.loads((BASELINE / "summary.json").read_text())
    raw_shift = pd.read_csv(TRIAL / "raw_shifts.csv")
    current_shift = pd.read_csv(BASELINE / "shifts.csv")
    shift2 = pd.read_csv(TRIAL / "low2_tail0p5/shifts.csv")
    shift4 = pd.read_csv(TRIAL / "low4_tail0p25/shifts.csv")
    datasets = {
        "Raw pool": raw_shift.rename(columns={"raw_secondary_ppm": "calc", "Exp_minus_raw_ppm": "resid"}),
        "Current fit": current_shift.rename(columns={"selected_secondary_ppm": "calc", "Exp_minus_selected_ppm": "resid"}),
        "2x low / 0.5x tail": shift2.rename(columns={"selected_secondary_ppm": "calc", "Exp_minus_selected_ppm": "resid"}),
        "4x low / 0.25x tail": shift4.rename(columns={"selected_secondary_ppm": "calc", "Exp_minus_selected_ppm": "resid"}),
    }
    atoms = [a for a in ["CA", "C", "CB", "N", "NH", "HA"] if a in set(raw_shift.atom)]
    rows = []
    for label, data in datasets.items():
        for atom in atoms:
            r = data.loc[data.atom == atom, "resid"].to_numpy()
            rows.append({"ensemble": label, "observable": atom, "RMSD": np.sqrt(np.mean(r * r)), "unit": "ppm"})
    for label, key in [("Raw pool", "raw"), ("2x low / 0.5x tail", "low2_tail0p5"), ("4x low / 0.25x tail", "low4_tail0p25")]:
        for metric in ["SAXS_chi2_all", "SAXS_chi2_low_q_lt_0p04", "SAXS_chi2_middle_0p04_to_0p10", "SAXS_chi2_tail_ge_0p10", "calculated_Guinier_Rg_A"]:
            rows.append({"ensemble": label, "observable": metric, "RMSD": trial[key][metric], "unit": "chi2/point" if "chi2" in metric else "A"})

    current_saxs = pd.read_csv(BASELINE / "saxs.csv")
    q_current = current_saxs.q_Ainv.to_numpy(); exp_current = current_saxs.I_exp.to_numpy(); err_current = current_saxs.sigma_exp.to_numpy()
    exp_i0_current, _ = guinier(q_current, exp_current, err_current)
    current_i = normalize_i0(q_current, current_saxs.I_selected.to_numpy(), exp_i0_current)
    z_current = (current_i - exp_current) / err_current
    current_saxs_metrics = {
        "SAXS_chi2_all": np.mean(z_current ** 2),
        "SAXS_chi2_low_q_lt_0p04": np.mean(z_current[q_current < .04] ** 2),
        "SAXS_chi2_middle_0p04_to_0p10": np.mean(z_current[(q_current >= .04) & (q_current < .10)] ** 2),
        "SAXS_chi2_tail_ge_0p10": np.mean(z_current[q_current >= .10] ** 2),
        "calculated_Guinier_Rg_A": guinier(q_current, current_i)[1],
    }
    for metric, value in current_saxs_metrics.items():
        rows.append({"ensemble": "Current fit", "observable": metric, "RMSD": value, "unit": "chi2/point" if "chi2" in metric else "A"})
    metrics = pd.DataFrame(rows)
    metrics.to_csv(TRIAL / "comparison_metrics.csv", index=False)

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(len(atoms), 1, figsize=(15, 2.5 * len(atoms)), sharex=True, layout="constrained")
    for ax, atom in zip(np.atleast_1d(axes), atoms):
        exp = raw_shift[raw_shift.atom == atom]
        ax.plot(exp.AR_residue, exp.experimental_secondary_ppm, "k.-", label="Experiment", lw=1.4)
        for label, data in datasets.items():
            d = data[data.atom == atom]
            ax.plot(d.AR_residue, d.calc, color=COLORS[label], label=label, lw=1.25)
        ax.set_ylabel(f"Δδ {atom} (ppm)")
        ax.legend(ncol=5, fontsize=8)
    axes[-1].set_xlabel("AR residue")
    fig.suptitle("WT expanded pool: experimental and calculated secondary shifts")
    fig.savefig(TRIAL / "secondary_shifts_comparison.png", dpi=180)
    fig.savefig(TRIAL / "secondary_shifts_comparison.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    x = np.arange(len(atoms)); width = 0.19
    for i, (label, data) in enumerate(datasets.items()):
        vals = [np.sqrt(np.mean(data.loc[data.atom == atom, "resid"] ** 2)) for atom in atoms]
        ax.bar(x + (i - 1.5) * width, vals, width, color=COLORS[label], label=label)
    ax.set_xticks(x, atoms); ax.set_ylabel("RMSD (ppm)"); ax.set_title("Chemical-shift RMSD by nucleus"); ax.legend()
    fig.savefig(TRIAL / "chemical_shift_RMSD.png", dpi=180)
    fig.savefig(TRIAL / "chemical_shift_RMSD.pdf")
    plt.close(fig)

    hcur = pd.read_csv(BASELINE / "helicity.csv"); h2 = pd.read_csv(TRIAL / "low2_tail0p5/helicity.csv"); h4 = pd.read_csv(TRIAL / "low4_tail0p25/helicity.csv")
    hcur = hcur[hcur.definition == "DSSP_H"]
    d2d_path = ROOT / "outputs/ASTEROIDS_master/data/WT_d2d.csv"
    d2d = pd.read_csv(d2d_path) if d2d_path.exists() else None
    fig, axes = plt.subplots(2, 1, figsize=(15, 7), layout="constrained")
    for ax, limits in zip(axes, [(330, 446), (391, 414)]):
        ax.plot(hcur.AR_residue, hcur.pool, color=COLORS["Raw pool"], label="Raw pool")
        ax.plot(hcur.AR_residue, hcur.selected, color=COLORS["Current fit"], label="Current fit")
        ax.plot(h2.AR_residue, h2.selected_DSSP_H, color=COLORS["2x low / 0.5x tail"], label="2x low / 0.5x tail")
        ax.plot(h4.AR_residue, h4.selected_DSSP_H, color=COLORS["4x low / 0.25x tail"], label="4x low / 0.25x tail")
        if d2d is not None:
            ax.plot(d2d.author_residue, d2d.helix, "k--", lw=1.6, label="δ2D")
        ax.set(xlim=limits, ylim=(0, 1), ylabel="DSSP-H population")
        ax.legend(ncol=5, fontsize=8)
    axes[-1].set_xlabel("AR residue")
    fig.suptitle("WT expanded pool: helicity before and after selection")
    fig.savefig(TRIAL / "helicity_comparison.png", dpi=180)
    fig.savefig(TRIAL / "helicity_comparison.pdf")
    plt.close(fig)

    raw_saxs = pd.read_csv(TRIAL / "raw_saxs.csv"); sx2 = pd.read_csv(TRIAL / "low2_tail0p5/saxs.csv"); sx4 = pd.read_csv(TRIAL / "low4_tail0p25/saxs.csv")
    q = raw_saxs.q_Ainv.to_numpy(); exp = raw_saxs.I_exp.to_numpy(); err = raw_saxs.sigma_exp.to_numpy(); exp_i0, exp_rg = guinier(q, exp, err)
    curves = {
        "Raw pool": raw_saxs.I_raw_I0_normalized.to_numpy(),
        "Current fit": normalize_i0(q, current_saxs.I_selected.to_numpy(), exp_i0),
        "2x low / 0.5x tail": sx2.I_selected_I0_normalized.to_numpy(),
        "4x low / 0.25x tail": sx4.I_selected_I0_normalized.to_numpy(),
    }
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), layout="constrained")
    axes[0].errorbar(q, exp, err, fmt="k.", ms=2, alpha=.7, label="Experiment")
    for label, curve in curves.items(): axes[0].plot(q, curve, color=COLORS[label], lw=1.25, label=label)
    axes[0].set(yscale="log", xlabel="q (Å⁻¹)", ylabel="I(q), common extrapolated I(0)"); axes[0].legend(fontsize=8)
    for label, curve in curves.items(): axes[1].plot(q, (exp - curve) / err, color=COLORS[label], lw=1, label=label)
    axes[1].axhline(0, color="k", lw=.6); axes[1].axvspan(q.min(), .04, color="gold", alpha=.15); axes[1].set(xlabel="q (Å⁻¹)", ylabel="[Iexp − Icalc] / σ")
    low = q <= .04
    axes[2].plot(q[low] ** 2, np.log(exp[low]), "k.", label=f"Experiment: Rg={exp_rg:.1f} Å")
    for label, curve in curves.items():
        _, rg = guinier(q, curve)
        axes[2].plot(q[low] ** 2, np.log(curve[low]), color=COLORS[label], label=f"{label}: {rg:.1f} Å")
    axes[2].set(xlabel="q² (Å⁻²)", ylabel="ln I(q)"); axes[2].legend(fontsize=8)
    fig.suptitle("WT expanded pool: SAXS comparison after multiplicative I(0) normalization")
    fig.savefig(TRIAL / "SAXS_comparison.png", dpi=180)
    fig.savefig(TRIAL / "SAXS_comparison.pdf")
    plt.close(fig)

    pivot = metrics.pivot(index="ensemble", columns="observable", values="RMSD")
    summary_lines = [
        "# WT expanded aggressive-fit trial", "",
        "All selected ensembles contain 500 distinct structures at equal weight. Shift scales are 0.05 ppm for all carbons, 0.25 ppm for N, and 0.05 ppm for HN/HA.", "",
        "SAXS curves are multiplicatively normalized to the experimental Guinier-extrapolated I(0); no additive offset or free post-selection scale is used. Low q is q<0.04 A^-1, the unchanged middle is 0.04-0.10 A^-1, and the tail is q>=0.10 A^-1.", "",
        "| Ensemble | C RMSD | CA RMSD | N RMSD | NH RMSD | low-q chi2/point | all-q chi2/point | Guinier Rg (A) |", "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for label in ["Raw pool", "Current fit", "2x low / 0.5x tail", "4x low / 0.25x tail"]:
        r = pivot.loc[label]
        summary_lines.append(f"| {label} | {r.C:.3f} | {r.CA:.3f} | {r.N:.3f} | {r.NH:.3f} | {r.SAXS_chi2_low_q_lt_0p04:.3f} | {r.SAXS_chi2_all:.3f} | {r.calculated_Guinier_Rg_A:.2f} |")
    summary_lines += ["", "Experimental diagnostic Guinier Rg: 27.82 A.", "", "The aggressive schedules improve C and CA only slightly beyond the current fit, do not improve N or NH, and leave the selected DSSP-H profile essentially unchanged. The 4x schedule gives the best low-q SAXS agreement, but sacrifices tail agreement. These results should not yet be propagated to every pool.", ""]
    (TRIAL / "README.md").write_text("\n".join(summary_lines))


if __name__ == "__main__":
    main()

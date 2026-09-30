"""Trial aggressive CS/SAXS equal-weight selections for WT expanded pool.

This is deliberately separate from the established ASTEROIDS-style outputs.
It compares two user-requested q-weight schedules and normalizes every SAXS
candidate to the experimental Guinier-extrapolated I(0) before selection.
"""
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "work/asteroids"))
from run_selection import load
from selection import select

OUT = ROOT / "outputs/Tau5_ASTEROIDS_trials/WT_expanded_aggressive_CS_lowq"
SCALES = {"CA": 0.05, "C": 0.05, "CB": 0.05, "N": 0.25, "NH": 0.05, "HA": 0.05}
SCHEDULES = {
    "low2_tail0p5": {"low": 2.0, "middle": 1.0, "tail": 0.5},
    "low4_tail0p25": {"low": 4.0, "middle": 1.0, "tail": 0.25},
}


def guinier_i0(q, intensity, sigma=None, cutoff=0.03):
    mask = (q <= cutoff + 1e-12) & (intensity > 0)
    x = q[mask] ** 2
    y = np.log(intensity[mask])
    if sigma is None:
        weights = np.ones(mask.sum())
    else:
        weights = intensity[mask] / sigma[mask]
    slope, intercept = np.polyfit(x, y, 1, w=weights)
    return float(np.exp(intercept)), float(np.sqrt(max(0.0, -3.0 * slope)))


class FixedI0WeightedObjective:
    def __init__(self, shifts, observed, shift_scales, saxs, intensity, errors, q_weights):
        self.ncs = shifts.shape[1]
        cs = (shifts - observed) / shift_scales
        sx = saxs / errors * np.sqrt(q_weights)
        self.target = intensity / errors * np.sqrt(q_weights)
        self.features = np.concatenate([cs, sx], axis=1)

    def evaluate(self, means):
        means = np.atleast_2d(means)
        cs = means[:, :self.ncs]
        sx = means[:, self.ncs:]
        return np.sum(cs * cs, axis=1) + np.sum((sx - self.target) ** 2, axis=1)


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def write_csv(path, header, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def metrics(ids, cs, observed, atoms, q, curves, intensity, errors, q_weights):
    pred_cs = cs[ids].mean(axis=0)
    pred_i = curves[ids].mean(axis=0)
    residual = pred_cs - observed
    z = (pred_i - intensity) / errors
    result = {
        "RMSD_ppm": {a: float(np.sqrt(np.mean(residual[atoms == a] ** 2))) for a in sorted(set(atoms))},
        "SAXS_chi2_all": float(np.mean(z * z)),
        "SAXS_chi2_low_q_lt_0p04": float(np.mean(z[q < 0.04] ** 2)),
        "SAXS_chi2_middle_0p04_to_0p10": float(np.mean(z[(q >= 0.04) & (q < 0.10)] ** 2)),
        "SAXS_chi2_tail_ge_0p10": float(np.mean(z[q >= 0.10] ** 2)),
        "weighted_SAXS_sum": float(np.sum(q_weights * z * z)),
    }
    result["calculated_I0"], result["calculated_Guinier_Rg_A"] = guinier_i0(q, pred_i)
    return result


def run(generations=3000, population=180, seed=20261229):
    OUT.mkdir(parents=True, exist_ok=True)
    records, restraints, cs, rc, observed, atoms, _, q, intensity, errors, curves, fingerprint = load("WT", 1, "expanded")
    scales = np.array([SCALES[a] for a in atoms])
    exp_i0, exp_rg = guinier_i0(q, intensity, errors)
    candidate_i0 = np.array([guinier_i0(q, curve)[0] for curve in curves])
    normalized = curves * (exp_i0 / candidate_i0)[:, None]
    raw_ids = np.arange(len(records))
    manifest = {
        "sample": "WT",
        "pool": "expanded",
        "n_candidates": len(records),
        "n_selected": 500,
        "shift_scales_ppm": SCALES,
        "q_regions_Ainv": {"low": "q < 0.04", "middle": "0.04 <= q < 0.10", "tail": "q >= 0.10"},
        "q_weight_schedules": SCHEDULES,
        "I0_normalization": "Each candidate multiplied by experimental Guinier I0 / candidate Guinier I0; intercept fit over q <= 0.03 A^-1. No additive background and no free post-selection scale.",
        "experimental_I0": exp_i0,
        "experimental_Guinier_Rg_A": exp_rg,
        "population": population,
        "generations": generations,
        "seed": seed,
        "source_input_sha256": fingerprint,
        "equal_weight": 1 / 500,
    }
    save_json(OUT / "run_manifest.json", manifest)
    write_csv(OUT / "candidates.csv", ["candidate_index", "source_pdb", "protonated_pdb", "candidate_I0_before_normalization"],
              ((i, r["pdb"], r["protonated_pdb"], candidate_i0[i]) for i, r in enumerate(records)))
    dssp_h = np.array([np.array(r["dssp"]) == "H" for r in records])
    raw_cs = cs.mean(axis=0)
    raw_i = normalized.mean(axis=0)
    results = {"raw": metrics(raw_ids, cs, observed, atoms, q, normalized, intensity, errors, np.ones_like(q))}
    write_csv(OUT / "raw_shifts.csv",
              ["AR_residue", "atom", "experimental_secondary_ppm", "raw_secondary_ppm", "Exp_minus_raw_ppm"],
              ((r["author_residue"], r["atom"], observed[j] - rc[j], raw_cs[j] - rc[j], observed[j] - raw_cs[j]) for j, r in enumerate(restraints)))
    write_csv(OUT / "raw_saxs.csv", ["q_Ainv", "I_exp", "sigma_exp", "I_raw_I0_normalized"], zip(q, intensity, errors, raw_i))

    for schedule_name, schedule in SCHEDULES.items():
        dest = OUT / schedule_name
        dest.mkdir(exist_ok=True)
        q_weights = np.where(q < 0.04, schedule["low"], np.where(q < 0.10, schedule["middle"], schedule["tail"]))
        objective = FixedI0WeightedObjective(cs, observed, scales, normalized, intensity, errors, q_weights)
        started = time.time()

        def checkpoint(generation, ids, score, trace):
            if generation % 100 == 0 or generation == generations:
                save_json(dest / "checkpoint.json", {
                    "generation": generation,
                    "selected_indices": ids.tolist(),
                    "objective": float(score),
                    "trace": trace,
                })
                print(schedule_name, generation, float(score), flush=True)

        ids, trace = select(objective, 500, seed=seed, population=population, generations=generations, callback=checkpoint)
        selected_cs = cs[ids].mean(axis=0)
        selected_i = normalized[ids].mean(axis=0)
        result = metrics(ids, cs, observed, atoms, q, normalized, intensity, errors, q_weights)
        result.update({"schedule": schedule, "final_objective": float(trace[-1]), "elapsed_seconds": time.time() - started})
        results[schedule_name] = result
        write_csv(dest / "selected.csv", ["candidate_index", "source_pdb", "protonated_pdb", "weight"],
                  ((int(i), records[i]["pdb"], records[i]["protonated_pdb"], 1 / 500) for i in ids))
        write_csv(dest / "shifts.csv",
                  ["AR_residue", "atom", "experimental_secondary_ppm", "raw_secondary_ppm", "selected_secondary_ppm", "Exp_minus_raw_ppm", "Exp_minus_selected_ppm"],
                  ((r["author_residue"], r["atom"], observed[j] - rc[j], raw_cs[j] - rc[j], selected_cs[j] - rc[j], observed[j] - raw_cs[j], observed[j] - selected_cs[j]) for j, r in enumerate(restraints)))
        write_csv(dest / "saxs.csv", ["q_Ainv", "I_exp", "sigma_exp", "I_raw_I0_normalized", "I_selected_I0_normalized"],
                  zip(q, intensity, errors, raw_i, selected_i))
        write_csv(dest / "helicity.csv", ["AR_residue", "raw_DSSP_H", "selected_DSSP_H"],
                  ((j + 328, dssp_h[:, j].mean(), dssp_h[ids, j].mean()) for j in range(2, 119)))
        save_json(dest / "summary.json", result)
    save_json(OUT / "comparison.json", results)


if __name__ == "__main__":
    run()

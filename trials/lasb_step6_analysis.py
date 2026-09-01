"""
Step 6: final analysis over rmsd_results.csv (+ secondary metrics if present).

6a. Primary comparison: ensemble (A) vs apo crystal (B) top-1 success rate,
    paired across ligands, with a Wilcoxon signed-rank test on RMSD and a
    bootstrap CI on the success-rate difference. Small n -> report wide CIs
    honestly rather than a misleadingly tight point estimate.
6b. Mean RMSD alongside success rate per condition (hides whether "misses"
    are near-misses at 2.1 A or garbage at 8 A).
6c. Ceiling fraction: how much of (Condition C - Condition B) gap does the
    ensemble (A) close.
6d. Cross-ligand coverage: does one conformer serve every ligand, or does
    the winning conformer vary by ligand (already computed in Step 5c).
"""
import csv
import numpy as np
from scipy import stats
from collections import defaultdict


def bootstrap_ci_diff_proportion(a_success, b_success, n=10000, seed=0):
    """Paired bootstrap CI on mean(a_success - b_success)."""
    rng = np.random.default_rng(seed)
    a, b = np.array(a_success, dtype=float), np.array(b_success, dtype=float)
    N = len(a)
    diffs = []
    for _ in range(n):
        idx = rng.integers(0, N, N)
        diffs.append(np.mean(a[idx]) - np.mean(b[idx]))
    return np.percentile(diffs, 2.5), np.percentile(diffs, 97.5)


def main():
    with open("results/lasb_ensemble_rmsd/rmsd_results.csv") as f:
        rows = list(csv.DictReader(f))

    by_ligand = defaultdict(dict)
    for r in rows:
        by_ligand[r["ligand"]][r["condition"]] = r

    paired_ligands = [lig for lig, d in by_ligand.items() if "A" in d and "B" in d and "C" in d]
    print(f"n = {len(paired_ligands)} ligands with all three conditions present\n")

    # 6a/6b: per-condition summary
    print(f"{'Condition':<10}{'n':>4}{'top1<2A':>10}{'top1<1A':>10}{'mean top1':>12}{'mean oracle':>13}")
    cond_rows = defaultdict(list)
    for lig in paired_ligands:
        for c in "ABC":
            cond_rows[c].append(by_ligand[lig][c])
    for c in "ABC":
        rs = cond_rows[c]
        n = len(rs)
        t2 = np.mean([r["top1_success_2A"] == "True" for r in rs])
        t1 = np.mean([r["top1_success_1A"] == "True" for r in rs])
        mt = np.mean([float(r["top1_rmsd"]) for r in rs])
        mo = np.mean([float(r["oracle_rmsd"]) for r in rs])
        print(f"{c:<10}{n:>4}{t2:>10.1%}{t1:>10.1%}{mt:>12.2f}{mo:>13.2f}")

    # 6a: paired A vs B
    a_top1 = np.array([float(by_ligand[lig]["A"]["top1_rmsd"]) for lig in paired_ligands])
    b_top1 = np.array([float(by_ligand[lig]["B"]["top1_rmsd"]) for lig in paired_ligands])
    a_success = [by_ligand[lig]["A"]["top1_success_2A"] == "True" for lig in paired_ligands]
    b_success = [by_ligand[lig]["B"]["top1_success_2A"] == "True" for lig in paired_ligands]

    print(f"\nPaired A (ensemble) vs B (apo crystal), top-1 RMSD, n={len(paired_ligands)}:")
    if len(set(a_top1 - b_top1)) > 1:
        w_stat, w_p = stats.wilcoxon(a_top1, b_top1)
        print(f"  Wilcoxon signed-rank: stat={w_stat:.2f}, p={w_p:.3f}")
    diff_ci = bootstrap_ci_diff_proportion(a_success, b_success)
    print(f"  Success-rate(2A) difference (A-B): {np.mean(a_success)-np.mean(b_success):+.1%}, "
          f"95% CI [{diff_ci[0]:+.1%}, {diff_ci[1]:+.1%}]")
    if diff_ci[0] < 0 < diff_ci[1]:
        print("  CI spans zero -- not distinguishable from no difference at this n.")

    # 6c: ceiling fraction (A - B) / (C - B), on RMSD scale (lower is better, so invert sign)
    c_top1 = np.array([float(by_ligand[lig]["C"]["top1_rmsd"]) for lig in paired_ligands])
    denom = b_top1 - c_top1  # positive if crystal self-docking beats apo, as expected
    numer = b_top1 - a_top1  # positive if ensemble beats apo
    valid = denom > 1e-6
    if valid.any():
        frac = numer[valid] / denom[valid]
        print(f"\nCeiling fraction captured (mean over {valid.sum()} ligands where C beats B): "
              f"{np.mean(frac):.1%}")
    else:
        print("\nCeiling fraction: condition C did not beat B on any ligand -- cannot compute meaningfully.")

    print(f"\nVerdict: {'ensemble beats apo crystal' if np.mean(a_success) > np.mean(b_success) else 'ensemble does NOT beat apo crystal'} "
          f"on top-1 success@2A at n={len(paired_ligands)} (see CI above for whether this is distinguishable from noise).")


if __name__ == "__main__":
    main()

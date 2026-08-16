import json
import glob
import os
from collections import Counter

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy import stats

import chem_metrics
import pocket_rmsd as pocket_rmsd_mod


def heavy_atom_count(smiles):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return mol.GetNumHeavyAtoms() if mol else None

def summarize_all_results():
    # Gather all results JSONs (handles generic results.json or specific ones like results/payload_results.json)
    json_files = glob.glob("**/*results.json", recursive=True)
    if os.path.exists("results.json") and "results.json" not in json_files:
        json_files.append("results.json")
        
    if not json_files:
        print("No results.json files found.")
        return

    os.makedirs("result_summaries", exist_ok=True)

    for file_path in json_files:
        print(f"Processing {file_path}...")
        
        # Extract a protein name (e.g., "3PTB_results.json" -> "3PTB")
        basename = os.path.basename(file_path)
        protein_name = basename.replace("_results.json", "").replace("results.json", "baseline").replace(".json", "")
        if not protein_name:
            protein_name = "protein"

        # 1. Load the JSON data
        try:
            with open(file_path, "r") as f:
                data = json.load(f)
        except Exception as e:
            print(f"Failed to process {file_path}: {e}")
            continue

        if not data or not isinstance(data, list):
            print(f"No valid list data in {file_path}. Skipping.")
            continue

        # Prepare output file
        out_path = os.path.join("result_summaries", f"{protein_name}_ligand_rank.txt")
        
        with open(out_path, "w") as out:
            # 2. The ConforMix MVP Analysis
            out.write("--- Winning Conformations Breakdown ---\n")
            conformations = [row.get('winning_conformation') for row in data if row.get('winning_conformation')]
            if conformations:
                counts = Counter(conformations)
                for conf, count in counts.most_common():
                    out.write(f"{conf:<30} {count}\n")
            else:
                out.write("No 'winning_conformation' data available.\n")
            out.write("\n\n")

            # 3. Define the "Goldilocks" Thresholds
            # Strong affinity, very easy to synthesize, and highly drug-like
            elite_candidates = []
            lowest_affinity = None
            
            for row in data:
                aff = row.get('true_affinity')
                sa = row.get('sa_score')
                qed = row.get('qed')

                # Track lowest affinity for debugging context if no candidates are found
                if aff is not None:
                    if lowest_affinity is None or aff < lowest_affinity:
                        lowest_affinity = aff

                # Safely check conditions
                passes_aff = (aff is not None and aff <= -7.0)
                passes_sa = (sa is not None and sa <= 4.0)
                passes_qed = (qed is not None and qed >= 0.5)

                # Check if elements are missing entirely vs failing threshold
                # If SA or QED are missing, we still filter on true_affinity as minimal fallback
                has_filter_data = ('sa_score' in row and 'qed' in row)
                
                if has_filter_data:
                    is_elite = passes_aff and passes_sa and passes_qed
                else:
                    is_elite = passes_aff

                if is_elite:
                    elite_candidates.append(row)

            # Output elite candidates
            out.write(f"--- Top {len(elite_candidates)} Elite Candidates ---\n")
            
            if elite_candidates:
                # Sort the elite candidates by affinity (lowest first)
                elite_candidates.sort(key=lambda r: (r.get('true_affinity') is None, r.get('true_affinity')))

                # Define the columns to print
                header = f"{'Idx':<4} {'SMILES':<60} {'Affinity':<10} {'QED':<6} {'SA Score':<8}"
                out.write(header + "\n")
                out.write("-" * len(header) + "\n")
                
                for row in elite_candidates:
                    idx = row.get('original_index', 'N/A')
                    smi = row.get('smiles', '')
                    if len(smi) > 58:
                        smi = smi[:55] + "..."
                    
                    aff = f"{row.get('true_affinity'):.2f}" if row.get('true_affinity') is not None else "N/A"
                    qed = f"{row.get('qed'):.3f}" if row.get('qed') is not None else "N/A"
                    sa = f"{row.get('sa_score'):.3f}" if row.get('sa_score') is not None else "N/A"
                    
                    out.write(f"{idx:<4} {smi:<60} {aff:<10} {qed:<6} {sa:<8}\n")
            else:
                out.write("No elite candidates found matching the affinity/QED/SA_Score filters.\n")
                if lowest_affinity is not None:
                    out.write(f"\n[Info] Lowest observed affinity in dataset was {lowest_affinity:.2f} kcal/mol\n")

            out.write("\n\n")

            # 4. Ligand efficiency ranking — raw affinity alone hides good binders that are
            # small (e.g. fragment-sized DiffSBDD output capped near -5 to -6 kcal/mol despite
            # a strong per-atom contribution). LE = |affinity| / heavy_atom_count.
            le_rows = []
            for row in data:
                aff = row.get('true_affinity')
                ha = heavy_atom_count(row.get('smiles'))
                if aff is not None and ha:
                    le_rows.append((row, ha, abs(aff) / ha))

            out.write("--- Full Ranking by Ligand Efficiency (|Affinity| / Heavy Atoms) ---\n")
            if le_rows:
                le_rows.sort(key=lambda t: -t[2])
                header = f"{'ID':<12}{'Affinity':>10}{'HeavyAtoms':>12}{'LE':>8}"
                out.write(header + "\n")
                out.write("-" * len(header) + "\n")
                for row, ha, le in le_rows:
                    ident = row.get('id', row.get('original_index', 'N/A'))
                    out.write(f"{str(ident):<12}{row.get('true_affinity'):>10.2f}{ha:>12d}{le:>8.2f}\n")
            else:
                out.write("No candidates with both an affinity and a valid SMILES.\n")

            out.write("\n")
            print(f"Saved ranking to {out_path}")

RESULTS_JSON = "results.json"
RIGID_CONTROL_CSV = "rigid_control.csv"
SUMMARY_MD = "results_summary.md"
BOOTSTRAP_RESAMPLES = 10000


def _bootstrap_ci_mean(values, n_resamples=BOOTSTRAP_RESAMPLES, ci=0.95, seed=42):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    n = len(values)
    means = np.empty(n_resamples)
    for i in range(n_resamples):
        sample = rng.choice(values, size=n, replace=True)
        means[i] = sample.mean()
    alpha = (1 - ci) / 2
    lo, hi = np.quantile(means, [alpha, 1 - alpha])
    return float(lo), float(hi)


def _cohens_d_paired(a, b):
    diff = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    sd = diff.std(ddof=1)
    return float(diff.mean() / sd) if sd > 0 else float("nan")


def generate_phase8_report(results_json=RESULTS_JSON, rigid_csv=RIGID_CONTROL_CSV, output_md=SUMMARY_MD):
    """
    Compute the Phase 8 statistics (paired tests on crystal-vs-ensemble AND the
    noise-corrected rigid-vs-ensemble comparison, bootstrap CI, effect size,
    pocket-displacement correlation, distribution summaries) and write a
    ready-to-paste markdown report. Every aggregate figure carries its N and,
    where applicable, a confidence interval -- per the plan, no number here is
    hand-transcribed elsewhere; it's all recomputed from results.json/rigid_control.csv.
    """
    if not os.path.exists(results_json):
        print(f"{results_json} not found -- run ensemble_auditor.py first.")
        return
    if not os.path.exists(rigid_csv):
        print(f"{rigid_csv} not found -- run rigid_control.py first.")
        return

    with open(results_json) as f:
        results = json.load(f)
    df = pd.DataFrame(results)

    rigid_df = pd.read_csv(rigid_csv)[["ID", "rigid_max", "noise_corrected_delta"]]
    df = df.merge(rigid_df, left_on="id", right_on="ID", how="left")

    # Chemistry / drug-likeness / pose-quality columns (Phase 3), computed fresh.
    df = chem_metrics.annotate(df)

    # Pocket-lining geometry for each candidate's winning conformer (Phase 5),
    # computed fresh rather than depending on a leftover conformer_geometry.csv.
    ref_structure = pocket_rmsd_mod.PDBParser(QUIET=True).get_structure(
        "ref", pocket_rmsd_mod.TRYPSIN_PDB
    )
    ligand_centroid = pocket_rmsd_mod.get_ben_centroid(pocket_rmsd_mod.TRYPSIN_PDB)
    lining_residues = pocket_rmsd_mod.find_pocket_lining_residues(
        ref_structure, ligand_centroid, pocket_rmsd_mod.DEFAULT_CUTOFF
    )
    ref_residues = pocket_rmsd_mod._standard_ca_residues(ref_structure)
    ref_res_ids = {id(r): i for i, r in enumerate(ref_residues)}
    lining_indices = {ref_res_ids[id(r)] for r in lining_residues}

    geometry_by_conformer = {}
    for conf_path in sorted(glob.glob(os.path.join(
        pocket_rmsd_mod.ALIGNED_RECEPTORS_DIR, "conformix_var_*.pdb"
    ))):
        basename = os.path.basename(conf_path)
        geometry_by_conformer[basename] = pocket_rmsd_mod.compute_conformer_geometry(
            ref_structure, conf_path, lining_indices
        )
    df["winner_pocket_ca_rmsd"] = df["ensemble_best_conformer"].map(
        lambda c: geometry_by_conformer.get(c, {}).get("pocket_ca_rmsd") if pd.notna(c) else None
    )

    # ── Complete-case subsets ───────────────────────────────────────────────
    crystal_vs_ensemble = df.dropna(subset=["crystal_affinity", "ensemble_best_affinity"])
    noise_corrected = df.dropna(subset=["rigid_max", "ensemble_best_affinity"])

    lines = ["# Results Summary (Phase 8)\n"]
    lines.append(
        f"N = {len(df)} candidates generated and docked; "
        f"{len(crystal_vs_ensemble)} with complete crystal/ensemble data, "
        f"{len(noise_corrected)} with complete rigid-control data.\n"
    )

    # ── 1. Crystal vs ensemble (raw comparison) ─────────────────────────────
    n1 = len(crystal_vs_ensemble)
    if n1 >= 2:
        t_stat, t_p = stats.ttest_rel(
            crystal_vs_ensemble["ensemble_best_affinity"], crystal_vs_ensemble["crystal_affinity"]
        )
        try:
            w_stat, w_p = stats.wilcoxon(
                crystal_vs_ensemble["ensemble_best_affinity"], crystal_vs_ensemble["crystal_affinity"]
            )
        except ValueError:
            w_stat, w_p = float("nan"), float("nan")
        mean_delta = crystal_vs_ensemble["delta_ensemble_vs_crystal"].mean()
        lines.append(
            f"## Crystal vs. ensemble (raw, not noise-corrected)\n\n"
            f"Paired t-test (ensemble_best_affinity vs crystal_affinity), N={n1}: "
            f"t={t_stat:.3f}, p={t_p:.4f}\n\n"
            f"Wilcoxon signed-rank, N={n1}: W={w_stat:.1f}, p={w_p:.4f}\n\n"
            f"Mean delta_ensemble_vs_crystal = {mean_delta:.3f} kcal/mol (N={n1})\n"
        )

    # ── 2. Noise-corrected comparison -- THE test that matters ─────────────
    n2 = len(noise_corrected)
    if n2 >= 2:
        t_stat2, t_p2 = stats.ttest_rel(
            noise_corrected["rigid_max"], noise_corrected["ensemble_best_affinity"]
        )
        try:
            w_stat2, w_p2 = stats.wilcoxon(
                noise_corrected["rigid_max"], noise_corrected["ensemble_best_affinity"]
            )
        except ValueError:
            w_stat2, w_p2 = float("nan"), float("nan")

        deltas = noise_corrected["noise_corrected_delta"].to_numpy()
        mean_delta2 = float(deltas.mean())
        median_delta2 = float(np.median(deltas))
        ci_lo, ci_hi = _bootstrap_ci_mean(deltas)
        d = _cohens_d_paired(noise_corrected["rigid_max"], noise_corrected["ensemble_best_affinity"])
        n_beat_noise = int((deltas > 0).sum())

        ci_includes_zero = ci_lo <= 0 <= ci_hi
        interpretation = (
            "The 95% CI on the mean includes zero -- **the ensemble effect is not "
            "distinguishable from sampling noise at this sample size.**"
            if ci_includes_zero else
            "The 95% CI on the mean excludes zero."
        )

        lines.append(
            f"## Noise-corrected comparison (rigid_max_over_seeds vs. ensemble_best_affinity) "
            f"-- the test that matters\n\n"
            f"Paired t-test, N={n2}: t={t_stat2:.3f}, p={t_p2:.4f}\n\n"
            f"Wilcoxon signed-rank, N={n2}: W={w_stat2:.1f}, p={w_p2:.4f}\n\n"
            f"Mean noise_corrected_delta = {mean_delta2:.3f} kcal/mol "
            f"(median {median_delta2:.3f}), 95% bootstrap CI [{ci_lo:.3f}, {ci_hi:.3f}] "
            f"({BOOTSTRAP_RESAMPLES} resamples), N={n2}\n\n"
            f"Cohen's d (paired) = {d:.3f}\n\n"
            f"{n_beat_noise}/{n2} candidates ({100 * n_beat_noise / n2:.0f}%) beat the "
            f"noise-matched rigid baseline.\n\n"
            f"{interpretation}\n"
        )

    # ── 3. Crystal beat every conformer ─────────────────────────────────────
    if n1 >= 1:
        crystal_won = int((crystal_vs_ensemble["delta_ensemble_vs_crystal"] < 0).sum())
        lines.append(
            f"## Crystal beat every ensemble conformer\n\n"
            f"{crystal_won}/{n1} candidates ({100 * crystal_won / n1:.0f}%), N={n1}\n"
        )

    # ── 4. Pocket-displacement correlation ──────────────────────────────────
    corr_df = df.dropna(subset=["delta_ensemble_vs_crystal", "winner_pocket_ca_rmsd"])
    n3 = len(corr_df)
    if n3 >= 3:
        rho, rho_p = stats.spearmanr(corr_df["delta_ensemble_vs_crystal"], corr_df["winner_pocket_ca_rmsd"])
        lines.append(
            f"## Spearman correlation: delta_ensemble_vs_crystal vs. winner_pocket_ca_rmsd\n\n"
            f"N={n3}: rho={rho:.3f}, p={rho_p:.4f}\n"
        )

    # ── 5. Distribution summaries ────────────────────────────────────────────
    def _dist_line(col, label, fmt="{:.3f}"):
        vals = df[col].dropna()
        if len(vals) == 0:
            return f"- {label}: no data"
        return (
            f"- {label} (N={len(vals)}): mean={fmt.format(vals.mean())}, "
            f"median={fmt.format(vals.median())}, sd={fmt.format(vals.std())}, "
            f"range=[{fmt.format(vals.min())}, {fmt.format(vals.max())}]"
        )

    lines.append("## Distribution summaries\n")
    lines.append(_dist_line("ligand_efficiency", "Ligand efficiency"))
    lines.append(_dist_line("qed", "QED"))
    lines.append(_dist_line("sa_score", "SA score"))
    lines.append(_dist_line("mw", "MW", fmt="{:.1f}"))
    lines.append("")

    # ── Interpretation guardrails (matches the plan's rules for generated text) ──
    lines.append(
        "## Notes\n\n"
        "- Every statistic above reports N, the test used, and the p-value where "
        "applicable -- none should be described as significant without those.\n"
        "- Raw affinity improvements are never reported without ligand efficiency "
        "alongside them (see the distribution summary above).\n"
    )

    with open(output_md, "w") as f:
        f.write("\n".join(lines))
    print(f"Phase 8 report written -> {output_md}")


if __name__ == "__main__":
    summarize_all_results()
    generate_phase8_report()
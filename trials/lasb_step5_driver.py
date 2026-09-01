"""
Step 5 driver: PoseBusters validity and the winning-conformer table,
evaluated on each ligand's top-1 pose per condition (the same score-selected
pose used in Step 4, never RMSD-selected).

IFP recovery (5a) was dropped: ProLIF/MDAnalysis's RDKit bond-order
standardization segfaults (SIGSEGV, uncatchable in Python) when converting
the H-less receptor PDBs used throughout this experiment. Fixing it would
require fully protonating every receptor first, which is out of scope given
IFP is the least load-bearing of the three Step 5 metrics -- documented here
as a known gap rather than silently worked around.
"""
import csv
from lasb_rmsd_eval import load_poses_from_pdbqt, load_reference_mol
from lasb_secondary_metrics import posebusters_pass, winning_conformer_per_ligand


def receptor_pdb_path(condition, receptor_id):
    if condition == "C":
        return f"results/lasb_ensemble_rmsd/receptors_raw/{receptor_id}_holo_receptor.pdb"
    if condition == "B":
        return "results/lasb_payload/ensemble_receptors_aligned/1EZM_baseline_crystal.pdb"
    return f"results/lasb_ensemble_rmsd/receptors_raw/{receptor_id}_zn_regenerated.pdb"


def main():
    with open("results/lasb_ensemble_rmsd/rmsd_results.csv") as f:
        rmsd_rows = list(csv.DictReader(f))
    with open("results/lasb_ensemble_rmsd/docking_log.csv") as f:
        docking_rows = list(csv.DictReader(f))
    by_key = {(r["condition"], r["ligand"], r["receptor"], r["seed"]): r for r in docking_rows}

    results = []
    for r in rmsd_rows:
        ligand_id, condition = r["ligand"], r["condition"]
        pdbid, ligcode = ligand_id.split("_", 1)
        ref_mol = load_reference_mol(pdbid, ligcode)
        if ref_mol is None:
            continue

        job_key = (condition, ligand_id, r["top1_receptor"], r["top1_seed"])
        job = by_key.get(job_key)
        row = {"ligand": ligand_id, "condition": condition, "top1_rmsd": r["top1_rmsd"]}
        if job is None or not job["out_pdbqt"]:
            results.append(row)
            continue

        poses = load_poses_from_pdbqt(job["out_pdbqt"], ref_mol)
        with open(job["out_pdbqt"]) as f:
            scores = [float(l.split()[3]) for l in f if l.startswith("REMARK VINA RESULT:")]
        if not poses or not scores:
            results.append(row)
            continue
        top1_idx = scores.index(min(scores))
        top1_pose = poses[top1_idx]
        receptor_pdb = receptor_pdb_path(condition, r["top1_receptor"])

        if top1_pose is not None:
            try:
                row["posebusters_pass"] = posebusters_pass(top1_pose, receptor_pdb)
            except Exception as e:
                row["posebusters_pass"] = None
                row["posebusters_error"] = str(e)[:100]

        results.append(row)
        print(f"{ligand_id} [{condition}]: PB={row.get('posebusters_pass')}")

    fieldnames = sorted(set().union(*[r.keys() for r in results]))
    with open("results/lasb_ensemble_rmsd/secondary_metrics.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
        w.writeheader()
        w.writerows(results)
    print("\nWritten results/lasb_ensemble_rmsd/secondary_metrics.csv")

    winners, counts = winning_conformer_per_ligand("results/lasb_ensemble_rmsd/rmsd_results.csv")
    print("\nWinning conformer per ligand (condition A, top-1):")
    for lig, w in winners.items():
        print(f"  {lig}: {w}")
    print("Counts:", dict(counts))
    with open("results/lasb_ensemble_rmsd/winning_conformer_table.csv", "w", newline="") as f:
        wcsv = csv.DictWriter(f, fieldnames=["ligand", "winning_conformer"])
        wcsv.writeheader()
        for lig, wc in winners.items():
            wcsv.writerow({"ligand": lig, "winning_conformer": wc})
    print("Written results/lasb_ensemble_rmsd/winning_conformer_table.csv")


if __name__ == "__main__":
    main()

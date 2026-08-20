"""
Build the tidy per-pose CSV deliverable: ligand, receptor, seed, pose_rank,
vina_score, rmsd (symmetry-corrected), posebusters_pass (top-1 pose only --
PoseBusters was only run on the score-selected pose per condition in Step 5,
not all 900 individual poses across the full experiment).
"""
import csv
from lasb_rmsd_eval import load_poses_from_pdbqt, load_reference_mol, symmetry_rmsd

with open("results/lasb_ensemble_rmsd/docking_log.csv") as f:
    docking_rows = list(csv.DictReader(f))
with open("results/lasb_ensemble_rmsd/secondary_metrics.csv") as f:
    pb_by_key = {(r["ligand"], r["condition"]): r.get("posebusters_pass") for r in csv.DictReader(f)}

ref_cache = {}
rows = []
for job in docking_rows:
    if job["ok"] != "True":
        continue
    ligand_id = job["ligand"]
    pdbid, ligcode = ligand_id.split("_", 1)
    if ligand_id not in ref_cache:
        ref_cache[ligand_id] = load_reference_mol(pdbid, ligcode)
    ref_mol = ref_cache[ligand_id]
    if ref_mol is None:
        continue

    poses = load_poses_from_pdbqt(job["out_pdbqt"], ref_mol)
    with open(job["out_pdbqt"]) as f:
        scores = [float(l.split()[3]) for l in f if l.startswith("REMARK VINA RESULT:")]

    ranked = sorted(zip(scores, poses), key=lambda x: x[0])
    for rank, (score, pose) in enumerate(ranked, start=1):
        rmsd = symmetry_rmsd(ref_mol, pose) if pose is not None else None
        rows.append({
            "ligand": ligand_id, "condition": job["condition"], "receptor": job["receptor"],
            "seed": job["seed"], "pose_rank": rank, "vina_score": score,
            "rmsd": round(rmsd, 3) if rmsd is not None else "",
            "posebusters_pass": pb_by_key.get((ligand_id, job["condition"]), "") if rank == 1 else "",
        })
    print(f"{ligand_id} [{job['condition']}] receptor={job['receptor']} seed={job['seed']}: {len(ranked)} poses")

with open("results/lasb_ensemble_rmsd/tidy_all_poses.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["ligand", "condition", "receptor", "seed", "pose_rank",
                                       "vina_score", "rmsd", "posebusters_pass"])
    w.writeheader()
    w.writerows(rows)
print(f"\nWritten results/lasb_ensemble_rmsd/tidy_all_poses.csv ({len(rows)} rows)")

"""Rerun Condition C only, after fixing the receptor-alignment bug in
lasb_holo_receptor_prep.py. Same box/exhaustiveness/seeds as lasb_dock_step3.py."""
import csv
from multiprocessing import Pool
from lasb_dock_step3 import job, SEEDS

with open("results/lasb_ensemble_rmsd/ligand_prep_log.csv") as f:
    ligand_rows = [
        {"id": f"{r['pdbid']}_{r['ligand']}", "pdbid": r["pdbid"], "pdbqt": r["pdbqt_path"]}
        for r in csv.DictReader(f) if r["pdbqt_ok"] == "True"
    ]

jobs = []
for lig in ligand_rows:
    holo_pdbqt = f"results/lasb_ensemble_rmsd/receptors_raw/{lig['pdbid']}_holo_receptor.pdbqt"
    for seed in SEEDS:
        jobs.append(("C", lig["id"], lig["pdbqt"], lig["pdbid"], holo_pdbqt, seed))

print(f"{len(jobs)} Condition C docking jobs")
# Pool(8) silently produced near-zero (degenerate) scores for most Condition C
# receptors -- confirmed the same job() call gives a normal score (-5.166 for
# 6F8B) when run in isolation but ~0.0 under 8-way parallelism, pointing at
# resource contention during Vina's per-receptor grid computation (15 distinct
# large receptor files here, vs. only 3 distinct/heavily-repeated files in
# conditions A/B). Reduced concurrency to avoid it.
with Pool(processes=2) as pool:
    results = pool.map(job, jobs)

n_ok = sum(1 for r in results if r["ok"])
print(f"{n_ok}/{len(results)} succeeded")

# Merge into docking_log.csv, replacing old Condition C rows
with open("results/lasb_ensemble_rmsd/docking_log.csv") as f:
    old_rows = [r for r in csv.DictReader(f) if r["condition"] != "C"]

all_rows = old_rows + results
with open("results/lasb_ensemble_rmsd/docking_log.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader()
    w.writerows(all_rows)
print("Updated results/lasb_ensemble_rmsd/docking_log.csv")

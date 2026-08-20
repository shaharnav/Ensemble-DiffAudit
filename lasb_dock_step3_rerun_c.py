"""Rerun Condition C only, after fixing the receptor-alignment bug in
lasb_holo_receptor_prep.py. Same box/exhaustiveness/seeds as lasb_dock_step3.py."""
import csv
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
# Multiprocessing.Pool proved unreliable on this machine for this rerun (the
# first near-zero-score run turned out to be two duplicate script instances
# racing on the same output files -- a process-management bug, not a Vina
# issue; a subsequent single clean Pool run then stalled with worker
# processes spawning and dying without ever invoking vina, plausibly due to
# system memory pressure). Running serially trades speed for reliability.
results = []
for i, j in enumerate(jobs):
    r = job(j)
    results.append(r)
    print(f"[{i+1}/{len(jobs)}] {r['ligand']} seed={r['seed']}: ok={r['ok']} score={r['best_score']}", flush=True)

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

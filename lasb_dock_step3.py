"""
Step 3: matched-compute docking across conditions A (ensemble conformers),
B (apo crystal, seed-matched to the surviving ensemble), C (self-docking
ceiling, each ligand's own holo crystal).

Box: center (55.521, 35.882, 20.807), size 24 A cube -- same box for every
condition. exhaustiveness 16. Retains all output poses (num_modes 20,
energy_range 5) for the oracle metric in Step 4.

Usage: set N_SEEDS below to the number of conformers that survive Step 1
validation, then run. Condition A gets exactly one seed per conformer;
Conditions B and C reuse that same seed list, so compute is matched exactly.
"""
import csv, os, subprocess, re, itertools
from multiprocessing import Pool

VINA_BIN = "./bin/vina_1.2.7_mac_aarch64"
BOX_CENTER = (55.521, 35.882, 20.807)
BOX_SIZE = 24.0
EXHAUSTIVENESS = 16
NUM_MODES = 20
ENERGY_RANGE = 5

# Filled in once Step 1 validation determines which conformers survive.
# Each conformer gets a distinct seed; conditions B and C reuse this exact list.
CONFORMER_RECEPTORS = {
    # "beta0.0": "results/lasb_ensemble_rmsd/receptors_raw/beta0.0_prepped.pdbqt",
    # ...
}
SEEDS = list(range(1, len(CONFORMER_RECEPTORS) + 1))  # matched across A/B/C

APO_RECEPTOR = "results/lasb_ensemble_rmsd/receptors_raw/1EZM_apo_prepped.pdbqt"

OUT_DIR = "results/lasb_ensemble_rmsd/docking"
os.makedirs(OUT_DIR, exist_ok=True)


def run_vina(receptor_pdbqt, ligand_pdbqt, out_pdbqt, seed, log_path):
    cmd = [
        VINA_BIN, "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
        "--center_x", str(BOX_CENTER[0]), "--center_y", str(BOX_CENTER[1]), "--center_z", str(BOX_CENTER[2]),
        "--size_x", str(BOX_SIZE), "--size_y", str(BOX_SIZE), "--size_z", str(BOX_SIZE),
        "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES),
        "--energy_range", str(ENERGY_RANGE), "--seed", str(seed),
        "--out", out_pdbqt, "--cpu", "1",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_path, "w") as f:
        f.write(r.stdout + "\n--STDERR--\n" + r.stderr)
    return r.returncode == 0 and os.path.exists(out_pdbqt)


def parse_vina_scores(out_pdbqt):
    """Return list of (mode, affinity_kcal_mol) parsed from a multi-model output pdbqt."""
    scores = []
    if not os.path.exists(out_pdbqt):
        return scores
    with open(out_pdbqt) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                scores.append(float(parts[3]))
    return scores


def job(args):
    condition, ligand_id, ligand_pdbqt, receptor_id, receptor_pdbqt, seed = args
    tag = f"{condition}_{ligand_id}_{receptor_id}_s{seed}"
    out_pdbqt = f"{OUT_DIR}/{tag}_out.pdbqt"
    log_path = f"{OUT_DIR}/{tag}.log"
    ok = run_vina(receptor_pdbqt, ligand_pdbqt, out_pdbqt, seed, log_path)
    scores = parse_vina_scores(out_pdbqt) if ok else []
    return {
        "condition": condition, "ligand": ligand_id, "receptor": receptor_id,
        "seed": seed, "ok": ok, "out_pdbqt": out_pdbqt if ok else "",
        "n_poses": len(scores), "best_score": min(scores) if scores else None,
    }


def build_jobs(ligand_rows, n_workers_hint=8):
    jobs = []
    for lig in ligand_rows:
        ligand_id, ligand_pdbqt = lig["id"], lig["pdbqt"]

        # Condition A: each conformer, 1 seed each
        for i, (conf_id, conf_pdbqt) in enumerate(CONFORMER_RECEPTORS.items()):
            jobs.append(("A", ligand_id, ligand_pdbqt, conf_id, conf_pdbqt, SEEDS[i]))

        # Condition B: apo crystal, seed-matched to condition A's seed count
        for seed in SEEDS:
            jobs.append(("B", ligand_id, ligand_pdbqt, "1EZM_apo", APO_RECEPTOR, seed))

        # Condition C: this ligand's own holo crystal, same seed count
        holo_pdbqt = f"results/lasb_ensemble_rmsd/receptors_raw/{lig['pdbid']}_holo_receptor.pdbqt"
        for seed in SEEDS:
            jobs.append(("C", ligand_id, ligand_pdbqt, lig["pdbid"], holo_pdbqt, seed))
    return jobs


if __name__ == "__main__":
    if not CONFORMER_RECEPTORS:
        raise SystemExit(
            "CONFORMER_RECEPTORS is empty -- fill in once Step 1 validation "
            "on the Zn-regenerated conformers determines which ones survive."
        )

    with open("results/lasb_ensemble_rmsd/ligand_prep_log.csv") as f:
        ligand_rows = [
            {"id": f"{r['pdbid']}_{r['ligand']}", "pdbid": r["pdbid"], "pdbqt": r["pdbqt_path"]}
            for r in csv.DictReader(f) if r["pdbqt_ok"] == "True"
        ]

    jobs = build_jobs(ligand_rows)
    print(f"{len(jobs)} docking jobs "
          f"({len(ligand_rows)} ligands x ({len(CONFORMER_RECEPTORS)} conformers + "
          f"{len(SEEDS)} apo seeds + {len(SEEDS)} holo seeds))")

    with Pool(processes=8) as pool:
        results = pool.map(job, jobs)

    n_ok = sum(1 for r in results if r["ok"])
    print(f"{n_ok}/{len(results)} docking jobs succeeded")

    with open("results/lasb_ensemble_rmsd/docking_log.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)
    print("Written results/lasb_ensemble_rmsd/docking_log.csv")

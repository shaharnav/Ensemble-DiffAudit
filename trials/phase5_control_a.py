"""
Phase 5 Control A — Seeded Crystal Noise Floor
===============================================
Docks all 20 Phase 4 candidates into the 1EZM crystal receptor 6 times,
each with a different explicit Vina seed. Identical receptor prep, exhaustiveness,
box, and center as Phase 4. Measures how much of the ensemble gain is pure
max-over-N sampling noise from re-seeding the same rigid receptor.

Reads candidates from: results/lasb_payload/valid_candidates.sdf
Receptor:               results/lasb_payload/ensemble_receptors_aligned/1EZM_baseline_crystal.pdb
Output:                 results/phase5_control_a_results.csv
"""

import csv
import json
import logging
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from docking_engine import run_docking

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger(__name__)

# ── Config (must match Phase 4 exactly) ────────────────────────────────────
CRYSTAL_PDB     = "results/lasb_payload/ensemble_receptors_aligned/1EZM_baseline_crystal.pdb"
SDF_PATH        = "results/lasb_payload/valid_candidates.sdf"
OUTPUT_DIR      = "results/phase5_control_a"
OUTPUT_CSV      = "results/phase5_control_a_results.csv"
POCKET_CENTER   = (55.521, 35.882, 20.807)
BOX_SIZE        = [24.0, 24.0, 24.0]
EXHAUSTIVENESS  = 16
SEEDS           = [42, 123, 456, 789, 1337, 2024]
N_CONCURRENT    = 4
VINA_CPU        = 2

os.makedirs(OUTPUT_DIR, exist_ok=True)


def read_smiles_from_sdf(sdf_path):
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    smiles_list = []
    for mol in suppl:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        smiles_list.append(smi)
    return smiles_list


def run_job(candidate_idx, smiles, seed):
    job_name = f"ctrl_a_c{candidate_idx:04d}_seed{seed}"
    result = run_docking(
        pdb_file=CRYSTAL_PDB,
        smiles=smiles,
        output_dir=OUTPUT_DIR,
        job_name=job_name,
        exhaustiveness=EXHAUSTIVENESS,
        center_coords=POCKET_CENTER,
        box_size=BOX_SIZE,
        seed=seed,
        num_modes=9,
        vina_cpu=VINA_CPU,
    )
    return candidate_idx, seed, smiles, result


def main():
    if not os.path.exists(CRYSTAL_PDB):
        log.error(f"Crystal PDB not found: {CRYSTAL_PDB}")
        sys.exit(1)

    smiles_list = read_smiles_from_sdf(SDF_PATH)
    log.info(f"Loaded {len(smiles_list)} candidates from {SDF_PATH}")
    log.info(f"Crystal receptor: {CRYSTAL_PDB}")
    log.info(f"Seeds: {SEEDS}")
    log.info(f"Total jobs: {len(smiles_list)} candidates × {len(SEEDS)} seeds = {len(smiles_list)*len(SEEDS)}")

    jobs = [(i, smi, seed) for i, smi in enumerate(smiles_list) for seed in SEEDS]
    rows = []
    completed = 0

    with ThreadPoolExecutor(max_workers=N_CONCURRENT) as pool:
        futures = {pool.submit(run_job, i, smi, seed): (i, seed) for i, smi, seed in jobs}
        for fut in as_completed(futures):
            try:
                cidx, seed, smi, result = fut.result()
                affinity = result.get("affinity") if result else None
                rows.append({
                    "candidate_idx": cidx,
                    "seed": seed,
                    "smiles": smi,
                    "affinity": affinity,
                    "h_bonds": result.get("h_bonds", 0) if result else 0,
                })
                completed += 1
                status = f"{affinity:.2f}" if affinity else "FAILED"
                log.info(f"  [{completed}/{len(jobs)}] c{cidx:04d} seed={seed}: {status} kcal/mol")
            except Exception as e:
                cidx, seed = futures[fut]
                log.error(f"  Job c{cidx:04d} seed={seed} raised: {e}")
                rows.append({"candidate_idx": cidx, "seed": seed, "smiles": "", "affinity": None, "h_bonds": 0})
                completed += 1

    # Sort and write CSV
    rows.sort(key=lambda r: (r["candidate_idx"], r["seed"]))
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["candidate_idx", "seed", "affinity", "h_bonds", "smiles"])
        writer.writeheader()
        writer.writerows(rows)
    log.info(f"Results written to {OUTPUT_CSV}")

    # Quick summary
    failures = sum(1 for r in rows if r["affinity"] is None)
    log.info(f"Completed: {len(rows) - failures}/{len(rows)} succeeded, {failures} failed")


if __name__ == "__main__":
    main()

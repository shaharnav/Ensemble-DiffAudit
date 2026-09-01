"""
Phase 6 — Dock Known LasB Actives
===================================
Docks 11 known LasB inhibitors (purified enzyme IC50/Ki) into:
  - 3DBK holo crystal (RDF+SO4 stripped, Zn/Ca/HOH retained)
  - 1EZM apo crystal (same as used in Phases 4-5)

3 seeds per compound per receptor = 66 total jobs.
Reports mean best score across seeds.

Identical settings to Phases 4-5 (exhaustiveness=16, same box/center approach).
"""

import csv, os, sys, statistics
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from docking_engine import run_docking

SEEDS = [42, 123, 456]

RECEPTORS = {
    "3DBK_holo": {
        "pdb": "results/3DBK_prepped.pdb",
        "center": (18.721, -5.093, 23.685),
    },
    "1EZM_apo": {
        "pdb": "results/lasb_payload/ensemble_receptors_aligned/1EZM_baseline_crystal.pdb",
        "center": (55.521, 35.882, 20.807),
    },
}

BOX_SIZE      = [24.0, 24.0, 24.0]
EXHAUSTIVENESS = 16
N_CONCURRENT   = 4
VINA_CPU       = 2

OUTPUT_DIR = "results/phase6_docking"
OUTPUT_CSV = "results/phase6_docking_scores.csv"
os.makedirs(OUTPUT_DIR, exist_ok=True)

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

# Load compounds
compounds = []
with open('results/phase6_compounds.csv') as f:
    for r in csv.DictReader(f):
        compounds.append(r)

log.info(f"{len(compounds)} compounds × {len(RECEPTORS)} receptors × {len(SEEDS)} seeds = "
         f"{len(compounds)*len(RECEPTORS)*len(SEEDS)} jobs")

def run_job(cmpd, rec_name, rec_cfg, seed):
    result = run_docking(
        pdb_file=rec_cfg["pdb"],
        smiles=cmpd["smiles"],
        output_dir=OUTPUT_DIR,
        job_name=f"p6_{cmpd['id']}_{rec_name}_s{seed}",
        exhaustiveness=EXHAUSTIVENESS,
        center_coords=rec_cfg["center"],
        box_size=BOX_SIZE,
        seed=seed,
        num_modes=9,
        vina_cpu=VINA_CPU,
    )
    return cmpd["id"], rec_name, seed, result

jobs = [
    (c, rn, rc, s)
    for c in compounds
    for rn, rc in RECEPTORS.items()
    for s in SEEDS
]

raw_rows = []
done = 0

with ThreadPoolExecutor(max_workers=N_CONCURRENT) as pool:
    futures = {pool.submit(run_job, c, rn, rc, s): (c["id"], rn, s)
               for c, rn, rc, s in jobs}
    for fut in as_completed(futures):
        cid, rec, s = futures[fut]
        try:
            cid, rec, seed, result = fut.result()
            aff = result.get("affinity") if result else None
            raw_rows.append({"id": cid, "receptor": rec, "seed": seed, "affinity": aff})
            done += 1
            status = f"{aff:.2f}" if aff else "FAILED"
            log.info(f"  [{done}/{len(jobs)}] {cid} | {rec} | seed={seed}: {status}")
        except Exception as e:
            log.error(f"  {cid}/{rec}/s{s}: {e}")
            raw_rows.append({"id": cid, "receptor": rec, "seed": s, "affinity": None})
            done += 1

# Aggregate: mean best score across seeds
from collections import defaultdict
agg = defaultdict(list)
for r in raw_rows:
    if r["affinity"] is not None:
        agg[(r["id"], r["receptor"])].append(r["affinity"])

# Load compound metadata
meta = {r["id"]: r for r in compounds}

with open(OUTPUT_CSV, "w", newline="") as f:
    fields = ["id", "receptor", "pIC50", "zbg_class", "mw", "heavy_atoms",
              "mean_affinity", "min_affinity", "n_seeds", "smiles", "source"]
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for (cid, rec), affs in sorted(agg.items()):
        m = meta[cid]
        w.writerow({
            "id": cid, "receptor": rec,
            "pIC50": m["pIC50"], "zbg_class": m["zbg_class"],
            "mw": m["mw"], "heavy_atoms": m["heavy_atoms"],
            "mean_affinity": round(statistics.mean(affs), 4),
            "min_affinity": round(min(affs), 4),
            "n_seeds": len(affs),
            "smiles": m["smiles"], "source": m["source"],
        })

failures = sum(1 for r in raw_rows if r["affinity"] is None)
log.info(f"Done. {len(raw_rows)-failures}/{len(raw_rows)} succeeded. Results → {OUTPUT_CSV}")

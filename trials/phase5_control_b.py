"""
Phase 5 Control B — Decoy Docking
===================================
Docks 20 property-matched PubChem decoys into the same 6 receptors used in
Phase 4 (5 ConforMix conformers + 1EZM crystal), with identical settings.

Decoys: results/phase5_decoys.sdf (20 PubChem compounds, Tanimoto < 0.3)
Receptors: results/lasb_payload/ensemble_receptors_aligned/*.pdb
Output:   results/phase5_control_b_results.csv
"""

import csv, os, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from docking_engine import run_docking

RECEPTOR_DIR  = "results/lasb_payload/ensemble_receptors_aligned"
DECOY_SDF     = "results/phase5_decoys.sdf"
OUTPUT_DIR    = "results/phase5_control_b"
OUTPUT_CSV    = "results/phase5_control_b_results.csv"
POCKET_CENTER = (55.521, 35.882, 20.807)
BOX_SIZE      = [24.0, 24.0, 24.0]
EXHAUSTIVENESS = 16
SEED           = 42
N_CONCURRENT   = 4
VINA_CPU       = 2

CONFORMERS = [
    "lasb_conformer_beta0.0.pdb",
    "lasb_conformer_beta0.8.pdb",
    "lasb_conformer_beta1.6.pdb",
    "lasb_conformer_beta2.4.pdb",
    "lasb_conformer_beta3.2.pdb",
]
CRYSTAL = "1EZM_baseline_crystal.pdb"

os.makedirs(OUTPUT_DIR, exist_ok=True)

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger(__name__)

# Load decoys
suppl = Chem.SDMolSupplier(DECOY_SDF, removeHs=False)
decoys = []
for mol in suppl:
    if mol is None:
        continue
    name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"decoy_{len(decoys)}"
    cand_idx = int(mol.GetProp('cand_idx')) if mol.HasProp('cand_idx') else len(decoys)
    smi = Chem.MolToSmiles(mol)
    decoys.append({'name': name, 'cand_idx': cand_idx, 'smiles': smi})

receptors = [CRYSTAL] + CONFORMERS
log.info(f"Decoys: {len(decoys)}, Receptors: {len(receptors)}, Jobs: {len(decoys)*len(receptors)}")

def run_job(d_idx, decoy, receptor):
    pdb_path = os.path.join(RECEPTOR_DIR, receptor)
    rec_label = receptor.replace('lasb_conformer_', 'beta').replace('.pdb', '').replace('1EZM_baseline_crystal', 'crystal')
    job_name = f"ctrlb_d{d_idx:04d}_{rec_label}"
    result = run_docking(
        pdb_file=pdb_path, smiles=decoy['smiles'],
        output_dir=OUTPUT_DIR, job_name=job_name,
        exhaustiveness=EXHAUSTIVENESS, center_coords=POCKET_CENTER,
        box_size=BOX_SIZE, seed=SEED, num_modes=9, vina_cpu=VINA_CPU,
    )
    return d_idx, receptor, decoy, result

jobs = [(i, d, r) for i, d in enumerate(decoys) for r in receptors]
rows = []
completed = 0

with ThreadPoolExecutor(max_workers=N_CONCURRENT) as pool:
    futures = {pool.submit(run_job, i, d, r): (i, r) for i, d, r in jobs}
    for fut in as_completed(futures):
        try:
            d_idx, receptor, decoy, result = fut.result()
            aff = result.get('affinity') if result else None
            rows.append({
                'decoy_idx': d_idx, 'pubchem_cid': decoy['name'],
                'cand_idx': decoy['cand_idx'], 'receptor': receptor,
                'smiles': decoy['smiles'],
                'affinity': aff,
                'h_bonds': result.get('h_bonds', 0) if result else 0,
            })
            completed += 1
            status = f"{aff:.2f}" if aff else "FAILED"
            log.info(f"  [{completed}/{len(jobs)}] {decoy['name']} vs {receptor}: {status}")
        except Exception as e:
            d_idx, r = futures[fut]
            log.error(f"  Job d{d_idx:04d}/{r} raised: {e}")
            rows.append({'decoy_idx': d_idx, 'receptor': r, 'affinity': None, 'h_bonds': 0,
                         'pubchem_cid': '', 'cand_idx': -1, 'smiles': ''})
            completed += 1

rows.sort(key=lambda r: (r['decoy_idx'], r['receptor']))
with open(OUTPUT_CSV, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['decoy_idx','pubchem_cid','cand_idx','receptor','affinity','h_bonds','smiles'])
    writer.writeheader()
    writer.writerows(rows)

failures = sum(1 for r in rows if r['affinity'] is None)
log.info(f"Done. {len(rows)-failures}/{len(rows)} succeeded. Results → {OUTPUT_CSV}")

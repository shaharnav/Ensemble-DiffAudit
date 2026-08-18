"""
Batch-extract (protein, ligand) for all 285 CASF-2016-derived complexes.
Parallelized with multiprocessing since CCD lookups + RDKit template matching
are independent per complex. Logs failures rather than silently dropping them.
"""
import csv, os, time
from multiprocessing import Pool
from casf_pipeline.extract import extract_ligand_mol, parse_pdb_hetatm_groups, pick_ligand_group, extract_protein

STRUCT_DIR = "results/casf2016/structures"
LIG_DIR = "results/casf2016/ligands"
PROT_DIR = "results/casf2016/proteins"
os.makedirs(LIG_DIR, exist_ok=True)
os.makedirs(PROT_DIR, exist_ok=True)

def process(pdbid):
    pdb_path = os.path.join(STRUCT_DIR, f"{pdbid}.pdb")
    sdf_path = os.path.join(LIG_DIR, f"{pdbid}_ligand.sdf")
    prot_path = os.path.join(PROT_DIR, f"{pdbid}_protein.pdb")
    groups = parse_pdb_hetatm_groups(pdb_path)
    picked = pick_ligand_group(groups)
    if picked is None:
        return pdbid, "FAIL_no_ligand_group", None
    res = extract_ligand_mol(pdb_path, sdf_path)
    if res is None:
        return pdbid, "FAIL_bond_order_assignment", picked[0][3]
    extract_protein(pdb_path, prot_path, picked[0])
    return pdbid, "ok", res[0]

if __name__ == "__main__":
    with open("results/casf2016/casf2016_core_pdbid_pkd.csv") as f:
        ids = [r["pdbid"] for r in csv.DictReader(f)]

    t0 = time.time()
    with Pool(processes=8) as pool:
        results = pool.map(process, ids)
    elapsed = time.time() - t0

    ok = [r for r in results if r[1] == "ok"]
    fail = [r for r in results if r[1] != "ok"]
    print(f"Extracted {len(ok)}/{len(ids)} in {elapsed:.1f}s ({len(fail)} failures)")

    with open("results/casf2016/extraction_log.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["pdbid","status","ligand_resname"]); w.writerows(results)

    print("Failure reasons:")
    from collections import Counter
    print(Counter(r[1] for r in fail))

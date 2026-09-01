"""
Stage 0 -- bias audit. Prep receptor+ligand pdbqt for all extracted complexes,
score with Vina and Vinardo (native vina 1.2.7 binary, --score_only), compute
HAC, and write a tidy scores CSV. Parallelized with multiprocessing.
"""
import os, csv, re, subprocess, sys, time
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from casf_pipeline.prep import prepare_receptor_pdbqt, prepare_ligand_pdbqt_from_sdf, ligand_bbox_center_size
from rdkit import Chem

VINA = os.path.abspath("./bin/vina_1.2.7_mac_aarch64")
LIG_DIR = "results/casf2016/ligands"
PROT_DIR = "results/casf2016/proteins"
PDBQT_DIR = "results/casf2016/pdbqt"
os.makedirs(PDBQT_DIR, exist_ok=True)

def score_only(receptor_pdbqt, ligand_pdbqt, center, size, scoring):
    cmd = [VINA, "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--scoring", scoring, "--score_only"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.\d+)", r.stdout)
    return float(m.group(1)) if m else None

def process(pdbid):
    sdf = os.path.join(LIG_DIR, f"{pdbid}_ligand.sdf")
    prot = os.path.join(PROT_DIR, f"{pdbid}_protein.pdb")
    if not (os.path.exists(sdf) and os.path.exists(prot)):
        return {"pdbid": pdbid, "status": "FAIL_missing_input"}

    rec_pdbqt = os.path.join(PDBQT_DIR, f"{pdbid}_receptor.pdbqt")
    lig_pdbqt = os.path.join(PDBQT_DIR, f"{pdbid}_ligand.pdbqt")

    try:
        if not prepare_receptor_pdbqt(prot, rec_pdbqt):
            return {"pdbid": pdbid, "status": "FAIL_receptor_prep"}
        if not prepare_ligand_pdbqt_from_sdf(sdf, lig_pdbqt):
            return {"pdbid": pdbid, "status": "FAIL_ligand_prep"}
    except Exception as e:
        return {"pdbid": pdbid, "status": f"FAIL_exception_prep: {e}"}

    center, size = ligand_bbox_center_size(sdf)
    try:
        vina_score = score_only(rec_pdbqt, lig_pdbqt, center, size, "vina")
        vinardo_score = score_only(rec_pdbqt, lig_pdbqt, center, size, "vinardo")
    except Exception as e:
        return {"pdbid": pdbid, "status": f"FAIL_exception_score: {e}"}

    if vina_score is None or vinardo_score is None:
        return {"pdbid": pdbid, "status": "FAIL_score_parse"}

    mol = Chem.SDMolSupplier(sdf, removeHs=False)[0]
    hac = mol.GetNumHeavyAtoms() if mol else None

    return {"pdbid": pdbid, "status": "ok", "vina": vina_score, "vinardo": vinardo_score, "hac": hac}

if __name__ == "__main__":
    with open("results/casf2016/extraction_log.csv") as f:
        ok_ids = [r["pdbid"] for r in csv.DictReader(f) if r["status"] == "ok"]

    print(f"Scoring {len(ok_ids)} complexes...")
    t0 = time.time()
    with Pool(processes=8) as pool:
        results = pool.map(process, ok_ids)
    elapsed = time.time() - t0

    ok = [r for r in results if r["status"] == "ok"]
    print(f"Scored {len(ok)}/{len(ok_ids)} in {elapsed:.1f}s ({elapsed/len(ok_ids):.2f}s/complex)")

    from collections import Counter
    print(Counter(r["status"] for r in results))

    pkd = {}
    with open("results/casf2016/casf2016_core_pdbid_pkd.csv") as f:
        for r in csv.DictReader(f):
            pkd[r["pdbid"]] = float(r["-logKd/Ki"])

    with open("results/casf2016/stage0_scores.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pdbid","pKd","vina","vinardo","hac","status"])
        w.writeheader()
        for r in results:
            w.writerow({
                "pdbid": r["pdbid"], "pKd": pkd.get(r["pdbid"]),
                "vina": r.get("vina"), "vinardo": r.get("vinardo"), "hac": r.get("hac"),
                "status": r["status"],
            })
    print("Written results/casf2016/stage0_scores.csv")

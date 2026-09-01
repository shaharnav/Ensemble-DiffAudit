"""
Phase 7 Step 8b -- PoseBusters validity on all existing poses. No new docking.
"""
import os, glob, subprocess, json, re
import pandas as pd
from posebusters import PoseBusters
from rdkit import Chem

DIRS = {
    "phase4_ensemble": "results/ensemble_audit",
    "control_a": "results/phase5_control_a",
    "control_b": "results/phase5_control_b",
    "phase6_actives": "results/phase6_docking",
}

TMP = "/tmp/p7_pb"
os.makedirs(TMP, exist_ok=True)

def extract_best_pose(out_pdbqt, dest):
    with open(out_pdbqt) as f:
        lines = f.readlines()
    keep, in_model1, seen = [], False, False
    for ln in lines:
        if ln.startswith("MODEL"):
            if not seen:
                in_model1 = True; seen = True
            else:
                break
            continue
        if ln.startswith("ENDMDL"):
            if in_model1: break
            continue
        if in_model1: keep.append(ln)
    if not keep:
        keep = [l for l in lines if not l.startswith(("MODEL","ENDMDL"))]
    with open(dest, "w") as f:
        f.writelines(keep)

pb = PoseBusters(config="dock")
rows = []

for label, d in DIRS.items():
    out_files = sorted(glob.glob(os.path.join(d, "*_out.pdbqt")))
    print(f"{label}: {len(out_files)} poses")
    for outf in out_files:
        base = outf[:-len("_out.pdbqt")]
        job = os.path.basename(base)
        recf = base + "_receptor.pdbqt"
        if not os.path.exists(recf):
            rows.append({"label": label, "job": job, "error": "missing_receptor"})
            continue
        lig_pdbqt = os.path.join(TMP, "lig.pdbqt")
        extract_best_pose(outf, lig_pdbqt)
        lig_sdf = os.path.join(TMP, "lig.sdf")
        rec_pdb = os.path.join(TMP, "rec.pdb")
        r1 = subprocess.run(["obabel", lig_pdbqt, "-O", lig_sdf], capture_output=True, text=True)
        r2 = subprocess.run(["obabel", recf, "-O", rec_pdb], capture_output=True, text=True)
        if not (os.path.exists(lig_sdf) and os.path.exists(rec_pdb)):
            rows.append({"label": label, "job": job, "error": "obabel_conversion_failed"})
            continue
        try:
            mol_pred = Chem.MolFromMolFile(lig_sdf, sanitize=False)
            mol_cond = Chem.MolFromPDBFile(rec_pdb, sanitize=False)
            if mol_pred is None or mol_cond is None:
                rows.append({"label": label, "job": job, "error": "rdkit_parse_failed"})
                continue
            res = pb.bust([mol_pred], mol_cond=mol_cond, full_report=False)
            d_res = res.iloc[0].to_dict()
            d_res["label"] = label
            d_res["job"] = job
            d_res["error"] = None
            rows.append(d_res)
        except Exception as e:
            rows.append({"label": label, "job": job, "error": f"posebusters_exception: {e}"})

df = pd.DataFrame(rows)
df.to_csv("results/phase7_posebusters.csv", index=False)
print(f"Done. {len(df)} rows written to results/phase7_posebusters.csv")
print("Errors:", df['error'].notna().sum() if 'error' in df else "n/a")

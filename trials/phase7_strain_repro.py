"""
Phase 7 Step 8c (torsional strain) + 8d (pose reproducibility).
No new docking -- uses existing Control A (6-seed) poses and phase4 ensemble poses.
"""
import os, glob, subprocess, re
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

TMP = "/tmp/p7_strain"
os.makedirs(TMP, exist_ok=True)

def extract_best_pose(out_pdbqt, dest):
    with open(out_pdbqt) as f:
        lines = f.readlines()
    keep, in_model1, seen = [], False, False
    for ln in lines:
        if ln.startswith("MODEL"):
            if not seen: in_model1=True; seen=True
            else: break
            continue
        if ln.startswith("ENDMDL"):
            if in_model1: break
            continue
        if in_model1: keep.append(ln)
    if not keep:
        keep = [l for l in lines if not l.startswith(("MODEL","ENDMDL"))]
    with open(dest,"w") as f: f.writelines(keep)

def pdbqt_to_mol(pdbqt_path):
    sdf = os.path.join(TMP, "tmp.sdf")
    r = subprocess.run(["obabel", pdbqt_path, "-O", sdf], capture_output=True, text=True)
    mol = Chem.MolFromMolFile(sdf, sanitize=True)
    return mol

def mmff_strain(mol):
    """docked-conformer MMFF94 energy minus freely-minimized conformer energy."""
    try:
        molH = Chem.AddHs(mol, addCoords=True)
        props = AllChem.MMFFGetMoleculeProperties(molH)
        if props is None:
            return None
        ff_docked = AllChem.MMFFGetMoleculeForceField(molH, props)
        if ff_docked is None:
            return None
        e_docked = ff_docked.CalcEnergy()
        mol_min = Chem.Mol(molH)
        ff_min = AllChem.MMFFGetMoleculeForceField(mol_min, props)
        ff_min.Minimize(maxIts=2000)
        e_min = ff_min.CalcEnergy()
        return e_docked - e_min
    except Exception:
        return None

# --- Step 8c: strain for phase4 ensemble (candidates) and control_b (decoys) ---
strain_rows = []
for label, d in [("phase4_ensemble","results/ensemble_audit"), ("control_b","results/phase5_control_b")]:
    out_files = sorted(glob.glob(os.path.join(d, "*_out.pdbqt")))
    for outf in out_files:
        base = outf[:-len("_out.pdbqt")]
        job = os.path.basename(base)
        pose_pdbqt = os.path.join(TMP,"pose.pdbqt")
        extract_best_pose(outf, pose_pdbqt)
        mol = pdbqt_to_mol(pose_pdbqt)
        strain = mmff_strain(mol) if mol is not None else None
        rec = re.sub(r'^ens_s\d+_|^ctrlb_d\d+_', '', job)
        strain_rows.append({"label": label, "job": job, "receptor": rec, "strain_kcal_mol": strain})

df_strain = pd.DataFrame(strain_rows)
df_strain.to_csv("results/phase7_strain.csv", index=False)
print("Strain done:", len(df_strain), "rows ->", df_strain['strain_kcal_mol'].notna().sum(), "succeeded")

# --- Step 8d: pose reproducibility (Control A, 6 seeds) ---
ca_files = sorted(glob.glob("results/phase5_control_a/*_out.pdbqt"))
ca_map = {}
for outf in ca_files:
    base = outf[:-len("_out.pdbqt")]
    job = os.path.basename(base)
    m = re.match(r'ctrl_a_c(\d+)_seed(\d+)', job)
    cand_idx, seed = int(m.group(1)), int(m.group(2))
    pose_pdbqt = os.path.join(TMP, f"ca_{cand_idx}_{seed}.pdbqt")
    extract_best_pose(outf, pose_pdbqt)
    ca_map.setdefault(cand_idx, {})[seed] = pose_pdbqt

def get_coords(pdbqt_path):
    coords = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM","HETATM")):
                try:
                    x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    coords.append((x,y,z))
                except ValueError:
                    pass
    return np.array(coords)

def rmsd(a, b):
    if a.shape != b.shape:
        return None
    return np.sqrt(np.mean(np.sum((a-b)**2, axis=1)))

repro_rows = []
for cand_idx, seeds in ca_map.items():
    seed_ids = sorted(seeds.keys())
    coords = {s: get_coords(seeds[s]) for s in seed_ids}
    pairwise = []
    for i in range(len(seed_ids)):
        for j in range(i+1, len(seed_ids)):
            r = rmsd(coords[seed_ids[i]], coords[seed_ids[j]])
            if r is not None:
                pairwise.append(r)
    if pairwise:
        repro_rows.append({"cand_idx": cand_idx, "mean_pairwise_rmsd": np.mean(pairwise), "max_pairwise_rmsd": np.max(pairwise)})

df_repro = pd.DataFrame(repro_rows)
df_repro.to_csv("results/phase7_reproducibility.csv", index=False)
print("Reproducibility done:", len(df_repro), "candidates")

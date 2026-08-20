"""
Step 4: symmetry-corrected RMSD evaluation of docking_log.csv (from
lasb_dock_step3.py) against the crystallographic reference pose for each
ligand.

Uses RDKit's GetBestRMS (symmetry-corrected: tries symmetry-equivalent atom
mappings, e.g. phenyl ring flips, so a physically-identical pose doesn't get
penalized for atom-index bookkeeping).

Per ligand x condition, reports two numbers:
  TOP-1  -- pose selected by Vina score (never by RMSD), then measured
  ORACLE -- best RMSD across every pose in every receptor/seed in that
            condition (the sampling ceiling)
Then success rates at 2.0 A and 1.0 A for both.
"""
import csv, os
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem


def load_reference_mol(pdbid, ligand):
    sdf = f"results/lasb_ensemble_rmsd/ligands_sdf/{pdbid}_{ligand}.sdf"
    return Chem.SDMolSupplier(sdf, removeHs=False)[0]


def load_poses_from_pdbqt(pdbqt_path, template_mol):
    """Parse a multi-model Vina output pdbqt into a list of RDKit mols with
    the docked coordinates, using template_mol for correct bond orders
    (Vina/meeko pdbqt loses bond-order info in the ROOT/torsion-tree format)."""
    mols = []
    if not os.path.exists(pdbqt_path):
        return mols
    with open(pdbqt_path) as f:
        blocks, current = [], []
        for line in f:
            if line.startswith("MODEL"):
                current = []
            elif line.startswith("ENDMDL"):
                blocks.append(current)
            elif line.startswith(("ATOM", "HETATM")):
                current.append(line)
    for block in blocks:
        coords = []
        for line in block:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            coords.append((x, y, z))
        if len(coords) != template_mol.GetNumAtoms():
            # pdbqt heavy-atom-only vs template with Hs -- filter template to heavy atoms
            heavy_idx = [a.GetIdx() for a in template_mol.GetAtoms() if a.GetSymbol() != "H"]
            if len(coords) != len(heavy_idx):
                mols.append(None)
                continue
            mol = Chem.RWMol(template_mol)
            for i in sorted([a.GetIdx() for a in template_mol.GetAtoms() if a.GetSymbol() == "H"], reverse=True):
                mol.RemoveAtom(i)
            mol = mol.GetMol()
        else:
            mol = Chem.Mol(template_mol)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, (x, y, z))
        mol.RemoveAllConformers()
        mol.AddConformer(conf)
        mols.append(mol)
    return mols


def symmetry_rmsd(ref_mol, pose_mol):
    ref_noH = Chem.RemoveHs(ref_mol)
    pose_noH = Chem.RemoveHs(pose_mol) if pose_mol.GetNumAtoms() != ref_noH.GetNumAtoms() else pose_mol
    try:
        return AllChem.GetBestRMS(pose_noH, ref_noH)
    except Exception:
        return None


def main():
    with open("results/lasb_ensemble_rmsd/docking_log.csv") as f:
        docking_rows = list(csv.DictReader(f))

    by_ligand_condition = defaultdict(list)
    for r in docking_rows:
        if r["ok"] != "True":
            continue
        by_ligand_condition[(r["ligand"], r["condition"])].append(r)

    ref_cache = {}
    results = []
    for (ligand_id, condition), rows in by_ligand_condition.items():
        pdbid, ligcode = ligand_id.split("_", 1)
        if ligand_id not in ref_cache:
            ref_cache[ligand_id] = load_reference_mol(pdbid, ligcode)
        ref_mol = ref_cache[ligand_id]
        if ref_mol is None:
            continue

        all_poses = []  # (rmsd, vina_score, receptor, seed)
        for r in rows:
            scores_path = r["out_pdbqt"]
            poses = load_poses_from_pdbqt(scores_path, ref_mol)
            with open(scores_path) as f:
                scores = [float(l.split()[3]) for l in f if l.startswith("REMARK VINA RESULT:")]
            for pose_mol, score in zip(poses, scores):
                if pose_mol is None:
                    continue
                rmsd = symmetry_rmsd(ref_mol, pose_mol)
                if rmsd is not None:
                    all_poses.append((rmsd, score, r["receptor"], r["seed"]))

        if not all_poses:
            continue

        top1 = min(all_poses, key=lambda p: p[1])  # best (most negative) Vina score
        oracle = min(all_poses, key=lambda p: p[0])  # best RMSD achievable

        results.append({
            "ligand": ligand_id, "condition": condition, "n_poses": len(all_poses),
            "top1_rmsd": round(top1[0], 3), "top1_score": top1[1],
            "top1_receptor": top1[2], "top1_seed": top1[3],
            "oracle_rmsd": round(oracle[0], 3), "oracle_score": oracle[1],
            "oracle_receptor": oracle[2], "oracle_seed": oracle[3],
            "top1_success_2A": top1[0] < 2.0, "top1_success_1A": top1[0] < 1.0,
            "oracle_success_2A": oracle[0] < 2.0, "oracle_success_1A": oracle[0] < 1.0,
        })

    with open("results/lasb_ensemble_rmsd/rmsd_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)
    print(f"Written results/lasb_ensemble_rmsd/rmsd_results.csv ({len(results)} ligand x condition rows)")

    # Summary success rates per condition
    by_cond = defaultdict(list)
    for r in results:
        by_cond[r["condition"]].append(r)
    print(f"\n{'Condition':<12}{'n':>4}{'top1<2A':>10}{'top1<1A':>10}{'oracle<2A':>11}{'oracle<1A':>11}{'mean top1':>12}{'mean oracle':>13}")
    for cond, rows in sorted(by_cond.items()):
        n = len(rows)
        t2 = sum(r["top1_success_2A"] for r in rows) / n
        t1 = sum(r["top1_success_1A"] for r in rows) / n
        o2 = sum(r["oracle_success_2A"] for r in rows) / n
        o1 = sum(r["oracle_success_1A"] for r in rows) / n
        mt = sum(r["top1_rmsd"] for r in rows) / n
        mo = sum(r["oracle_rmsd"] for r in rows) / n
        print(f"{cond:<12}{n:>4}{t2:>10.1%}{t1:>10.1%}{o2:>11.1%}{o1:>11.1%}{mt:>12.2f}{mo:>13.2f}")


if __name__ == "__main__":
    main()

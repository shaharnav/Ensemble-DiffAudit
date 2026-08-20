"""
Step 5: secondary metrics on the docking poses already produced by Step 3.

5a. Interaction fingerprint (ProLIF) recovery: native holo IFP (computed from
    the crystallographic ligand + its own holo receptor) vs each docked pose's
    IFP, Tanimoto similarity, reported per condition.
5b. PoseBusters validity pass rate per condition.
5c. Per-ligand winning conformer (which beta value produces the best top-1
    pose in condition A) -- one winner throughout vs different winners per
    ligand distinguishes "ensemble does real work" from "one better receptor."
"""
import csv, os
from collections import defaultdict
import MDAnalysis as mda
import prolif as plf
from posebusters import PoseBusters
from rdkit import Chem


def native_ifp(holo_receptor_pdb, ligand_sdf):
    u = mda.Universe(holo_receptor_pdb)
    protein = u.select_atoms("protein")
    lig_mol = Chem.SDMolSupplier(ligand_sdf, removeHs=False)[0]
    lig_plf = plf.Molecule.from_rdkit(lig_mol)
    prot_plf = plf.Molecule.from_mda(protein)
    fp = plf.Fingerprint()
    fp.run_from_iterable([lig_plf], prot_plf)
    return fp.to_dataframe().iloc[0]


def pose_ifp(receptor_pdbqt_or_pdb, pose_mol):
    u = mda.Universe(receptor_pdbqt_or_pdb)
    protein = u.select_atoms("protein")
    lig_plf = plf.Molecule.from_rdkit(pose_mol)
    prot_plf = plf.Molecule.from_mda(protein)
    fp = plf.Fingerprint()
    fp.run_from_iterable([lig_plf], prot_plf)
    return fp.to_dataframe().iloc[0]


def ifp_tanimoto(ref_series, pose_series):
    ref_keys = set(k for k, v in ref_series.items() if v)
    pose_keys = set(k for k, v in pose_series.items() if v)
    if not ref_keys and not pose_keys:
        return 1.0
    if not ref_keys or not pose_keys:
        return 0.0
    inter = len(ref_keys & pose_keys)
    union = len(ref_keys | pose_keys)
    return inter / union if union else 0.0


def posebusters_pass(pose_mol, receptor_pdb):
    pb = PoseBusters(config="dock")
    try:
        df = pb.bust([pose_mol], None, receptor_pdb)
        row = df.iloc[0]
        return bool(row.drop(labels=[c for c in row.index if c in ("molecule", "file")]).all())
    except Exception:
        return None


def winning_conformer_per_ligand(rmsd_csv):
    """From condition A's top-1 rows only, which receptor (beta value) wins per ligand."""
    with open(rmsd_csv) as f:
        rows = [r for r in csv.DictReader(f) if r["condition"] == "A"]
    winners = {r["ligand"]: r["top1_receptor"] for r in rows}
    counts = defaultdict(int)
    for w in winners.values():
        counts[w] += 1
    return winners, counts


if __name__ == "__main__":
    print("This module provides reusable functions for Step 5. Run after Step 3/4")
    print("produce docking_log.csv and rmsd_results.csv, with a driver script that")
    print("loads poses via lasb_rmsd_eval.load_poses_from_pdbqt and calls these.")

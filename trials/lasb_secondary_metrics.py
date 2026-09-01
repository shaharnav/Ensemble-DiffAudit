"""
Step 5: secondary metrics on the docking poses already produced by Step 3.

5b. PoseBusters validity pass rate per condition.
5c. Per-ligand winning conformer (which beta value produces the best top-1
    pose in condition A) -- one winner throughout vs different winners per
    ligand distinguishes "ensemble does real work" from "one better receptor."

5a (ProLIF IFP recovery) was dropped -- see lasb_step5_driver.py docstring:
ProLIF/MDAnalysis's RDKit bond-order standardization segfaults on the H-less
receptor PDBs used throughout this experiment.
"""
import csv
from collections import defaultdict
from posebusters import PoseBusters


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

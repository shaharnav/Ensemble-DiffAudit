"""
Pocket-lining conformational displacement (Phase 5).

If the ensemble effect is real, the per-ligand gain should correlate with how much the
winning conformer's pocket actually moved. If it doesn't, the "gain" is noise dressed up
as induced fit. This computes, for every conformer already aligned to the crystal by
ensemble_auditor.py, the CA-only and all-heavy-atom RMSD over pocket-lining residues, plus
the single largest CA displacement -- so that covariate can be checked against
delta_ensemble_vs_crystal per ligand.

Pocket-lining residues are defined as those with any heavy atom within --cutoff (default
8.0 A) of the co-crystallized ligand's centroid in the crystal structure.

Usage:
    ./venv/bin/python pocket_rmsd.py
"""
import argparse
import csv
import glob
import json
import logging
import os
import sys

import numpy as np
from Bio.PDB import PDBParser

from calibrate import get_ben_centroid

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

TRYPSIN_PDB = "pdbs/3PTB.pdb"
LIGAND_CODE = "BEN"
ALIGNED_RECEPTORS_DIR = os.path.join("results", "payload_unpacked", "aligned_receptors")
OUTPUT_CSV = "conformer_geometry.csv"
RESULTS_JSON = "results.json"
DEFAULT_CUTOFF = 8.0


def find_pocket_lining_residues(structure, ligand_centroid: np.ndarray, cutoff: float) -> list:
    """Standard (non-HETATM) residues with any heavy atom within *cutoff* A of the
    ligand centroid, ordered by sequence position."""
    lining = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " " or not residue.has_id("CA"):
                    continue
                for atom in residue:
                    if atom.element == "H":
                        continue
                    if np.linalg.norm(atom.get_coord() - ligand_centroid) <= cutoff:
                        lining.append(residue)
                        break
    return lining


def _standard_ca_residues(structure):
    return [
        res
        for model in structure
        for chain in model
        for res in chain
        if res.id[0] == " " and res.has_id("CA")
    ]


def compute_conformer_geometry(ref_structure, conformer_path: str, lining_indices: set[int]) -> dict:
    """
    ref_structure and the conformer are both already aligned into the same frame
    (ensemble_auditor.py superimposed conformers onto the crystal by sequence-index CA
    pairing before writing them to ALIGNED_RECEPTORS_DIR) -- so raw coordinate distances
    directly reflect conformational displacement, no further superposition needed here.
    """
    parser = PDBParser(QUIET=True)
    conf_structure = parser.get_structure("conf", conformer_path)

    ref_residues = _standard_ca_residues(ref_structure)
    conf_residues = _standard_ca_residues(conf_structure)

    if len(ref_residues) != len(conf_residues):
        raise ValueError(
            f"{conformer_path}: residue count mismatch vs. reference "
            f"({len(conf_residues)} vs {len(ref_residues)})"
        )

    ca_deltas = []
    allatom_deltas = []
    max_ca_displacement = 0.0

    for idx, (ref_res, conf_res) in enumerate(zip(ref_residues, conf_residues)):
        if idx not in lining_indices:
            continue
        if ref_res.get_resname() != conf_res.get_resname():
            raise ValueError(
                f"{conformer_path}: residue-name mismatch at sequence position {idx} "
                f"({ref_res.get_resname()} vs {conf_res.get_resname()})"
            )

        ca_disp = np.linalg.norm(ref_res["CA"].get_coord() - conf_res["CA"].get_coord())
        ca_deltas.append(ca_disp)
        max_ca_displacement = max(max_ca_displacement, ca_disp)

        conf_atoms_by_name = {a.get_name(): a for a in conf_res if a.element != "H"}
        for ref_atom in ref_res:
            if ref_atom.element == "H":
                continue
            conf_atom = conf_atoms_by_name.get(ref_atom.get_name())
            if conf_atom is None:
                continue
            allatom_deltas.append(
                np.linalg.norm(ref_atom.get_coord() - conf_atom.get_coord())
            )

    pocket_ca_rmsd = float(np.sqrt(np.mean(np.square(ca_deltas)))) if ca_deltas else None
    pocket_allatom_rmsd = (
        float(np.sqrt(np.mean(np.square(allatom_deltas)))) if allatom_deltas else None
    )

    return {
        "pocket_ca_rmsd": pocket_ca_rmsd,
        "pocket_allatom_rmsd": pocket_allatom_rmsd,
        "max_ca_displacement": max_ca_displacement if ca_deltas else None,
        "n_pocket_residues": len(ca_deltas),
    }


def join_geometry_to_results(results: list[dict], geometry_by_conformer: dict[str, dict]) -> list[dict]:
    """Attach each ligand's ensemble_best_conformer geometry onto its results row."""
    joined = []
    for r in results:
        winner = r.get("ensemble_best_conformer")
        geometry = geometry_by_conformer.get(winner, {})
        joined.append({
            **r,
            "winner_pocket_ca_rmsd": geometry.get("pocket_ca_rmsd"),
            "winner_pocket_allatom_rmsd": geometry.get("pocket_allatom_rmsd"),
            "winner_max_ca_displacement": geometry.get("max_ca_displacement"),
        })
    return joined


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF,
                         help=f"Pocket-lining distance cutoff in A (default {DEFAULT_CUTOFF})")
    args = parser.parse_args()

    if not os.path.exists(TRYPSIN_PDB):
        logger.error(f"{TRYPSIN_PDB} not found.")
        return 1

    bio_parser = PDBParser(QUIET=True)
    ref_structure = bio_parser.get_structure("ref", TRYPSIN_PDB)
    ligand_centroid = get_ben_centroid(TRYPSIN_PDB)
    logger.info(f"Ligand centroid: {ligand_centroid}")

    lining_residues = find_pocket_lining_residues(ref_structure, ligand_centroid, args.cutoff)
    ref_residues = _standard_ca_residues(ref_structure)
    ref_res_ids = {id(r): i for i, r in enumerate(ref_residues)}
    lining_indices = {ref_res_ids[id(r)] for r in lining_residues}
    logger.info(
        f"{len(lining_indices)} pocket-lining residues within {args.cutoff} A of ligand centroid"
    )

    conformer_paths = sorted(glob.glob(os.path.join(ALIGNED_RECEPTORS_DIR, "conformix_var_*.pdb")))
    if not conformer_paths:
        logger.error(f"No conformers found in {ALIGNED_RECEPTORS_DIR}. Run ensemble_auditor.py first.")
        return 1

    rows = []
    for conf_path in conformer_paths:
        basename = os.path.basename(conf_path)
        try:
            geometry = compute_conformer_geometry(ref_structure, conf_path, lining_indices)
        except ValueError as e:
            logger.warning(f"  ⚠ {e}")
            continue
        rows.append({"conformer": basename, **geometry})
        logger.info(
            f"  {basename}: pocket_ca_rmsd={geometry['pocket_ca_rmsd']:.3f} A, "
            f"pocket_allatom_rmsd={geometry['pocket_allatom_rmsd']:.3f} A, "
            f"max_ca_displacement={geometry['max_ca_displacement']:.3f} A"
        )

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "conformer", "pocket_ca_rmsd", "pocket_allatom_rmsd",
            "max_ca_displacement", "n_pocket_residues",
        ])
        writer.writeheader()
        writer.writerows(rows)
    logger.info(f"Conformer geometry written -> {OUTPUT_CSV}")

    if os.path.exists(RESULTS_JSON):
        with open(RESULTS_JSON) as f:
            results = json.load(f)
        geometry_by_conformer = {row["conformer"]: row for row in rows}
        joined = join_geometry_to_results(results, geometry_by_conformer)
        logger.info("Per-ligand winning-conformer geometry (for Phase 8's correlation check):")
        for r in joined:
            logger.info(
                f"  {r.get('id', '?')}: delta_ensemble_vs_crystal="
                f"{r.get('delta_ensemble_vs_crystal')}, "
                f"winner_pocket_ca_rmsd={r.get('winner_pocket_ca_rmsd')}"
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())

"""
Positive-control calibration for the docking engine.

Redocks the native benzamidine (BEN) ligand of 3PTB into its own waters/BEN-stripped
receptor and checks that the top-scoring pose recovers the crystallographic pose
(symmetry-corrected heavy-atom RMSD), across a fixed set of Vina seeds. Pose recovery
is the field-standard positive control -- a docking engine that reproduces the correct
binding score for the wrong pose is not actually validated.

Usage:
    ./venv/bin/python calibrate.py
"""
import json
import logging
import os
import statistics
import sys

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from docking_engine import run_docking

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

TRYPSIN_PDB = "pdbs/3PTB.pdb"
BENZAMIDINE_SMILES = "NC(=N)c1ccccc1"
LIGAND_CODE = "BEN"
WORK_DIR = "results/calibration"
STRIPPED_RECEPTOR = os.path.join(WORK_DIR, "3PTB_stripped.pdb")
REPORT_PATH = "validation/redocking_report.json"

# Same seed list Phase 4's rigid control will use, so calibration and the rigid
# control are drawn from an identical sampling regime.
SEEDS = [1, 2, 3, 4, 5, 6]
EXHAUSTIVENESS = 16
BOX_SIZE = [22.5, 22.5, 22.5]
RMSD_PASS_THRESHOLD = 2.0


class KeepProteinAndMetals(Select):
    """Drop crystallographic waters and the native BEN ligand; keep the protein and Ca2+."""

    def accept_residue(self, residue):
        resname = residue.get_resname().strip().upper()
        if resname in ("HOH", "WAT", "BEN"):
            return False
        return True


def get_ben_centroid(pdb_file: str) -> np.ndarray:
    """Read the native benzamidine (BEN) HETATM coordinates and return their centroid."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("3ptb", pdb_file)
    coords = []
    for atom in structure.get_atoms():
        residue = atom.get_parent()
        if residue.get_resname().strip().upper() == "BEN":
            coords.append(atom.get_coord())
    if not coords:
        raise RuntimeError(f"No BEN (benzamidine) HETATM records found in {pdb_file}")
    return np.mean(np.array(coords), axis=0)


def strip_receptor(pdb_file: str, output_path: str) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("3ptb", pdb_file)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, KeepProteinAndMetals())


def extract_reference_ligand(pdb_file: str, resname: str, template_smiles: str) -> Chem.Mol:
    """
    Build an RDKit mol (heavy atoms only, correct bond orders) for the co-crystallized
    ligand from its HETATM coordinates, since the PDB itself carries no bond information.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ref", pdb_file)

    lines = []
    serial = 1
    for atom in structure.get_atoms():
        residue = atom.get_parent()
        if residue.get_resname().strip().upper() != resname:
            continue
        x, y, z = atom.get_coord()
        name = atom.get_name()
        elem = (atom.element or name[0]).strip()
        lines.append(
            f"HETATM{serial:5d} {name:<4s} {resname} A 999    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}"
        )
        serial += 1
    if not lines:
        raise RuntimeError(f"No {resname} HETATM records found in {pdb_file}")
    block = "\n".join(lines) + "\nEND\n"

    raw_mol = Chem.MolFromPDBBlock(block, removeHs=False)
    if raw_mol is None:
        raise RuntimeError(f"RDKit failed to parse extracted {resname} HETATM block")

    template = Chem.MolFromSmiles(template_smiles)
    mol = AllChem.AssignBondOrdersFromTemplate(template, raw_mol)
    return Chem.RemoveHs(mol)


def redock_top_pose_rmsd(docked_pdbqt: str, reference_mol: Chem.Mol) -> float:
    """Symmetry-corrected heavy-atom RMSD between the top-scoring docked pose and the
    reference crystal pose. Both molecules live in the receptor's coordinate frame
    already -- no additional superposition is performed, since we're measuring
    positional accuracy in the binding site, not shape similarity."""
    pdbqt_mol = PDBQTMolecule.from_file(docked_pdbqt, skip_typing=True)
    docked_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    docked = Chem.RemoveHs(docked_mols[0])
    # Conformer 0 is Vina's top-scoring (first) binding mode.
    return rdMolAlign.CalcRMS(docked, reference_mol, prbId=0, refId=0)


def main() -> int:
    if not os.path.exists(TRYPSIN_PDB):
        logger.error(
            f"{TRYPSIN_PDB} not found. Download it with:\n"
            f"  curl -L https://files.rcsb.org/download/3PTB.pdb -o {TRYPSIN_PDB}"
        )
        return 1

    os.makedirs(WORK_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)

    center = get_ben_centroid(TRYPSIN_PDB)
    logger.info(f"Native benzamidine centroid (true S1 pocket center): {center}")

    strip_receptor(TRYPSIN_PDB, STRIPPED_RECEPTOR)
    logger.info(f"Stripped waters + native ligand -> {STRIPPED_RECEPTOR}")

    reference_mol = extract_reference_ligand(TRYPSIN_PDB, LIGAND_CODE, BENZAMIDINE_SMILES)
    logger.info(f"Reference crystal pose reconstructed ({reference_mol.GetNumAtoms()} heavy atoms).")

    rmsd_per_seed: dict[int, float] = {}
    score_per_seed: dict[int, float] = {}

    for seed in SEEDS:
        job_name = f"calibration_benzamidine_seed{seed}"
        result = run_docking(
            pdb_file=STRIPPED_RECEPTOR,
            smiles=BENZAMIDINE_SMILES,
            output_dir=WORK_DIR,
            job_name=job_name,
            exhaustiveness=EXHAUSTIVENESS,
            center_coords=center.tolist(),
            box_size=BOX_SIZE,
            seed=seed,
        )

        if result is None or result.get("affinity") is None:
            logger.error(f"Docking failed for seed {seed}.")
            continue

        affinity = result["affinity"]
        docked_pdbqt = os.path.join(WORK_DIR, f"{job_name}_out.pdbqt")
        rmsd = redock_top_pose_rmsd(docked_pdbqt, reference_mol)

        rmsd_per_seed[seed] = rmsd
        score_per_seed[seed] = affinity
        logger.info(f"  seed={seed}: affinity={affinity:.3f} kcal/mol, top-pose RMSD={rmsd:.3f} A")

    if not rmsd_per_seed:
        logger.error("All redocking runs failed -- cannot calibrate.")
        return 1

    rmsd_values = list(rmsd_per_seed.values())
    median_rmsd = statistics.median(rmsd_values)
    pass_count = sum(1 for r in rmsd_values if r <= RMSD_PASS_THRESHOLD)
    pass_rate = pass_count / len(rmsd_values)

    report = {
        "target": "3PTB",
        "ligand_code": LIGAND_CODE,
        "rmsd_per_seed": rmsd_per_seed,
        "median_rmsd": median_rmsd,
        "score_per_seed": score_per_seed,
        "pass_rate_2A": pass_rate,
        "vina_version": "1.2.7",
        "exhaustiveness": EXHAUSTIVENESS,
    }
    with open(REPORT_PATH, "w") as fh:
        json.dump(report, fh, indent=2)
    logger.info(f"Redocking report written -> {REPORT_PATH}")

    logger.info(
        f"Median top-pose RMSD across {len(rmsd_values)} seeds: {median_rmsd:.3f} A "
        f"(pass rate @ {RMSD_PASS_THRESHOLD} A: {pass_rate:.0%})"
    )

    if median_rmsd <= RMSD_PASS_THRESHOLD:
        logger.info(f"PASS: median RMSD {median_rmsd:.3f} A is within {RMSD_PASS_THRESHOLD} A.")
        return 0
    else:
        logger.error(
            f"FAIL: median RMSD {median_rmsd:.3f} A exceeds {RMSD_PASS_THRESHOLD} A. "
            f"The docking engine is not reliably recovering the correct pose -- "
            f"investigate before trusting any candidate affinity."
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())

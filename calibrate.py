"""
Positive-control calibration for the docking engine.

Docks benzamidine (the native S1-pocket ligand of 3PTB) into the real, waters/BEN-stripped
3PTB crystal structure and checks that the score lands in the experimentally reasonable
range. If this fails, the engine is wrong -- not the generated molecules.

Usage:
    ./venv/bin/python calibrate.py
"""
import logging
import os
import sys

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

from docking_engine import run_docking

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

TRYPSIN_PDB = "pdbs/3PTB.pdb"
BENZAMIDINE_SMILES = "NC(=N)c1ccccc1"
WORK_DIR = "results/calibration"
STRIPPED_RECEPTOR = os.path.join(WORK_DIR, "3PTB_stripped.pdb")

EXPECTED_LOW, EXPECTED_HIGH = -7.0, -6.0
BOX_SIZE = [22.5, 22.5, 22.5]


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


def main() -> int:
    if not os.path.exists(TRYPSIN_PDB):
        logger.error(
            f"{TRYPSIN_PDB} not found. Download it with:\n"
            f"  curl -L https://files.rcsb.org/download/3PTB.pdb -o {TRYPSIN_PDB}"
        )
        return 1

    os.makedirs(WORK_DIR, exist_ok=True)

    center = get_ben_centroid(TRYPSIN_PDB)
    logger.info(f"Native benzamidine centroid (true S1 pocket center): {center}")

    strip_receptor(TRYPSIN_PDB, STRIPPED_RECEPTOR)
    logger.info(f"Stripped waters + native ligand -> {STRIPPED_RECEPTOR}")

    result = run_docking(
        pdb_file=STRIPPED_RECEPTOR,
        smiles=BENZAMIDINE_SMILES,
        output_dir=WORK_DIR,
        job_name="calibration_benzamidine",
        exhaustiveness=16,
        center_coords=center.tolist(),
        box_size=BOX_SIZE,
    )

    if result is None or result.get("affinity") is None:
        logger.error("Calibration docking run failed to produce an affinity.")
        return 1

    affinity = result["affinity"]
    logger.info(f"Benzamidine affinity vs. real 3PTB: {affinity:.3f} kcal/mol")

    if EXPECTED_LOW <= affinity <= EXPECTED_HIGH:
        logger.info(
            f"PASS: {affinity:.3f} kcal/mol is within the expected "
            f"[{EXPECTED_LOW}, {EXPECTED_HIGH}] range."
        )
        return 0
    else:
        logger.error(
            f"FAIL: {affinity:.3f} kcal/mol is outside the expected "
            f"[{EXPECTED_LOW}, {EXPECTED_HIGH}] range. The docking engine is miscalibrated "
            f"-- investigate before trusting any candidate affinity."
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())

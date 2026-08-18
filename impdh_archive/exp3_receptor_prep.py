"""
Experiment 3, Phase 1c: receptor preparation for CLas IMPDH (PDB 6KCF).

Builds R_apo: chain A only, waters stripped, protein atoms only.

Single-chain convention: IMPDH is a homotetramer, but the catalytic site (IMP- and
NAD-binding subsites) is intra-subunit in every characterized IMPDH structure -- the
CBS/Bateman domains (deleted in this construct) mediate tetramer contacts and allosteric
regulation, not catalysis. verify_6kcf.py measured the catalytic Cys303 SG to be 8.09 A
from the nearest atom of the next-closest chain (B) -- too far to be pocket-lining at the
8 A cutoff used throughout this project, and consistent with an intra-subunit active site
rather than an interface one. Single-chain R_apo follows the convention already used for
Experiments 1/2 and Xylellain.

Altloc handling: verify_6kcf.py's structure scan found zero disordered (altloc) residues in
6KCF, so no altloc-resolution logic is exercised here -- the Select subclass below still
implements deterministic altloc handling (higher occupancy, tiebreak 'A') defensively, per
the project's standing rule, in case a future re-run against a re-refined deposition differs.

Usage:
    ./venv/bin/python exp3_receptor_prep.py
"""
import logging
import os
import sys

from Bio.PDB import PDBParser, PDBIO, Select

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

RAW_PDB = os.path.join("pdbs", "exp3_impdh", "6KCF.pdb")
CHAIN_ID = "A"
OUT_DIR = os.path.join("results", "experiment3", "receptors")
R_APO_PATH = os.path.join(OUT_DIR, "R_apo.pdb")


class DeterministicAltlocSelect(Select):
    """Chain A, standard (ATOM) residues only -- waters and any heteroatom dropped. Altlocs
    resolved deterministically: highest occupancy first, tiebreak alphabetically (so 'A' wins
    over 'B' on an exact tie), matching the rule used for Xylellain (3OIS)."""

    def __init__(self):
        self._chosen_altloc = {}

    def accept_chain(self, chain):
        return chain.id == CHAIN_ID

    def accept_residue(self, residue):
        return residue.id[0] == " "

    def accept_atom(self, atom):
        if not atom.is_disordered():
            return True
        residue = atom.get_parent()
        key = (residue.get_parent().id, residue.id, atom.get_name())
        if key not in self._chosen_altloc:
            candidates = sorted(
                atom.disordered_get_id_list(),
                key=lambda loc: (-atom.disordered_get(loc).get_occupancy(), loc),
            )
            self._chosen_altloc[key] = candidates[0]
        chosen = self._chosen_altloc[key]
        if atom.get_altloc() != chosen:
            return False
        atom.set_altloc(" ")
        atom.set_occupancy(1.0)
        return True


def main() -> int:
    if not os.path.isfile(RAW_PDB):
        logger.error(f"{RAW_PDB} not found.")
        return 1
    os.makedirs(OUT_DIR, exist_ok=True)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("6kcf", RAW_PDB)

    io = PDBIO()
    io.set_structure(structure)
    io.save(R_APO_PATH, DeterministicAltlocSelect())

    check = PDBParser(QUIET=True).get_structure("r_apo", R_APO_PATH)
    residues = [r for r in check[0][CHAIN_ID] if r.id[0] == " "]
    logger.info(f"R_apo written -> {R_APO_PATH}: chain {CHAIN_ID}, {len(residues)} residues, "
                f"range {residues[0].id[1]}-{residues[-1].id[1]}")
    het = [r for r in check[0].get_residues() if r.id[0] != " "]
    if het:
        logger.error(f"R_apo still contains {len(het)} heteroatom residues -- selection logic failed.")
        return 1
    logger.info("No heteroatoms in R_apo (waters/ligands correctly stripped).")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""
Experiment 3, Phase 1: receptor preparation for the Xylellain (3OIS) campaign.

Builds two clean receptor PDBs from chain A of the deposited zymogen structure, no waters,
no HETATM (the UDP ribonucleotide is peripheral -- see verify_3ois.py):

  R1 -- zymogen (occluded): residues 23-290 as deposited. The only experimentally observed
        state. Floor.
  R2 -- propeptide removed (modeled ceiling proxy): residues 67-290 only. Not remodeled,
        minimized, or relaxed -- an unrelaxed deletion, per the plan, since that is a more
        honest proxy than an energy-minimized structure that introduces its own artifacts.
        No explicit N-terminal cap is added; residue 67 (Phe) is left as a normal residue,
        exactly as every other receptor's terminus is treated by this pipeline's existing
        receptor-prep step (Meeko/mk_prepare_receptor in docking_engine.py).

Propeptide boundary (residues 23-66) is not assumed from the abstract -- it's the contiguous
N-terminal segment ending where the conserved mature-domain motif immediately preceding the
catalytic cysteine (residues 67-78, "...QGRIGSC...") begins folding into the active-site
cleft. See verify_3ois.py's output for the distance evidence: residues 26-36 are the
occluding loop proper (3.9-9.8 A from the catalytic-triad centroid); residues 37-66 sit far
from the site (12-31 A) but are still N-terminal to, and contiguous with, that occluding
loop and the conserved motif -- consistent with the paper's own description ("part of the
N-terminal sequence blocks the active site," Leite et al. 2013 abstract) rather than the
full propeptide occluding uniformly.

Usage:
    ./venv/bin/python exp3_receptor_prep.py
"""
import logging
import os
import sys

from Bio.PDB import PDBParser, PDBIO, Select

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

SOURCE_PDB = "pdbs/3OIS.pdb"
CHAIN_ID = "A"
PROPEPTIDE_RANGE = (23, 66)  # inclusive
OUT_DIR = os.path.join("results", "experiment3", "receptors")
R1_PATH = os.path.join(OUT_DIR, "R1_zymogen.pdb")
R2_PATH = os.path.join(OUT_DIR, "R2_mature_proxy.pdb")


class ChainAProteinOnly(Select):
    """Chain A, standard residues only -- drops waters, UDP, and the other 3 crystallographic copies."""

    def accept_residue(self, residue):
        return residue.get_parent().id == CHAIN_ID and residue.id[0] == " "

    def accept_chain(self, chain):
        return chain.id == CHAIN_ID


class ChainAMatureOnly(ChainAProteinOnly):
    """Chain A, standard residues, propeptide (residues 23-66) removed."""

    def accept_residue(self, residue):
        if not super().accept_residue(residue):
            return False
        return not (PROPEPTIDE_RANGE[0] <= residue.id[1] <= PROPEPTIDE_RANGE[1])


def main() -> int:
    if not os.path.exists(SOURCE_PDB):
        logger.error(f"{SOURCE_PDB} not found. Fetch via "
                      f"curl -L https://files.rcsb.org/download/3OIS.pdb -o {SOURCE_PDB}")
        return 1

    os.makedirs(OUT_DIR, exist_ok=True)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("3ois", SOURCE_PDB)

    io = PDBIO()
    io.set_structure(structure)
    io.save(R1_PATH, ChainAProteinOnly())
    logger.info(f"R1 (zymogen, chain A, no waters/HETATM) -> {R1_PATH}")

    io.save(R2_PATH, ChainAMatureOnly())
    logger.info(f"R2 (propeptide residues {PROPEPTIDE_RANGE[0]}-{PROPEPTIDE_RANGE[1]} "
                f"removed) -> {R2_PATH}")

    for path, label in [(R1_PATH, "R1"), (R2_PATH, "R2")]:
        p2 = PDBParser(QUIET=True)
        s2 = p2.get_structure(label, path)
        residues = [r for r in s2[0]["A"] if r.has_id("CA")]
        logger.info(f"  {label}: {len(residues)} residues, "
                    f"{residues[0].id[1]}-{residues[-1].id[1]}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

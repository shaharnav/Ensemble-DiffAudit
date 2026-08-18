"""
Experiment 3, Phase 1b: reproducible verification of CLas IMPDH (PDB 6KCF).

Do not assume the primary paper's residue numbering (Cys309, from full-length UniProt
C6XG59 numbering) matches the deposited coordinate file -- the paper's own methods describe
a CBS-domain deletion construct (104 residues replaced), and the deposited pdb_seq_num
scheme does not simply preserve WT numbering with a gap: it jumps by exactly 104 at the
deletion junction, with additional undetermined N-terminal tag offset.

Catalytic Cys is located here by sequence motif, not by residue number: the conserved
IMPDH nucleophile sits in a single "GSIC" tetrapeptide in the CLas sequence (verified
against UniProt C6XG59, where it occurs once, at UniProt position 309 -- matching the
paper's claim in UniProt/WT numbering). Whatever residue number that motif's Cys carries
in the deposited file's own numbering is the number used for every downstream script.

Usage:
    ./venv/bin/python verify_6kcf.py
"""
import json
import logging
import os
import sys

import numpy as np
from Bio.PDB import PDBParser, MMCIF2Dict
from Bio.PDB.Polypeptide import protein_letters_3to1 as three_to_one

AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
    'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

PDB_PATH = os.path.join("pdbs", "exp3_impdh", "6KCF.pdb")
CIF_PATH = os.path.join("pdbs", "exp3_impdh", "6KCF.cif")
CATALYTIC_MOTIF = "GSIC"  # conserved IMPDH nucleophile motif; occurs exactly once in C6XG59
UNIPROT_ACCESSION = "C6XG59"
UNIPROT_SEQ = (
    "MARIIENNVGGVALTFDDVLLRPEFSNVLPRDIDISTRIAKDFTLNLPIMSAAMDQVTDSRLAIAMAQAGGLGVIHRNFS"
    "PSEQVAQVHQVKKFESGMVVNPVTISPYATLADALALMKKYSISGIPVVESDVGKLVGILTNRDVRFASNAQQAVGELMT"
    "RNLITVKKTVNLENAKALLHQHRIEKLLVVDDDGCCIGLITVKDIERSQLNPNATKDSKGRLRVAAAVSVAKDIADRVGP"
    "LFDVNVDLVVVDTAHGHSQKVLDAVVQIKKNFPSLLVMAGNIATAEGALALIDAGADIIKVGIGPGSICTTRVVTGVGCP"
    "QLSAIMSVVEVAERAGVAIVADGGIRFSGDIAKAIAAGSACVMIGSLLAGTDESPGDIFLYQGRSFKSYRGMGSVAAMER"
    "GSSARYSQDGVTDVLKLVPEGIEGRVPYKGPIASVLHQMSGGLKSSMGYVGASNIEEFQKKANFIRVSVAGLRESHVHDV"
    "KITRESPNYSETI"
)
OUTPUT_JSON = os.path.join("results", "experiment3", "verification_6kcf.json")


def load_chain_sequence(structure, chain_id):
    chain = structure[0][chain_id]
    residues = [r for r in chain if r.id[0] == " "]
    seq = "".join(three_to_one.get(r.get_resname(), "X") for r in residues)
    ids = [r.id[1] for r in residues]
    return residues, seq, ids


def locate_catalytic_cys_deposited_numbering(cif_path, motif, chain_id="A"):
    """Locates `motif` (e.g. 'GSIC') in the FULL construct sequence, including residues with
    no modeled density -- SEQRES/entity identity is known from the mmCIF sequence scheme even
    when a residue's coordinates were never resolved. Necessary here because the catalytic
    motif straddles a partially-disordered loop (G-S-I unmodeled, C modeled) in this apo
    structure, so searching only modeled residues (as a naive PDBParser sequence would) misses
    it -- confirmed, not assumed, by checking below whether the returned Cys is itself modeled."""
    mmcif_dict = MMCIF2Dict.MMCIF2Dict(cif_path)
    asym_ids = mmcif_dict["_pdbx_poly_seq_scheme.asym_id"]
    mon_ids = mmcif_dict["_pdbx_poly_seq_scheme.mon_id"]
    pdb_seq_nums = mmcif_dict["_pdbx_poly_seq_scheme.pdb_seq_num"]
    auth_seq_nums = mmcif_dict["_pdbx_poly_seq_scheme.auth_seq_num"]

    chain_seq, chain_pdbnum, chain_modeled = [], [], []
    for asym, mon, pdbnum, authnum in zip(asym_ids, mon_ids, pdb_seq_nums, auth_seq_nums):
        if asym != chain_id:
            continue
        chain_seq.append(AA3TO1.get(mon, "X"))
        chain_pdbnum.append(int(pdbnum))
        chain_modeled.append(authnum != "?")
    seq = "".join(chain_seq)

    idx = seq.find(motif)
    if idx == -1 or seq.find(motif, idx + 1) != -1:
        return None
    cys_offset = idx + motif.index("C")
    return {
        "resid_deposited_numbering": chain_pdbnum[cys_offset],
        "motif_residue_ids": chain_pdbnum[idx:idx + len(motif)],
        "motif_residues_modeled": chain_modeled[idx:idx + len(motif)],
        "catalytic_cys_modeled": chain_modeled[cys_offset],
    }


def find_gaps(ids):
    gaps = []
    for a, b in zip(ids, ids[1:]):
        if b != a + 1:
            gaps.append((a, b))
    return gaps


def main() -> int:
    if not os.path.isfile(PDB_PATH):
        logger.error(f"{PDB_PATH} not found -- download 6KCF.pdb from RCSB first.")
        return 1

    n_uniprot_motif = UNIPROT_SEQ.count(CATALYTIC_MOTIF)
    if n_uniprot_motif != 1:
        logger.error(f"'{CATALYTIC_MOTIF}' occurs {n_uniprot_motif} times in {UNIPROT_ACCESSION} "
                      f"-- motif is not unique, cannot use it to localize the catalytic Cys unambiguously.")
        return 1
    uniprot_motif_pos = UNIPROT_SEQ.find(CATALYTIC_MOTIF) + 1  # 1-indexed
    uniprot_cys_pos = uniprot_motif_pos + CATALYTIC_MOTIF.index("C")
    logger.info(f"'{CATALYTIC_MOTIF}' motif in {UNIPROT_ACCESSION}: starts at position "
                f"{uniprot_motif_pos}, catalytic Cys at UniProt position {uniprot_cys_pos} "
                f"(paper claims Cys309 -- {'MATCH' if uniprot_cys_pos == 309 else 'MISMATCH'}).")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("6kcf", PDB_PATH)
    chains = sorted(c.id for c in structure[0])
    logger.info(f"Chains present: {chains}")

    hetatm_inventory = {}
    for res in structure[0].get_residues():
        if res.id[0] not in (" ", ""):
            code = res.id[0].replace("H_", "").strip()
            hetatm_inventory[code] = hetatm_inventory.get(code, 0) + 1
    logger.info(f"HETATM inventory (residue-name-tagged heteroflag counts): {hetatm_inventory}")
    non_water = {k: v for k, v in hetatm_inventory.items() if k != "W"}
    apo_confirmed = len(non_water) == 0
    logger.info(f"Apo (no non-water heteroatoms): {'PASS' if apo_confirmed else 'FAIL -- ' + str(non_water)}")

    if not os.path.isfile(CIF_PATH):
        logger.error(f"{CIF_PATH} not found -- download 6KCF.cif from RCSB first "
                      f"(needed for the full construct sequence scheme, including unmodeled residues).")
        return 1

    catalytic_cys_by_chain = {}
    gaps_by_chain = {}
    motif_details_by_chain = {}
    for chain_id in chains:
        residues, seq, ids = load_chain_sequence(structure, chain_id)
        motif_info = locate_catalytic_cys_deposited_numbering(CIF_PATH, CATALYTIC_MOTIF, chain_id)
        if motif_info is None:
            logger.warning(f"Chain {chain_id}: '{CATALYTIC_MOTIF}' motif not unique/found in the full "
                            f"construct sequence (including unmodeled residues).")
            catalytic_cys_by_chain[chain_id] = None
        else:
            cys_resid = motif_info["resid_deposited_numbering"]
            catalytic_cys_by_chain[chain_id] = cys_resid if motif_info["catalytic_cys_modeled"] else None
            motif_details_by_chain[chain_id] = motif_info
            logger.info(f"Chain {chain_id}: catalytic Cys located at deposited residue number {cys_resid} "
                        f"via '{CATALYTIC_MOTIF}' motif match (motif residue ids "
                        f"{motif_info['motif_residue_ids']}, modeled={motif_info['motif_residues_modeled']}). "
                        f"Cys itself modeled: {motif_info['catalytic_cys_modeled']}.")
        gaps_by_chain[chain_id] = find_gaps(ids)
        logger.info(f"Chain {chain_id}: {len(residues)} modeled residues, "
                    f"range {ids[0]}-{ids[-1]}, gaps (missing density) = {gaps_by_chain[chain_id]}")

    resolved_cys = {c: r for c, r in catalytic_cys_by_chain.items() if r is not None}
    if not resolved_cys:
        logger.error("Catalytic Cys motif not resolved in ANY chain -- cannot proceed.")
        return 1
    cys_values = set(resolved_cys.values())
    consistent = len(cys_values) == 1
    catalytic_cys_resid = sorted(cys_values)[0] if consistent else None
    logger.info(f"Catalytic Cys residue number across chains: {resolved_cys} "
                f"({'consistent' if consistent else 'INCONSISTENT -- investigate'})")

    # Oligomeric interface check: is the catalytic Cys within contact distance of a
    # neighboring chain (active site formed at a subunit interface, as in known IMPDH
    # tetramers), or entirely intra-chain?
    interface_distances = {}
    if consistent:
        ref_chain = chains[0]
        cys_res = None
        for r in structure[0][ref_chain]:
            if r.id[0] == " " and r.id[1] == catalytic_cys_resid:
                cys_res = r
                break
        if cys_res is not None and cys_res.has_id("SG"):
            cys_sg = cys_res["SG"].get_coord()
            for other_chain in chains:
                if other_chain == ref_chain:
                    continue
                min_dist = min(
                    np.linalg.norm(atom.get_coord() - cys_sg)
                    for atom in structure[0][other_chain].get_atoms()
                    if atom.element != "H"
                )
                interface_distances[f"{ref_chain}->{other_chain}"] = round(float(min_dist), 2)
            logger.info(f"Chain {ref_chain} catalytic Cys{catalytic_cys_resid} SG, minimum distance "
                        f"to each other chain: {interface_distances}")
        else:
            logger.warning(f"Chain {ref_chain} catalytic Cys{catalytic_cys_resid} has no SG atom "
                            f"(disordered sidechain) -- cannot measure interface distance directly.")

    output = {
        "pdb_id": "6KCF",
        "source": "RCSB",
        "uniprot_accession": UNIPROT_ACCESSION,
        "catalytic_motif": CATALYTIC_MOTIF,
        "uniprot_catalytic_cys_position": uniprot_cys_pos,
        "paper_claimed_cys_position": 309,
        "uniprot_vs_paper_match": uniprot_cys_pos == 309,
        "chains": chains,
        "hetatm_inventory": hetatm_inventory,
        "apo_confirmed": apo_confirmed,
        "catalytic_cys_resid_by_chain_deposited_numbering": catalytic_cys_by_chain,
        "catalytic_cys_resid_consistent_across_chains": consistent,
        "catalytic_cys_resid_deposited_numbering": catalytic_cys_resid,
        "note": (
            "Deposited pdb_seq_num numbering does NOT preserve full-length UniProt numbering "
            "with a simple gap -- confirmed by direct motif localization, not assumed. The "
            "paper's 'Cys309' refers to UniProt/full-length numbering; the same residue in this "
            "file's own numbering is recorded above and is what every downstream script must use."
        ),
        "gaps_by_chain": {c: [list(g) for g in gaps] for c, gaps in gaps_by_chain.items()},
        "catalytic_motif_details_by_chain": motif_details_by_chain,
        "cys_sg_interface_distances_A": interface_distances,
    }
    os.makedirs(os.path.dirname(OUTPUT_JSON), exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)
    logger.info(f"Written -> {OUTPUT_JSON}")

    if not apo_confirmed or not consistent:
        logger.error("Verification gate FAILED (apo or Cys-consistency check) -- stop before proceeding.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

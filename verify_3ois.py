"""
Experiment 3, Phase 1a: structural verification of 3OIS (Xylellain).

No assumptions from the abstract -- everything here is computed from the deposited
coordinates:

1. Enumerate HETATM records (chains, residue names).
2. Locate the ribonucleotide (UDP) and compute its distance to the catalytic-triad centroid.
   Gate: if within 8 A, the site is not apo and the target premise fails.
3. Identify the catalytic triad by sequence motif (papain-family Cys in a "GxCxxx" context)
   plus 3D proximity (His within H-bonding/ion-pair distance of the Cys thiol, Asn within
   H-bonding distance of that His) -- not by assuming papain's own numbering transfers.
4. Characterize the N-terminal occluding segment: per-residue minimum distance from the
   catalytic-triad centroid, for residues N-terminal to the catalytic Cys, to find the
   actual occluding loop rather than assume a residue count from the abstract.

Usage:
    ./venv/bin/python verify_3ois.py
"""
import json
import logging
import sys

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as three_to_one

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

SOURCE_PDB = "pdbs/3OIS.pdb"
CHAIN_ID = "A"
CATALYTIC_MOTIF = "GSCTANA"  # papain-family Gly-x-Cys lead-in, located by sequence scan
ACTIVE_SITE_GATE_ANGSTROM = 8.0
OUTPUT_JSON = "results/experiment3/verification_3ois.json"


def main() -> int:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("3ois", SOURCE_PDB)
    chain = structure[0][CHAIN_ID]

    # 1. HETATM inventory (chain A only; the other 3 chains are crystallographic copies).
    hetatm = {}
    for res in chain:
        if res.id[0] not in (" ", "W") and res.get_resname() != "HOH":
            hetatm.setdefault(res.get_resname(), []).append(res.id)
    logger.info(f"HETATM in chain A: { {k: len(v) for k, v in hetatm.items()} }")

    # 2. Catalytic Cys via sequence motif scan.
    residues = [r for r in chain if r.id[0] == " " and r.has_id("CA")]
    seq = "".join(three_to_one.get(r.get_resname(), "X") for r in residues)
    motif_idx = seq.find(CATALYTIC_MOTIF)
    if motif_idx < 0:
        logger.error(f"Catalytic motif {CATALYTIC_MOTIF} not found -- verification failed.")
        return 1
    cys_res = residues[motif_idx + 2]  # 'C' is the 3rd character of "GSCTANA"
    assert cys_res.get_resname() == "CYS"
    cys_id = cys_res.id[1]
    logger.info(f"Catalytic Cys candidate: Cys{cys_id} (motif at seq index {motif_idx})")

    res_by_id = {r.id[1]: r for r in residues}
    sg = cys_res["SG"].get_coord()

    # His within ion-pair/H-bond range of the Cys thiol.
    his_candidates = [r for r in residues if r.get_resname() == "HIS"]
    his_dists = []
    for h in his_candidates:
        for atomname in ("NE2", "ND1"):
            if h.has_id(atomname):
                his_dists.append((np.linalg.norm(sg - h[atomname].get_coord()), h, atomname))
    his_dists.sort(key=lambda t: t[0])
    best_his_dist, his_res, his_atom = his_dists[0]
    his_id = his_res.id[1]
    logger.info(f"Nearest His to Cys{cys_id} SG: His{his_id} {his_atom} at "
                f"{best_his_dist:.2f} A")

    # Asn within H-bond range of that His.
    asn_candidates = [r for r in residues if r.get_resname() == "ASN"]
    his_nd1 = his_res["ND1"].get_coord() if his_res.has_id("ND1") else his_res[his_atom].get_coord()
    asn_dists = []
    for a in asn_candidates:
        if a.has_id("OD1"):
            asn_dists.append((np.linalg.norm(his_nd1 - a["OD1"].get_coord()), a))
    asn_dists.sort(key=lambda t: t[0])
    best_asn_dist, asn_res = asn_dists[0]
    asn_id = asn_res.id[1]
    logger.info(f"Nearest Asn to His{his_id} ND1: Asn{asn_id} OD1 at {best_asn_dist:.2f} A")

    triad = {"cys": cys_id, "his": his_id, "asn": asn_id}
    active_site_centroid = np.mean([
        sg, his_res["NE2"].get_coord(), asn_res["OD1"].get_coord()
    ], axis=0)

    # 3. UDP distance to the active site.
    udp_res = next((r for r in chain if r.get_resname() == "UDP"), None)
    udp_result = None
    if udp_res is not None:
        udp_coords = np.array([a.get_coord() for a in udp_res if a.element != "H"])
        d_centroid = float(np.linalg.norm(udp_coords.mean(axis=0) - active_site_centroid))
        d_min = float(min(np.linalg.norm(c - sg) for c in udp_coords))
        udp_result = {"centroid_distance_A": round(d_centroid, 2),
                      "closest_atom_to_cys_sg_A": round(d_min, 2)}
        logger.info(f"UDP: centroid {d_centroid:.2f} A from active-site centroid, "
                    f"closest atom {d_min:.2f} A from Cys{cys_id} SG")
        if d_min < ACTIVE_SITE_GATE_ANGSTROM:
            logger.error(f"UDP is within {ACTIVE_SITE_GATE_ANGSTROM} A of the active site -- "
                          f"the site is NOT apo. Target premise fails. Stopping.")
            return 1
        logger.info(f"GATE PASSED: UDP is peripheral (> {ACTIVE_SITE_GATE_ANGSTROM} A). "
                    f"Active site is apo.")
    else:
        logger.warning("No UDP found in chain A.")

    # 4. N-terminal occlusion profile.
    occlusion = []
    for r in residues:
        rn = r.id[1]
        if rn >= cys_id:
            break
        dmin = min(np.linalg.norm(a.get_coord() - active_site_centroid)
                   for a in r if a.element != "H")
        occlusion.append({"resid": rn, "resname": r.get_resname(), "dist_to_active_site_A": round(float(dmin), 2)})

    occluding = [o for o in occlusion if o["dist_to_active_site_A"] <= 10.0]
    logger.info(f"Residues within 10 A of the active-site centroid (N-terminal to Cys{cys_id}): "
                f"{[o['resid'] for o in occluding]}")

    result = {
        "hetatm_inventory": {k: len(v) for k, v in hetatm.items()},
        "catalytic_triad": triad,
        "udp": udp_result,
        "n_terminal_occlusion_profile": occlusion,
        "occluding_residues_within_10A": [o["resid"] for o in occluding],
    }
    import os
    os.makedirs("results/experiment3", exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    logger.info(f"Written -> {OUTPUT_JSON}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

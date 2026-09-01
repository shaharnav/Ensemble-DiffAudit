"""
Phase 1b: structural screen over the PocketMiner apo/holo pair set.

No GPU, no docking -- PDB downloads and RMSD math, reusing Phase 2a's geometry code
(apo_holo_geometry.compute_apo_holo_geometry) so the screen and the final target's Phase 2
numbers come from the identical method.

Source: Meller et al., Nat Commun 2023 (PocketMiner), supplementary table
"validation_and_test_sets" from https://github.com/Mickdub/gvp/tree/pocket_pred
(data/pm-dataset/supplementary-tables.xlsx), reduced to `pocketminer_pairs.csv`: 43 apo/holo
pairs with a resolved holo ligand (6 apo-only "highly rigid protein" rows excluded -- there is
no holo pocket to measure).

Usage:
    ./venv/bin/python target_screen.py
"""
import csv
import logging
import os
import re
import sys
import urllib.request

import numpy as np
from Bio.PDB import PDBParser

from apo_holo_geometry import compute_apo_holo_geometry

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

PAIRS_CSV = "pocketminer_pairs.csv"
PDB_DIR = os.path.join("pdbs", "pocketminer_screen")
OUTPUT_CSV = "target_screen.csv"

# Waters, common cryoprotectants/buffer components, and monatomic ions -- never disqualifying
# on their own, and excluded from ligand-centroid definition.
COMMON_IGNORE = {
    "HOH", "WAT", "HTO", "GOL", "EDO", "PEG", "SO4", "PO4", "DMS", "ACT", "ACY", "TRS", "CL",
    "NA", "MG", "ZN", "CA", "K", "CS", "YB", "FMT", "MPD", "IMD", "BME", "1PE", "12P", "PGE",
    "EPE", "BOG", "MN", "CU", "NI", "CD", "CO", "IOD", "BR", "UNK", "SCN", "NO3", "MRD", "CIT",
    "TAM", "LDA", "MB3", "BTB",
    "MSE",  # selenomethionine -- a modified protein residue (SeMet phasing), not a pocket ligand
}

# Selection band (criterion 3): clear of trypsin's ~0.2 A noise floor, inside the ConforMix
# --twist-target-stop 2.0 A shell.
CA_RMSD_MIN = 1.0
CA_RMSD_MAX = 2.5


def _download(pdb_id: str) -> str:
    path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    if not os.path.exists(path) or os.path.getsize(path) < 500:
        os.makedirs(PDB_DIR, exist_ok=True)
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.pdb", path)
    return path


def parse_ligand_codes(field: str) -> set[str]:
    """'NDP,TOP' -> {'NDP','TOP'}; '2xCHD' -> {'CHD'}; 'AHK:403' -> {'AHK'};
    'GDP-Ca,ERY' -> {'GDP','CA','ERY'}; '2xCIT-2xFe' -> {'CIT','FE'}."""
    codes = set()
    for token in field.split(","):
        for part in token.split("-"):
            part = re.sub(r"^\d+x", "", part.strip(), flags=re.IGNORECASE)
            part = part.split(":")[0].strip().upper()
            if part:
                codes.add(part)
    return codes


def hetero_species(pdb_path: str, chain_id: str, exclude_codes: set[str]) -> set[str]:
    """Non-water, non-solvent, non-ion HETATM residue codes present in a chain, other than
    the annotated ligand -- a proxy for 'something else is bound here' (e.g. a cofactor)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    exclude = exclude_codes | COMMON_IGNORE
    species = set()
    try:
        chain = structure[0][chain_id]
    except KeyError:
        return species
    for res in chain:
        if res.id[0].startswith("H_"):
            code = res.id[0][2:]
            if code not in exclude:
                species.add(code)
    return species


def screen_one(row: dict) -> dict:
    apo_pdb_id, apo_chain = row["apo_pdb"], row["apo_chain"]
    holo_pdb_id, holo_chain = row["holo_pdb"], row["holo_chain"]
    ligand_field = row["ligand_field"]

    out = {
        "apo_pdb": apo_pdb_id, "apo_chain": apo_chain,
        "holo_pdb": holo_pdb_id, "holo_chain": holo_chain,
        "ligand_field": ligand_field, "motion_type": row["motion_type"],
        "status": "ok",
    }

    try:
        apo_path = _download(apo_pdb_id)
        holo_path = _download(holo_pdb_id)
        ligand_codes_all = parse_ligand_codes(ligand_field)
        ligand_codes_for_centroid = ligand_codes_all - COMMON_IGNORE
        if not ligand_codes_for_centroid:
            out["status"] = "skip: no non-ion ligand code"
            return out

        geometry, pocket_pairs, rot, tran = compute_apo_holo_geometry(
            apo_pdb=apo_path, holo_pdb=holo_path,
            apo_chain=apo_chain, holo_chain=holo_chain,
            ligand_codes=ligand_codes_for_centroid,
        )
        out.update({
            "n_residues_matched": geometry["n_residues_matched"],
            "n_pocket_residues": geometry["n_pocket_residues"],
            "global_ca_rmsd": round(geometry["global_ca_rmsd"], 3),
            "pocket_ca_rmsd": round(geometry["apo_holo_pocket_ca_rmsd"], 3),
            "pocket_allatom_rmsd": round(geometry["apo_holo_pocket_allatom_rmsd"], 3),
            "ca_to_allatom_ratio": round(
                geometry["apo_holo_pocket_ca_rmsd"] / geometry["apo_holo_pocket_allatom_rmsd"], 3
            ) if geometry["apo_holo_pocket_allatom_rmsd"] > 0 else None,
            "max_residue_displacement": round(geometry["max_residue_displacement"], 3),
            "max_residue_displacement_holo_resid": geometry["max_residue_displacement_holo_resid"],
        })

        apo_extra = hetero_species(apo_path, apo_chain, ligand_codes_all)
        holo_extra = hetero_species(holo_path, holo_chain, ligand_codes_all)
        # A multi-component annotated ligand field (e.g. "NDP,TOP") is ambiguous about
        # which code is the drug-like ligand under study and which is a second bound
        # species (often a cofactor) -- so it's treated as a cofactor-difference risk on
        # its own, not just when an *unlisted* extra species shows up. This is what
        # catches DHFR (2W9T/2W9S, "NDP,TOP"): NDP is annotated as part of "the ligand"
        # but is a cofactor absent from apo, which a naive "extra species" diff misses
        # because NDP is excluded from that comparison as a documented ligand component.
        cofactor_diff = bool(apo_extra) or bool(holo_extra) or len(ligand_codes_for_centroid) > 1
        out["cofactor_difference"] = cofactor_diff
        out["apo_extra_hetero"] = ",".join(sorted(apo_extra)) or None
        out["holo_extra_hetero"] = ",".join(sorted(holo_extra)) or None

        in_band = CA_RMSD_MIN <= geometry["apo_holo_pocket_ca_rmsd"] <= CA_RMSD_MAX
        backbone_driven = (
            out["ca_to_allatom_ratio"] is not None and out["ca_to_allatom_ratio"] >= 0.5
        )
        out["qualifies"] = bool(in_band and backbone_driven and not cofactor_diff)

    except Exception as e:
        out["status"] = f"error: {e}"

    return out


def main() -> int:
    with open(PAIRS_CSV) as f:
        rows = list(csv.DictReader(f))
    logger.info(f"Screening {len(rows)} PocketMiner apo/holo pairs "
                f"(band: {CA_RMSD_MIN}-{CA_RMSD_MAX} A pocket CA RMSD, "
                f"no cofactor difference, backbone-driven).")

    results = []
    for i, row in enumerate(rows, 1):
        logger.info(f"[{i}/{len(rows)}] {row['apo_pdb']} -> {row['holo_pdb']} "
                     f"({row['ligand_field']}, {row['motion_type']})")
        out = screen_one(row)
        results.append(out)
        if out["status"] != "ok":
            logger.warning(f"  {out['status']}")
        else:
            logger.info(f"  pocket_ca_rmsd={out['pocket_ca_rmsd']} A, "
                         f"ratio={out['ca_to_allatom_ratio']}, "
                         f"cofactor_diff={out['cofactor_difference']}, "
                         f"qualifies={out['qualifies']}")

    fieldnames = [
        "apo_pdb", "apo_chain", "holo_pdb", "holo_chain", "ligand_field", "motion_type",
        "status", "n_residues_matched", "n_pocket_residues", "global_ca_rmsd",
        "pocket_ca_rmsd", "pocket_allatom_rmsd", "ca_to_allatom_ratio",
        "max_residue_displacement", "max_residue_displacement_holo_resid",
        "cofactor_difference", "apo_extra_hetero", "holo_extra_hetero", "qualifies",
    ]
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r.get(k) for k in fieldnames})
    logger.info(f"Written -> {OUTPUT_CSV}")

    qualifying = [r for r in results if r.get("qualifies")]
    qualifying.sort(key=lambda r: r["ca_to_allatom_ratio"], reverse=True)
    logger.info(f"{len(qualifying)}/{len(results)} pairs qualify on all criteria.")
    for r in qualifying:
        logger.info(f"  {r['apo_pdb']}/{r['apo_chain']} -> {r['holo_pdb']}/{r['holo_chain']} "
                     f"({r['ligand_field']}): pocket_ca_rmsd={r['pocket_ca_rmsd']}, "
                     f"ratio={r['ca_to_allatom_ratio']}, motion={r['motion_type']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

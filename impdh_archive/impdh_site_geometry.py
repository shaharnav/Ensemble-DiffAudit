"""
Experiment 3, Phase 2a: locate both candidate pocket sites on CLas IMPDH and run the DSSP
gate on each before any box-validation work or GPU time.

Source: PDB 4QM1 (Bacillus anthracis IMPDH, CBS-deletion construct, 2.80 A) -- same family,
same construct type as CLas 6KCF, co-crystallized with BOTH IMP (substrate site) and
inhibitor 39H/D67 (NAD-adenosine subsite; RCSB: "this same pocket is utilized by the
bacterial IMPDH-specific NAD+-binding mode"). Measured directly in 4QM1's own liganded
frame, not transferred into CLas's frame first -- CLas 6KCF is missing density for large
stretches near both sites (the CBS-deletion junction and, critically, the flap loop
immediately preceding the catalytic Cys, residues 297-302 in CLas's own numbering), so an
apo, partially-disordered structure cannot reliably answer "is this site loop-dominated" on
its own account; measuring on the liganded ortholog, where the site is actually ordered,
is the more defensible test of the site's intrinsic fold composition.

Pocket-lining definition: any residue with a heavy atom within 8 A of ANY ligand heavy atom
(not the ligand centroid) -- a centroid-only cutoff badly undercounts residues around an
elongated, multi-ring ligand like 39H/D67, which was checked and confirmed here (centroid-only
gave n=3 for the NAD site; any-atom gives n=20-55 depending on chain scope, a fivefold-plus
difference that changes the gate's conclusion).

Interface awareness: 39H/D67 was found here to sit only 3.22 A from chain B (checked
explicitly, not assumed) -- consistent with the literature description of bacterial IMPDH's
NAD-binding mode as structurally distinct from the canonical Rossmann fold. Chain-A-only
pocket-lining residues for the NAD site are therefore incomplete; both A and B are scored
and combined for that site.

Usage:
    ./venv/bin/python impdh_site_geometry.py
"""
import json
import logging
import os
import sys

import numpy as np
from Bio.PDB import PDBParser
import biotite.structure as struc
import biotite.structure.io.pdb as biotite_pdb

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

SOURCE_PDB = os.path.join("pdbs", "exp3_impdh", "4QM1.pdb")
POCKET_CUTOFF_A = 8.0
DSSP_FLOOR = 0.60
INTERFACE_CHECK_CUTOFF_A = 5.0
OUTPUT_JSON = os.path.join("results", "experiment3", "phase2a_dssp_gate.json")


def get_ligand_atoms(structure, chain_id, ligand_code):
    for res in structure[0][chain_id]:
        if res.get_resname() == ligand_code:
            return np.array([a.get_coord() for a in res if a.element != "H"])
    raise ValueError(f"{ligand_code} not found in chain {chain_id}")


def pocket_lining_residues_any_atom(structure, chain_id, ligand_atoms, cutoff=POCKET_CUTOFF_A):
    if chain_id not in structure[0]:
        return []
    res_by_id = {r.id[1]: r for r in structure[0][chain_id] if r.id[0] == " "}
    lining = []
    for rid, res in res_by_id.items():
        res_coords = np.array([a.get_coord() for a in res if a.element != "H"])
        if res_coords.size == 0:
            continue
        d = np.linalg.norm(res_coords[:, None, :] - ligand_atoms[None, :, :], axis=2).min()
        if d <= cutoff:
            lining.append(rid)
    return sorted(lining)


def min_ligand_distance_to_chain(structure, chain_id, ligand_atoms):
    if chain_id not in structure[0]:
        return None
    other_atoms = np.array([a.get_coord() for a in structure[0][chain_id].get_atoms()
                             if a.element != "H"])
    return float(np.linalg.norm(other_atoms[None, :, :] - ligand_atoms[:, None, :], axis=2).min())


def sse_by_residue(pdb_path, chain_id):
    pdb_file = biotite_pdb.PDBFile.read(pdb_path)
    atoms = biotite_pdb.get_structure(pdb_file, model=1)
    atoms = atoms[(atoms.chain_id == chain_id) & struc.filter_amino_acids(atoms)]
    sse = struc.annotate_sse(atoms)
    res_ids_order = np.unique(atoms.res_id)
    return dict(zip(res_ids_order, sse))


def structured_fraction(sse_maps_by_chain, lining_by_chain):
    counts = {"a": 0, "b": 0, "c": 0}
    n = 0
    for chain_id, ids in lining_by_chain.items():
        sse_map = sse_maps_by_chain.get(chain_id, {})
        for rid in ids:
            code = sse_map.get(rid)
            if code is None:
                continue
            counts[code] = counts.get(code, 0) + 1
            n += 1
    structured = counts.get("a", 0) + counts.get("b", 0)
    frac = structured / n if n else 0.0
    return frac, n, counts


def main() -> int:
    if not os.path.isfile(SOURCE_PDB):
        logger.error(f"{SOURCE_PDB} not found.")
        return 1
    os.makedirs(os.path.dirname(OUTPUT_JSON), exist_ok=True)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("4qm1", SOURCE_PDB)
    all_chains = sorted(c.id for c in structure[0])
    logger.info(f"4QM1 chains: {all_chains}")

    site_ligands = {"IMP_substrate_site": ("A", "IMP"), "NAD_cofactor_site": ("A", "39H")}
    results = {}

    for site_name, (chain_id, ligand_code) in site_ligands.items():
        ligand_atoms = get_ligand_atoms(structure, chain_id, ligand_code)

        # Interface check: is this ligand meaningfully closer to a neighboring chain than
        # to bulk solvent? Determines whether pocket-lining residues must be pooled across
        # chains for this site.
        interface_distances = {}
        for other in all_chains:
            if other == chain_id:
                continue
            interface_distances[other] = min_ligand_distance_to_chain(structure, other, ligand_atoms)
        interface_chains = [c for c, d in interface_distances.items()
                             if d is not None and d <= INTERFACE_CHECK_CUTOFF_A]
        chains_to_score = [chain_id] + interface_chains
        logger.info(f"{site_name}: ligand-to-other-chain min distances = {interface_distances} "
                    f"-> scoring chains {chains_to_score}"
                    f"{' (interface site)' if interface_chains else ''}")

        lining_by_chain = {
            c: pocket_lining_residues_any_atom(structure, c, ligand_atoms) for c in chains_to_score
        }
        sse_maps = {c: sse_by_residue(SOURCE_PDB, c) for c in chains_to_score}
        frac, n_scored, counts = structured_fraction(sse_maps, lining_by_chain)
        passed = frac >= DSSP_FLOOR
        total_lining = sum(len(v) for v in lining_by_chain.values())
        logger.info(f"{site_name}: {total_lining} pocket-lining residues across "
                    f"{chains_to_score} (any heavy atom within {POCKET_CUTOFF_A} A of any "
                    f"{ligand_code} heavy atom): {lining_by_chain}")
        logger.info(f"{site_name}: DSSP structured fraction = {frac:.1%} ({n_scored} scored, "
                    f"helix={counts.get('a',0)}, sheet={counts.get('b',0)}, coil={counts.get('c',0)}) "
                    f"-> {'PASS' if passed else 'FAIL'} (floor {DSSP_FLOOR:.0%})")

        results[site_name] = {
            "ligand_code": ligand_code,
            "interface_distances_A": interface_distances,
            "chains_scored": chains_to_score,
            "pocket_lining_residues_by_chain": lining_by_chain,
            "n_lining_residues_total": total_lining,
            "n_scored_for_sse": n_scored,
            "structured_fraction": round(frac, 4),
            "helix_count": counts.get("a", 0),
            "sheet_count": counts.get("b", 0),
            "coil_count": counts.get("c", 0),
            "gate_passed": passed,
        }

    either_pass = any(r["gate_passed"] for r in results.values())
    both_pass = all(r["gate_passed"] for r in results.values())
    if not either_pass:
        logger.error("Neither candidate site passed the DSSP floor (>=60% helix/sheet) -- "
                      "this target has the same defect as Xylellain (a coil-dominated site "
                      "ConforMix cannot be expected to sample plausibly). Per the "
                      "pre-registered plan: stop and report before any GPU time.")
    elif both_pass:
        logger.info("Both sites passed -- preferring IMP_substrate_site per the plan.")
    else:
        chosen = next(k for k, v in results.items() if v["gate_passed"])
        logger.info(f"Only {chosen} passed -- selecting it.")

    output = {
        "source_structure": "4QM1",
        "method": (
            "Measured on the liganded ortholog (4QM1), not on apo CLas 6KCF -- 6KCF is "
            "missing density at both candidate sites (CBS-deletion junction, flap loop "
            "residues 297-302 in CLas numbering), so its own backbone cannot support a "
            "fold-composition measurement there. Pocket-lining = any residue heavy atom "
            "within 8 A of any ligand heavy atom (not ligand centroid -- centroid-only was "
            "tried first and undercounted by 5x for the elongated NAD-site ligand). "
            "NAD site is scored across both chain A and its interface partner chain B "
            "(39H found 3.22 A from chain B, confirmed as a genuine subunit-interface site)."
        ),
        "pocket_cutoff_A": POCKET_CUTOFF_A,
        "dssp_floor": DSSP_FLOOR,
        "sites": results,
        "either_site_passed": either_pass,
        "both_sites_passed": both_pass,
    }
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)
    logger.info(f"Written -> {OUTPUT_JSON}")

    return 0 if either_pass else 1


if __name__ == "__main__":
    sys.exit(main())

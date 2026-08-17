"""
Experiment 3, Phase 3: per-conformer geometry for the 6 ConforMix conformers
generated from R1 (the zymogen).

R2 is R1 with residues 23-66 deleted, no coordinates moved, so R1 and R2 share the
same absolute frame and numbering. The conformers do not: despite the notebook's
description of ConforMix output as "locked to the CIF reference coordinates,"
Boltz's actual output (a) is in its own arbitrary frame, not R1's, and (b)
renumbers residues from 1 regardless of R1's native numbering (verified:
conformer residue 1's sequence matches R1 residue 23's) -- both corrected here
(`load_residues(..., renumber_offset=...)`, then a CA superposition per
conformer) before any distance is measured.

The superposition itself anchors on the mature domain only (residues >= 67), not
the whole chain: a whole-structure fit average-fits two subdomains (the
propeptide and the mature domain) that can move independently, which produces a
poor and misleading global RMSD if the propeptide really did move a lot. Anchor
residues are stationary by definition of the fit, so this only gets used to ask
"how far did the (non-anchor) propeptide move relative to the (fixed) mature
domain" -- exactly the question that matters here.

**Note on pocket_ca_rmsd_to_R2, as literally specified in the plan:** R2 is an
*unrelaxed* deletion of R1 (Phase 1c) -- it does not carry a moved/relaxed
geometry for the residues it retains, only the absence of the propeptide.
For every pocket-lining residue outside the propeptide range (23-66), R2's
coordinates are byte-identical to R1's, so pocket_ca_rmsd_to_R2 for those
residues is mathematically identical to pocket_ca_rmsd_to_R1 -- not a second,
independent signal. It's still computed and reported (per the plan), with
this explained rather than silently glossed over. The genuinely informative
questions -- did the propeptide move, and did the cavity open -- are answered
by `propeptide_ca_rmsd_to_R1` and `cavity_volume_A3` instead.

Usage:
    ./venv/bin/python exp3_conformer_geometry.py
"""
import csv
import glob
import json
import logging
import os
import sys

import numpy as np
import yaml
from Bio.PDB import PDBParser, Superimposer

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

R1_PATH = os.path.join("results", "experiment3", "receptors", "R1_zymogen.pdb")
R2_PATH = os.path.join("results", "experiment3", "receptors", "R2_mature_proxy.pdb")
CONFORMER_DIR = os.path.join("results", "experiment3", "apo_ensemble")
POCKET_CONFIG_YAML = "pocket_config.yaml"
CHAIN_ID = "A"
PROPEPTIDE_RANGE = (23, 66)
OCCLUDING_LOOP_RANGE = (26, 36)
POCKET_CUTOFF = 8.0
GRID_SPACING = 0.5
POCKET_RADIUS = 10.0
CLASH_MIN = 1.4
OUTPUT_CSV = "experiment3_conformer_geometry.csv"

# Phase 3 plausibility gate -- same status as the Phase 2 volume gate: a conformer
# failing any of these is excluded before it can silently corrupt C_mid selection or
# Phase 4 conditioning. Thresholds per plan: propeptide CA RMSD > 15 A, propeptide
# COM-to-domain-surface distance increase > 10 A vs. R1, or cavity volume > R2's.
GATE_PROPEPTIDE_CA_RMSD_MAX = 15.0
GATE_COM_SURFACE_DELTA_MAX = 10.0
# Cavity-volume gate threshold is R2's own measured volume, computed at runtime
# (there's no free-standing constant for it -- it depends on which structure ran).


class _FakeAtom:
    def __init__(self, coord):
        self.coord = np.asarray(coord, dtype=float)

    def get_coord(self):
        return self.coord


def load_residues(pdb_path, renumber_offset=0):
    """renumber_offset: Boltz/ConforMix output renumbers from 1 regardless of the input
    structure's own numbering (confirmed: conformer residue 1 has R1 residue 23's
    sequence -- identical 268-residue sequence, just reindexed). Passing
    renumber_offset=22 (R1's first residue id, 23, minus 1) maps conformer ids back onto
    R1/R2's native numbering so residue-id-keyed comparisons are actually comparing the
    same residue."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(pdb_path), pdb_path)
    return {
        r.id[1] + renumber_offset: r
        for r in structure[0][CHAIN_ID] if r.id[0] == " " and r.has_id("CA")
    }


def find_pocket_lining_residue_ids(res_by_id, centroid, cutoff=POCKET_CUTOFF):
    lining = []
    for rid, res in res_by_id.items():
        for atom in res:
            if atom.element == "H":
                continue
            if np.linalg.norm(atom.get_coord() - centroid) <= cutoff:
                lining.append(rid)
                break
    return sorted(lining)


def superimpose(mobile_res_by_id, fixed_res_by_id, anchor_ids=None):
    """CA-only superposition. anchor_ids restricts which residues define the fit --
    default (None) uses every residue id present in both, but a whole-structure fit is
    only meaningful if the whole structure actually moved as one rigid body. Here it
    doesn't: the propeptide can swing independently of the mature domain, so a
    same-frame comparison should anchor on the (expected-stationary) mature domain,
    not average-fit two subdomains that may have moved relative to each other."""
    common_ids = sorted(set(mobile_res_by_id) & set(fixed_res_by_id))
    if anchor_ids is not None:
        common_ids = [i for i in common_ids if i in anchor_ids]
    mobile_ca = np.array([mobile_res_by_id[i]["CA"].get_coord() for i in common_ids])
    fixed_ca = np.array([fixed_res_by_id[i]["CA"].get_coord() for i in common_ids])
    sup = Superimposer()
    sup.set_atoms([_FakeAtom(c) for c in fixed_ca], [_FakeAtom(c) for c in mobile_ca])
    rot, tran = sup.rotran
    return rot, tran, sup.rms, common_ids


def ca_rmsd(mobile_res_by_id, fixed_res_by_id, rot, tran, residue_ids):
    ids = [i for i in residue_ids if i in mobile_res_by_id and i in fixed_res_by_id]
    if not ids:
        return None, 0
    mobile = np.array([mobile_res_by_id[i]["CA"].get_coord() for i in ids]) @ rot.T + tran
    fixed = np.array([fixed_res_by_id[i]["CA"].get_coord() for i in ids])
    d = np.linalg.norm(mobile - fixed, axis=1)
    return float(np.sqrt(np.mean(d ** 2))), len(ids)


def allatom_rmsd(mobile_res_by_id, fixed_res_by_id, rot, tran, residue_ids):
    deltas = []
    for i in residue_ids:
        if i not in mobile_res_by_id or i not in fixed_res_by_id:
            continue
        mres, fres = mobile_res_by_id[i], fixed_res_by_id[i]
        fatoms = {a.get_name(): a for a in fres if a.element != "H"}
        for matom in mres:
            if matom.element == "H":
                continue
            fatom = fatoms.get(matom.get_name())
            if fatom is None:
                continue
            mcoord = matom.get_coord() @ rot.T + tran
            deltas.append(np.linalg.norm(mcoord - fatom.get_coord()))
    return float(np.sqrt(np.mean(np.square(deltas)))) if deltas else None


def cavity_volume_from_atoms(atoms, centroid):
    nearby = atoms[np.linalg.norm(atoms - centroid, axis=1) <= POCKET_RADIUS + CLASH_MIN]
    if len(nearby) == 0:
        # No protein atoms anywhere near the fixed center at all -- physically
        # implausible for a folded mature domain that stayed put (mature_anchor_rmsd
        # is checked separately), almost certainly a garbage/outlier structure rather
        # than a genuinely wide-open pocket. Flag it rather than silently reporting
        # the trivial full-sphere volume as if it were a real measurement.
        logger.warning("No atoms within range of the pocket center -- suspicious, "
                        "not a real cavity-volume measurement. Returning None.")
        return None
    axis = np.arange(-POCKET_RADIUS, POCKET_RADIUS + 1e-6, GRID_SPACING)
    n_open = 0
    for x in axis:
        for y in axis:
            pts = np.stack([
                np.full_like(axis, centroid[0] + x),
                np.full_like(axis, centroid[1] + y),
                centroid[2] + axis,
            ], axis=1)
            in_sphere = np.linalg.norm(pts - centroid, axis=1) <= POCKET_RADIUS
            if not np.any(in_sphere):
                continue
            pts = pts[in_sphere]
            d = np.linalg.norm(nearby[None, :, :] - pts[:, None, :], axis=2).min(axis=1)
            n_open += int(np.sum(d >= CLASH_MIN))
    return round(n_open * GRID_SPACING ** 3, 1)


def cavity_volume(pdb_path, centroid):
    """R1/R2 are already in pocket_config.yaml's frame -- no transform needed."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    atoms = np.array([a.get_coord() for a in structure.get_atoms() if a.element != "H"])
    return cavity_volume_from_atoms(atoms, centroid)


def propeptide_com_to_surface_distance(propeptide_res_by_id, mature_atoms):
    """Center of mass of the propeptide's CA atoms, then its nearest-atom distance to
    the mature domain -- a proxy for whether the propeptide is still resting against
    the folded domain (near 0) or has drifted away into free space."""
    com = np.mean([r["CA"].get_coord() for r in propeptide_res_by_id.values()], axis=0)
    return float(np.min(np.linalg.norm(mature_atoms - com, axis=1)))


def plausibility_gate(propeptide_ca_rmsd, com_surface_delta, cavity_vol, cavity_max):
    """Returns (passed: bool, reasons: list[str])."""
    reasons = []
    if propeptide_ca_rmsd is not None and propeptide_ca_rmsd > GATE_PROPEPTIDE_CA_RMSD_MAX:
        reasons.append(f"propeptide CA RMSD {propeptide_ca_rmsd:.1f} A > "
                        f"{GATE_PROPEPTIDE_CA_RMSD_MAX} A")
    if com_surface_delta is not None and com_surface_delta > GATE_COM_SURFACE_DELTA_MAX:
        reasons.append(f"propeptide COM-to-surface distance increased by "
                        f"{com_surface_delta:.1f} A > {GATE_COM_SURFACE_DELTA_MAX} A vs. R1")
    if cavity_vol is not None and cavity_max is not None and cavity_vol > cavity_max:
        reasons.append(f"cavity volume {cavity_vol:.1f} A^3 exceeds R2's {cavity_max:.1f} A^3")
    return len(reasons) == 0, reasons


def main() -> int:
    conformer_paths = sorted(glob.glob(os.path.join(CONFORMER_DIR, "conformix_var_*.pdb")))
    if not conformer_paths:
        logger.error(f"No conformers found in {CONFORMER_DIR}.")
        return 1
    for path in (R1_PATH, R2_PATH, POCKET_CONFIG_YAML):
        if not os.path.exists(path):
            logger.error(f"{path} not found.")
            return 1

    with open(POCKET_CONFIG_YAML) as f:
        pocket_config = yaml.safe_load(f)
    centroid = np.array(pocket_config["center"])

    r1_res = load_residues(R1_PATH)
    r2_res = load_residues(R2_PATH)
    conformer_offset = min(r1_res) - 1  # Boltz/ConforMix renumbers output from 1
    pocket_ids = find_pocket_lining_residue_ids(r2_res, centroid)
    propeptide_ids = [i for i in r1_res if PROPEPTIDE_RANGE[0] <= i <= PROPEPTIDE_RANGE[1]]
    occluding_ids = [i for i in r1_res if OCCLUDING_LOOP_RANGE[0] <= i <= OCCLUDING_LOOP_RANGE[1]]
    mature_ids = [i for i in r1_res if i >= PROPEPTIDE_RANGE[1] + 1]
    logger.info(f"{len(pocket_ids)} pocket-lining residues (8 A of pocket_config.yaml center): {pocket_ids}")
    logger.info(f"{len(propeptide_ids)} propeptide residues, {len(occluding_ids)} in the occluding loop proper")

    r1_volume = cavity_volume(R1_PATH, centroid)
    r2_volume = cavity_volume(R2_PATH, centroid)
    logger.info(f"R1 cavity: {r1_volume} A^3, R2 cavity: {r2_volume} A^3 (from receptors.yaml/pocket_volume.json)")

    # Baseline for the COM-to-surface gate: R1's own propeptide resting distance
    # (expected near 0 -- it's occluding the site, touching the domain) and the fixed
    # mature-domain atom set every conformer's propeptide COM is measured against.
    r1_propeptide_res = {i: r1_res[i] for i in propeptide_ids}
    r1_parser = PDBParser(QUIET=True)
    r1_structure = r1_parser.get_structure("r1", R1_PATH)
    r1_mature_atoms = np.array([
        a.get_coord() for a in r1_structure[0][CHAIN_ID].get_atoms()
        if a.get_parent().id[1] >= PROPEPTIDE_RANGE[1] + 1 and a.element != "H"
    ])
    r1_com_surface_dist = propeptide_com_to_surface_distance(r1_propeptide_res, r1_mature_atoms)
    logger.info(f"R1 baseline: propeptide COM-to-mature-domain-surface distance = "
                f"{r1_com_surface_dist:.2f} A")

    rows = []
    for conf_path in conformer_paths:
        name = os.path.basename(conf_path)
        conf_res = load_residues(conf_path, renumber_offset=conformer_offset)

        rot, tran, global_rmsd, common_ids = superimpose(conf_res, r1_res, anchor_ids=mature_ids)

        pocket_ca_r1, n_pocket = ca_rmsd(conf_res, r1_res, rot, tran, pocket_ids)
        pocket_allatom_r1 = allatom_rmsd(conf_res, r1_res, rot, tran, pocket_ids)
        pocket_ca_r2, n_pocket_r2 = ca_rmsd(conf_res, r2_res, rot, tran, pocket_ids)
        propeptide_ca_r1, n_prop = ca_rmsd(conf_res, r1_res, rot, tran, propeptide_ids)
        occluding_ca_r1, n_occ = ca_rmsd(conf_res, r1_res, rot, tran, occluding_ids)

        # Conformers come out of Boltz in their own frame, not R1's absolute frame --
        # transform every atom by the same CA-derived (rot, tran) before measuring
        # cavity volume against pocket_config.yaml's fixed center.
        conf_parser = PDBParser(QUIET=True)
        conf_structure = conf_parser.get_structure("c", conf_path)
        conf_atoms = np.array([a.get_coord() for a in conf_structure.get_atoms() if a.element != "H"])
        conf_atoms_transformed = conf_atoms @ rot.T + tran
        vol = cavity_volume_from_atoms(conf_atoms_transformed, centroid)

        # COM-to-surface gate input: this conformer's propeptide, transformed into
        # R1's frame, against R1's own (fixed) mature-domain atoms.
        conf_propeptide_res = {i: conf_res[i] for i in propeptide_ids if i in conf_res}
        conf_com = np.mean(
            [r["CA"].get_coord() @ rot.T + tran for r in conf_propeptide_res.values()], axis=0
        )
        conf_com_surface_dist = float(np.min(np.linalg.norm(r1_mature_atoms - conf_com, axis=1)))
        com_surface_delta = conf_com_surface_dist - r1_com_surface_dist

        gate_passed, gate_reasons = plausibility_gate(
            propeptide_ca_r1, com_surface_delta, vol, r2_volume
        )

        row = {
            "conformer": name,
            "mature_domain_anchor_rmsd_to_R1": round(global_rmsd, 3),
            "n_anchor_residues": len(common_ids),
            "pocket_ca_rmsd_to_R1": round(pocket_ca_r1, 3) if pocket_ca_r1 is not None else None,
            "pocket_allatom_rmsd_to_R1": round(pocket_allatom_r1, 3) if pocket_allatom_r1 is not None else None,
            "pocket_ca_rmsd_to_R2": round(pocket_ca_r2, 3) if pocket_ca_r2 is not None else None,
            "n_pocket_residues": n_pocket,
            "propeptide_ca_rmsd_to_R1": round(propeptide_ca_r1, 3) if propeptide_ca_r1 is not None else None,
            "occluding_loop_ca_rmsd_to_R1": round(occluding_ca_r1, 3) if occluding_ca_r1 is not None else None,
            "cavity_volume_A3": vol,
            "cavity_opening_vs_R1_A3": round(vol - r1_volume, 1) if vol is not None else None,
            "cavity_pct_of_R1_to_R2_gap": (
                round(100 * (vol - r1_volume) / (r2_volume - r1_volume), 1)
                if vol is not None and r2_volume != r1_volume else None
            ),
            "com_surface_delta_A": round(com_surface_delta, 2),
            "plausibility_gate_passed": gate_passed,
            "plausibility_gate_reasons": "; ".join(gate_reasons) if gate_reasons else None,
        }
        rows.append(row)
        if not gate_passed:
            logger.warning(f"{name}: PLAUSIBILITY GATE FAILED -- {'; '.join(gate_reasons)}")
        logger.info(
            f"{name}: mature_anchor_rmsd={row['mature_domain_anchor_rmsd_to_R1']}, "
            f"pocket_ca_to_R1={row['pocket_ca_rmsd_to_R1']}, "
            f"propeptide_ca={row['propeptide_ca_rmsd_to_R1']}, "
            f"occluding_loop_ca={row['occluding_loop_ca_rmsd_to_R1']}, "
            f"cavity={row['cavity_volume_A3']} A^3 "
            f"({row['cavity_pct_of_R1_to_R2_gap']}% of R1->R2 gap)"
        )

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    logger.info(f"Written -> {OUTPUT_CSV}")

    # C_mid: conformer with the median cavity volume among those that pass the
    # plausibility gate -- same status as the Phase 2 volume gate. A conformer that
    # fails is excluded, not just flagged, so it can't silently become C_mid or feed
    # Phase 4 conditioning.
    n_failed = sum(1 for r in rows if not r["plausibility_gate_passed"])
    if n_failed:
        logger.warning(f"{n_failed}/{len(rows)} conformer(s) FAILED the plausibility gate "
                        f"and are excluded from C_mid selection: "
                        f"{[r['conformer'] for r in rows if not r['plausibility_gate_passed']]}")
    valid_rows = [
        r for r in rows
        if r["plausibility_gate_passed"] and r["cavity_volume_A3"] is not None
    ]
    if not valid_rows:
        logger.error("No conformer passed the plausibility gate -- cannot select C_mid. "
                      "Regenerate Phase 3 before proceeding to Phase 4.")
        return 1
    sorted_by_vol = sorted(valid_rows, key=lambda r: r["cavity_volume_A3"])
    mid_idx = len(sorted_by_vol) // 2
    c_mid = sorted_by_vol[mid_idx]
    logger.info(f"C_mid (median cavity volume among {len(valid_rows)} gate-passing "
                f"conformers): {c_mid['conformer']} ({c_mid['cavity_volume_A3']} A^3)")

    return 0


if __name__ == "__main__":
    sys.exit(main())

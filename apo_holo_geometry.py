"""
Structural measurement and directional test for the Experiment 2 pilot (Phase 2).

Runs before any docking, so the independent variable (how far apo DHFR is from holo, and
whether the ConforMix ensemble moves toward holo) is fixed before any docking outcome is
seen.

2a — ground truth: aligns apo (1RTC) to holo (1BR6) chain A on the matched sequence region
(residues paired by sequence-alignment index, not PDB residue number -- see
`_sequence_align_residues`), and reports pocket displacement.

2b — directional test: given a 6-conformer ConforMix ensemble generated from the apo
structure (see targets.yaml / README "Generation (Colab)"), computes `directional_gain`
against an isotropic-displacement null. Requires
`results/experiment2/apo_ensemble/conformix_var_*.pdb` to exist; if that directory is
empty, 2b is skipped and 2a's results are reported alone.

Usage:
    ./venv/bin/python apo_holo_geometry.py
"""
import csv
import glob
import logging
import os
import sys

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import protein_letters_3to1 as _three_to_one

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

APO_PDB = "pdbs/1RTC.pdb"
HOLO_PDB = "pdbs/1BR6.pdb"
CHAIN_ID = "A"
LIGAND_CODE = "PT1"
DEFAULT_CUTOFF = 8.0
CONFORMER_DIR = os.path.join("results", "experiment2", "apo_ensemble")
GEOMETRY_CSV = "apo_holo_geometry.csv"
CONFORMER_GEOMETRY_CSV = "experiment2_conformer_geometry.csv"
N_NULL_SAMPLES = 10000
RNG_SEED = 42


def _standard_ca_residues(structure, chain_id):
    return [
        res
        for res in structure[0][chain_id]
        if res.id[0] == " " and res.has_id("CA")
    ]


def _sequence_align_residues(apo_residues, holo_residues):
    """Pair residues by sequence-alignment index, not PDB residue number.

    Both chains are the same 161-residue construct with the same numbering convention
    (verified via RCSB: identical sequence, no mutations) but different modeled extents
    at the termini, so a global alignment on the one-letter sequence -- rather than a
    naive positional zip -- is the correct correspondence and catches any residue that
    doesn't actually match.
    """
    apo_seq = "".join(_three_to_one[r.get_resname()] for r in apo_residues)
    holo_seq = "".join(_three_to_one[r.get_resname()] for r in holo_residues)

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(apo_seq, holo_seq)[0]

    pairs = []
    for (a_start, a_end), (h_start, h_end) in zip(*alignment.aligned):
        for offset in range(a_end - a_start):
            apo_res = apo_residues[a_start + offset]
            holo_res = holo_residues[h_start + offset]
            if apo_res.get_resname() == holo_res.get_resname():
                pairs.append((apo_res, holo_res))
    return pairs


def get_ligand_centroid(structure, chain_id, ligand_codes):
    """ligand_codes: iterable of PDB HETATM residue codes (e.g. {"TOP"} or {"NDP", "TOP"})."""
    het_ids = {f"H_{code}" for code in ligand_codes}
    coords = []
    for res in structure[0][chain_id]:
        if res.id[0] in het_ids:
            coords.extend(atom.get_coord() for atom in res if atom.element != "H")
    if not coords:
        raise ValueError(f"None of {ligand_codes} found in chain {chain_id}")
    return np.mean(np.array(coords), axis=0)


def find_pocket_lining_residues(structure, ligand_centroid, cutoff, chain_id):
    lining = []
    for res in structure[0][chain_id]:
        if res.id[0] != " " or not res.has_id("CA"):
            continue
        for atom in res:
            if atom.element == "H":
                continue
            if np.linalg.norm(atom.get_coord() - ligand_centroid) <= cutoff:
                lining.append(res)
                break
    return lining


def compute_apo_holo_geometry(apo_pdb=None, holo_pdb=None, apo_chain=CHAIN_ID,
                               holo_chain=CHAIN_ID, ligand_codes=None, cutoff=DEFAULT_CUTOFF):
    apo_pdb = apo_pdb or APO_PDB
    holo_pdb = holo_pdb or HOLO_PDB
    ligand_codes = ligand_codes or {LIGAND_CODE}

    parser = PDBParser(QUIET=True)
    apo_structure = parser.get_structure("apo", apo_pdb)
    holo_structure = parser.get_structure("holo", holo_pdb)

    apo_residues = _standard_ca_residues(apo_structure, apo_chain)
    holo_residues = _standard_ca_residues(holo_structure, holo_chain)
    pairs = _sequence_align_residues(apo_residues, holo_residues)
    logger.info(f"{len(pairs)} residues matched by sequence-alignment index "
                f"({len(apo_residues)} apo, {len(holo_residues)} holo)")
    if len(pairs) < 20:
        raise ValueError(f"Only {len(pairs)} residues matched -- apo/holo chains likely "
                          f"aren't the same construct.")

    ligand_centroid = get_ligand_centroid(holo_structure, holo_chain, ligand_codes)
    logger.info(f"Holo ligand ({ligand_codes}) centroid: {ligand_centroid}")

    lining_holo = find_pocket_lining_residues(holo_structure, ligand_centroid, cutoff, holo_chain)
    lining_holo_ids = {id(r) for r in lining_holo}
    pocket_pairs = [(a, h) for a, h in pairs if id(h) in lining_holo_ids]
    logger.info(f"{len(pocket_pairs)} pocket-lining residues within {cutoff} A of ligand "
                f"centroid, present in both apo and holo")
    if len(pocket_pairs) < 3:
        raise ValueError(f"Only {len(pocket_pairs)} pocket-lining residues resolved in both "
                          f"structures -- too few for a stable RMSD.")

    # Global superposition: align apo onto holo using ALL matched CA atoms (not just
    # pocket), so pocket displacement is measured in a frame fixed by the whole protein,
    # not biased by the pocket itself.
    apo_ca_all = np.array([a["CA"].get_coord() for a, h in pairs])
    holo_ca_all = np.array([h["CA"].get_coord() for a, h in pairs])
    sup = Superimposer()
    sup.set_atoms(
        [_FakeAtom(c) for c in holo_ca_all],
        [_FakeAtom(c) for c in apo_ca_all],
    )
    global_ca_rmsd = sup.rms

    # Apply the same rotation/translation to every apo atom so pocket-only distances
    # are measured post-superposition.
    rot, tran = sup.rotran
    apo_pocket_ca = np.array([a["CA"].get_coord() for a, h in pocket_pairs]) @ rot + tran
    holo_pocket_ca = np.array([h["CA"].get_coord() for a, h in pocket_pairs])
    ca_deltas = np.linalg.norm(apo_pocket_ca - holo_pocket_ca, axis=1)
    apo_holo_pocket_ca_rmsd = float(np.sqrt(np.mean(ca_deltas ** 2)))

    allatom_deltas = []
    for a, h in pocket_pairs:
        h_atoms_by_name = {atom.get_name(): atom for atom in h if atom.element != "H"}
        for a_atom in a:
            if a_atom.element == "H":
                continue
            h_atom = h_atoms_by_name.get(a_atom.get_name())
            if h_atom is None:
                continue
            a_coord = a_atom.get_coord() @ rot + tran
            allatom_deltas.append(np.linalg.norm(a_coord - h_atom.get_coord()))
    apo_holo_pocket_allatom_rmsd = float(np.sqrt(np.mean(np.square(allatom_deltas))))

    max_idx = int(np.argmax(ca_deltas))
    max_residue_displacement = float(ca_deltas[max_idx])
    max_residue_id = pocket_pairs[max_idx][1].id

    result = {
        "apo_pdb": apo_pdb,
        "holo_pdb": holo_pdb,
        "n_residues_matched": len(pairs),
        "n_pocket_residues": len(pocket_pairs),
        "cutoff_angstrom": cutoff,
        "global_ca_rmsd": global_ca_rmsd,
        "apo_holo_pocket_ca_rmsd": apo_holo_pocket_ca_rmsd,
        "apo_holo_pocket_allatom_rmsd": apo_holo_pocket_allatom_rmsd,
        "max_residue_displacement": max_residue_displacement,
        "max_residue_displacement_holo_resid": max_residue_id[1],
    }
    return result, pocket_pairs, rot, tran


class _FakeAtom:
    """Minimal Bio.PDB Atom stand-in -- Superimposer only needs .get_coord()/.coord."""

    def __init__(self, coord):
        self.coord = np.asarray(coord, dtype=float)

    def get_coord(self):
        return self.coord


def _isotropic_null_best_of_n(displacement_magnitude, distance_to_target, n_conformers,
                               n_resamples=N_NULL_SAMPLES, seed=RNG_SEED):
    """Monte Carlo null: if a point moves a fixed distance `displacement_magnitude` in a
    uniformly random 3D direction, what's the expected distance to a target `distance_to_target`
    away, taking the best (closest) of `n_conformers` independent draws?

    This is the null `directional_gain` must beat -- random displacement of the observed
    magnitude, with no preference for the holo direction, still occasionally lands closer to
    holo by chance when you take the best of several draws.
    """
    rng = np.random.default_rng(seed)
    origin_to_target = np.array([distance_to_target, 0.0, 0.0])
    best_distances = np.empty(n_resamples)
    for i in range(n_resamples):
        directions = rng.normal(size=(n_conformers, 3))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        points = directions * displacement_magnitude
        dists = np.linalg.norm(points - origin_to_target, axis=1)
        best_distances[i] = dists.min()
    return best_distances


def compute_directional_gain(pocket_pairs, rot, tran, apo_holo_pocket_ca_rmsd, cutoff=DEFAULT_CUTOFF):
    conformer_paths = sorted(glob.glob(os.path.join(CONFORMER_DIR, "conformix_var_*.pdb")))
    if not conformer_paths:
        logger.warning(
            f"No conformers in {CONFORMER_DIR} -- Phase 2b (directional test) needs the "
            f"6-conformer ConforMix ensemble generated from {APO_PDB} on a GPU Colab "
            f"runtime (see README 'Generation (Colab)'). Skipping 2b; 2a's ground-truth "
            f"geometry above stands alone for now."
        )
        return None

    parser = PDBParser(QUIET=True)
    holo_structure = parser.get_structure("holo", HOLO_PDB)
    apo_structure = parser.get_structure("apo", APO_PDB)
    apo_residues = _standard_ca_residues(apo_structure, CHAIN_ID)

    rows = []
    conformer_holo_rmsds = []
    apo_conformer_rmsds = []
    for conf_path in conformer_paths:
        conf_structure = parser.get_structure("conf", conf_path)
        conf_residues = _standard_ca_residues(conf_structure, CHAIN_ID)
        if len(conf_residues) != len(apo_residues):
            raise ValueError(
                f"{conf_path}: {len(conf_residues)} residues vs apo's {len(apo_residues)} "
                f"-- ConforMix should preserve residue count from the apo input."
            )

        # conf <-> holo, restricted to the same pocket-lining residues (via apo index,
        # since conf shares apo's numbering).
        apo_idx_by_id = {id(r): i for i, r in enumerate(apo_residues)}
        pocket_apo_idx = {apo_idx_by_id[id(a)] for a, h in pocket_pairs}
        holo_by_apo_idx = {apo_idx_by_id[id(a)]: h for a, h in pocket_pairs}

        conf_pocket_ca = np.array([
            conf_residues[i]["CA"].get_coord() for i in sorted(pocket_apo_idx)
        ])
        holo_pocket_ca = np.array([
            holo_by_apo_idx[i]["CA"].get_coord() for i in sorted(pocket_apo_idx)
        ])
        conf_pocket_ca_in_holo_frame = conf_pocket_ca @ rot + tran
        holo_rmsd = float(np.sqrt(np.mean(
            np.sum((conf_pocket_ca_in_holo_frame - holo_pocket_ca) ** 2, axis=1)
        )))
        conformer_holo_rmsds.append(holo_rmsd)

        apo_pocket_ca_ref = np.array([apo_residues[i]["CA"].get_coord() for i in sorted(pocket_apo_idx)])
        apo_rmsd = float(np.sqrt(np.mean(
            np.sum((conf_pocket_ca - apo_pocket_ca_ref) ** 2, axis=1)
        )))
        apo_conformer_rmsds.append(apo_rmsd)

        rows.append({
            "conformer": os.path.basename(conf_path),
            "conformer_holo_rmsd": holo_rmsd,
            "apo_conformer_rmsd": apo_rmsd,
        })
        logger.info(f"  {os.path.basename(conf_path)}: holo_rmsd={holo_rmsd:.3f} A, "
                    f"apo_rmsd={apo_rmsd:.3f} A")

    best_conformer_holo_rmsd = min(conformer_holo_rmsds)
    directional_gain = apo_holo_pocket_ca_rmsd - best_conformer_holo_rmsd

    mean_displacement_from_apo = float(np.mean(apo_conformer_rmsds))
    if mean_displacement_from_apo == 0.0:
        logger.error("Conformers show zero displacement from apo -- ConforMix run failed.")
        null_p = None
    else:
        null_distances = _isotropic_null_best_of_n(
            displacement_magnitude=mean_displacement_from_apo,
            distance_to_target=apo_holo_pocket_ca_rmsd,
            n_conformers=len(conformer_paths),
        )
        null_mean = float(np.mean(null_distances))
        # One-sided: is the observed best-of-N distance to holo smaller than the null
        # distribution's, more often than chance predicts?
        null_p = float(np.mean(null_distances <= best_conformer_holo_rmsd))
        logger.info(f"Isotropic null: mean best-of-{len(conformer_paths)} distance to holo "
                    f"= {null_mean:.3f} A (observed: {best_conformer_holo_rmsd:.3f} A)")

    with open(CONFORMER_GEOMETRY_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["conformer", "conformer_holo_rmsd", "apo_conformer_rmsd"])
        writer.writeheader()
        writer.writerows(rows)

    return {
        "n_conformers": len(conformer_paths),
        "best_conformer_holo_rmsd": best_conformer_holo_rmsd,
        "directional_gain": directional_gain,
        "mean_displacement_from_apo": mean_displacement_from_apo,
        "directional_gain_null_mean": null_mean if mean_displacement_from_apo else None,
        "directional_gain_null_p_one_sided": null_p,
    }


def main() -> int:
    if not (os.path.exists(APO_PDB) and os.path.exists(HOLO_PDB)):
        logger.error(f"{APO_PDB} and/or {HOLO_PDB} not found. Fetch via "
                      f"curl -L https://files.rcsb.org/download/1RTC.pdb -o {APO_PDB} "
                      f"(and 1BR6 similarly).")
        return 1

    geometry, pocket_pairs, rot, tran = compute_apo_holo_geometry()
    logger.info("=== Phase 2a: apo/holo ground truth ===")
    for k, v in geometry.items():
        logger.info(f"  {k}: {v}")

    with open(GEOMETRY_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(geometry.keys()))
        writer.writeheader()
        writer.writerow(geometry)
    logger.info(f"Written -> {GEOMETRY_CSV}")

    logger.info("=== Phase 2b: directional test ===")
    directional = compute_directional_gain(
        pocket_pairs, rot, tran, geometry["apo_holo_pocket_ca_rmsd"]
    )
    if directional:
        for k, v in directional.items():
            logger.info(f"  {k}: {v}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

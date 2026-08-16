"""
Attrition ledger (Phase 7).

The current candidate set has IDs 0000, 0003, 0004, 0006-0009 -- three of ten
generated candidates (0001, 0002, 0005) silently disappeared before docking. Where
candidates die, and why, is a reportable metric in its own right.

Funnel stages actually observable from the payload + local pipeline outputs:

    generated -> rdkit_valid (in valid_candidates.sdf) -> pdbqt_conversion_succeeded
    -> docking_completed -> in_final_results

The plan's "PoseBusters valid" and "3D embedding succeeded" checks aren't separately
logged by the generation pipeline -- the payload only tells us whether a candidate made
it into valid_candidates.sdf or not, not which specific sub-check it failed. For anyone
dropped before that stage, this script attempts a best-effort LOCAL diagnosis from the
final denoising-trajectory frame (fraction of atoms with no neighbor within a plausible
bond length) -- but this is explicitly labeled as inferred, not the original
generation-time error, since a naive RDKit bond-order reconstruction on this data
produces false positives even on candidates we *know* are valid (confirmed by testing
it against mol_0000, which passed generation but still fails naive bond perception).

Usage:
    ./venv/bin/python attrition.py
"""
import glob
import json
import logging
import os
import sys

import numpy as np

from ensemble_auditor import extract_smiles_from_sdf

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

UNPACK_DIR = os.path.join("results", "payload_unpacked")
TRAJECTORIES_DIR = os.path.join(UNPACK_DIR, "valid_trajectories")
SDF_PATH = os.path.join(UNPACK_DIR, "valid_candidates.sdf")
METADATA_PATH = os.path.join(UNPACK_DIR, "metadata.json")
ENSEMBLE_WORK_DIR = os.path.join("results", "ensemble_audit")
RESULTS_JSON = "results.json"
OUTPUT_JSON = "attrition.json"

BOND_DISTANCE_CUTOFF = 2.0  # A -- generous upper bound for a plausible covalent bond
BONDED_FRACTION_THRESHOLD = 0.9  # below this, call it non-convergent


def parse_last_frame(xyz_path: str) -> list[tuple]:
    """The trajectory format has no atom-count header -- each 'Frame N / M | mol_X'
    line starts a new block of 'element x y z' lines. Return the final frame."""
    frames = []
    current = []
    with open(xyz_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("Frame"):
                if current:
                    frames.append(current)
                current = []
            else:
                parts = line.split()
                if len(parts) == 4:
                    current.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    if current:
        frames.append(current)
    if not frames:
        raise ValueError(f"No frames parsed from {xyz_path}")
    return frames[-1]


def diagnose_geometry(xyz_path: str) -> str:
    """Best-effort local reason a candidate never made it into valid_candidates.sdf."""
    try:
        frame = parse_last_frame(xyz_path)
    except Exception as e:
        return f"could not parse trajectory file: {e}"

    coords = np.array([(x, y, z) for _, x, y, z in frame])
    n = len(coords)
    if n < 2:
        return f"final frame has only {n} atom(s)"

    min_dists = np.array([
        np.delete(np.linalg.norm(coords - coords[i], axis=1), i).min()
        for i in range(n)
    ])
    frac_bonded = float((min_dists < BOND_DISTANCE_CUTOFF).mean())

    if frac_bonded < BONDED_FRACTION_THRESHOLD:
        return (
            f"geometry did not converge: only {frac_bonded:.0%} of {n} atoms have a "
            f"neighbor within {BOND_DISTANCE_CUTOFF} A in the final denoising frame "
            f"(max nearest-neighbor distance {min_dists.max():.2f} A)"
        )
    return (
        f"final frame geometry looks locally bonded ({frac_bonded:.0%} of {n} atoms "
        f"have a near neighbor) -- likely failed a downstream chemistry check "
        f"(valence/sanitization/RDKit parse) not preserved in the payload"
    )


def main() -> int:
    if not os.path.exists(TRAJECTORIES_DIR):
        logger.error(f"{TRAJECTORIES_DIR} not found. Run ensemble_auditor.py against a payload first.")
        return 1

    trajectory_paths = sorted(glob.glob(os.path.join(TRAJECTORIES_DIR, "mol_*.xyz")))
    generated_indices = sorted(
        int(os.path.splitext(os.path.basename(p))[0].split("_")[1]) for p in trajectory_paths
    )
    n_generated = len(generated_indices)
    logger.info(f"Stage 'generated': {n_generated} candidates ({TRAJECTORIES_DIR})")

    if os.path.exists(METADATA_PATH):
        with open(METADATA_PATH) as f:
            metadata = json.load(f)
        declared = metadata.get("n_samples_generated")
        if declared is not None and declared != n_generated:
            logger.warning(
                f"metadata.json declares n_samples_generated={declared}, but found "
                f"{n_generated} trajectory files -- using the file count as ground truth."
            )

    sdf_entries = extract_smiles_from_sdf(SDF_PATH) if os.path.exists(SDF_PATH) else []
    rdkit_valid_indices = sorted(
        int(e["OriginalIndex"]) for e in sdf_entries if e.get("OriginalIndex") is not None
    )
    logger.info(f"Stage 'rdkit_valid' (in valid_candidates.sdf): {len(rdkit_valid_indices)} candidates")

    dropped_at_generation = sorted(set(generated_indices) - set(rdkit_valid_indices))
    dropped_reasons = {}
    for idx in dropped_at_generation:
        xyz_path = os.path.join(TRAJECTORIES_DIR, f"mol_{idx:04d}.xyz")
        reason = diagnose_geometry(xyz_path) if os.path.exists(xyz_path) else "trajectory file missing"
        dropped_reasons[f"mol_{idx:04d}"] = reason
        logger.info(f"  ✗ mol_{idx:04d} dropped before rdkit_valid: {reason}")

    # Everything past rdkit_valid: did it reach docking / final results?
    results = []
    if os.path.exists(RESULTS_JSON):
        with open(RESULTS_JSON) as f:
            results = json.load(f)
    docked_smiles = {r["smiles"] for r in results if r.get("overall_best_affinity") is not None}

    sdf_by_index = {int(e["OriginalIndex"]): e for e in sdf_entries if e.get("OriginalIndex") is not None}
    # A candidate that reached docking_engine.run_docking successfully (PDBQT conversion
    # included -- run_docking returns None on a failed prepare_ligand/prepare_receptor
    # step) shows up with a non-null affinity in results.json. This payload draws no
    # distinction between "PDBQT conversion failed" and "Vina itself failed"; both would
    # simply be absent from docked_smiles.
    docking_ok_indices = [
        idx for idx in rdkit_valid_indices if sdf_by_index[idx]["smiles"] in docked_smiles
    ]

    n_rdkit_valid = len(rdkit_valid_indices)
    n_docking_completed = len(docking_ok_indices)
    n_final_results = len(results)

    stages = [
        {"stage": "generated", "count": n_generated, "survival_pct": 100.0},
        {
            "stage": "rdkit_valid",
            "count": n_rdkit_valid,
            "survival_pct": round(100.0 * n_rdkit_valid / n_generated, 1) if n_generated else None,
        },
        {
            "stage": "docking_completed",
            "count": n_docking_completed,
            "survival_pct": round(100.0 * n_docking_completed / n_generated, 1) if n_generated else None,
        },
        {
            "stage": "in_final_results",
            "count": n_final_results,
            "survival_pct": round(100.0 * n_final_results / n_generated, 1) if n_generated else None,
        },
    ]

    report = {
        "n_generated": n_generated,
        "generated_indices": generated_indices,
        "rdkit_valid_indices": rdkit_valid_indices,
        "dropped_before_rdkit_valid": dropped_reasons,
        "stages": stages,
    }
    with open(OUTPUT_JSON, "w") as f:
        json.dump(report, f, indent=2)
    logger.info(f"Attrition ledger written -> {OUTPUT_JSON}")

    logger.info("\n%-22s %8s %12s", "Stage", "Count", "Survival %")
    for s in stages:
        logger.info("%-22s %8d %11s%%", s["stage"], s["count"], s["survival_pct"])

    # Reconciliation check: every generated candidate accounted for exactly once.
    accounted = set(rdkit_valid_indices) | {int(k.split("_")[1]) for k in dropped_reasons}
    assert set(generated_indices) == accounted, (
        f"Reconciliation failure: generated={set(generated_indices)}, "
        f"accounted={accounted}"
    )
    logger.info("Reconciliation OK: every generated candidate accounted for.")

    return 0


if __name__ == "__main__":
    sys.exit(main())

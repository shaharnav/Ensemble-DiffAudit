"""
Noise-corrected rigid control (Phase 4 -- the critical-path experiment).

Taking the best of M conformer scores inflates the result even when all M
conformers are structurally identical, purely because you're maximizing over M
noisy Vina samples. To measure that inflation, this script docks each ligand
against the *unmodified* rigid crystal M times with different seeds and takes
the max. That's the score max-over-M would produce with zero conformational
change. Any real induced-fit effect in ensemble_best_affinity must beat it.

Requires ensemble_auditor.py to have already been run against the payload --
this reuses its aligned baseline-crystal receptor and results.json, and
asserts its own docking parameters (box, exhaustiveness) match what the
ensemble run actually used, so the comparison is apples-to-apples.

Usage:
    ./venv/bin/python rigid_control.py
"""
import argparse
import csv
import glob
import json
import logging
import os
import statistics
import sys

from docking_engine import run_docking

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

ALIGNED_RECEPTORS_DIR = os.path.join("results", "payload_unpacked", "aligned_receptors")
ENSEMBLE_WORK_DIR = os.path.join("results", "ensemble_audit")
RIGID_WORK_DIR = os.path.join("results", "rigid_control")
RESULTS_JSON = "results.json"
OUTPUT_CSV = "rigid_control.csv"

# Same seed list Phase 2's calibration used.
DEFAULT_SEEDS = [1, 2, 3, 4, 5, 6]


def find_baseline_receptor() -> str:
    matches = glob.glob(os.path.join(ALIGNED_RECEPTORS_DIR, "*_baseline_crystal.pdb"))
    if not matches:
        raise RuntimeError(
            f"No baseline crystal receptor found in {ALIGNED_RECEPTORS_DIR}. "
            "Run ensemble_auditor.py against the payload first."
        )
    if len(matches) > 1:
        raise RuntimeError(f"Ambiguous: multiple baseline crystal receptors found: {matches}")
    return matches[0]


def load_ensemble_docking_params(baseline_receptor: str, n_candidates: int) -> dict:
    """
    Read back the center/size/exhaustiveness the ensemble run actually used against
    this baseline receptor (from docking_engine's per-job params sidecars), and
    confirm every ligand's job agreed on them. This is the source of truth the rigid
    control must match -- any drift here invalidates the noise-corrected comparison.

    results/ensemble_audit/ accumulates sidecars across every ensemble_auditor.py run
    ever performed in this repo (including older runs with a different SMILES set or
    pocket center, before job_name index N stops being reused). Only sidecars whose
    index falls within the *current* candidate count are trusted -- anything beyond
    that is cruft from a prior, larger/differently-configured run.
    """
    stem = os.path.splitext(os.path.basename(baseline_receptor))[0]
    sidecars = [
        os.path.join(ENSEMBLE_WORK_DIR, f"ens_s{i:04d}_{stem}_params.json")
        for i in range(n_candidates)
    ]
    sidecars = [p for p in sidecars if os.path.exists(p)]
    if not sidecars:
        raise RuntimeError(
            f"No ensemble docking param sidecars found for '{stem}' in {ENSEMBLE_WORK_DIR}. "
            "Run ensemble_auditor.py against the payload first."
        )

    params_seen = []
    for path in sidecars:
        with open(path) as f:
            params_seen.append(json.load(f))

    reference = params_seen[0]
    for p in params_seen[1:]:
        if (p["center"], p["size"], p["exhaustiveness"]) != (
            reference["center"], reference["size"], reference["exhaustiveness"]
        ):
            raise RuntimeError(
                f"Ensemble run's own docking parameters are inconsistent across ligands "
                f"for '{stem}' -- cannot establish a single ground truth to match. "
                f"{reference} != {p}"
            )

    return {
        "center": reference["center"],
        "size": reference["size"],
        "exhaustiveness": reference["exhaustiveness"],
    }


def summarize_rigid_scores(scores: list[float], ensemble_best_affinity: float | None) -> dict:
    """
    Pure aggregation: given the rigid-control scores for one ligand's M seeded runs
    and that ligand's ensemble_best_affinity, compute rigid_max/mean/sd/range and the
    headline noise_corrected_delta = rigid_max_over_seeds - ensemble_best_affinity.

    Vina affinities are negative-is-better, so "the best score across M seeds" (what
    the plan calls rigid_max_over_seeds) is the MINIMUM value, not Python's max().
    """
    rigid_max = min(scores)
    rigid_mean = statistics.mean(scores)
    rigid_sd = statistics.stdev(scores) if len(scores) > 1 else 0.0
    rigid_range = max(scores) - min(scores)
    noise_corrected_delta = (
        rigid_max - ensemble_best_affinity if ensemble_best_affinity is not None else None
    )
    return {
        "rigid_max": rigid_max,
        "rigid_mean": rigid_mean,
        "rigid_sd": rigid_sd,
        "rigid_range": rigid_range,
        "noise_corrected_delta": noise_corrected_delta,
    }


def load_candidates() -> list[dict]:
    if not os.path.exists(RESULTS_JSON):
        raise RuntimeError(f"{RESULTS_JSON} not found. Run ensemble_auditor.py first.")
    with open(RESULTS_JSON) as f:
        return json.load(f)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--seeds", type=int, nargs="+", default=None,
        help="Vina seeds to use (default: as many of [1,2,3,4,5,6] as there are conformers)",
    )
    args = parser.parse_args()

    baseline_receptor = find_baseline_receptor()
    logger.info(f"Rigid control receptor: {baseline_receptor}")

    candidates = load_candidates()
    logger.info(f"Loaded {len(candidates)} candidates from {RESULTS_JSON}")

    docking_params = load_ensemble_docking_params(baseline_receptor, len(candidates))
    logger.info(
        f"Matched ensemble run's docking parameters: center={docking_params['center']}, "
        f"size={docking_params['size']}, exhaustiveness={docking_params['exhaustiveness']}"
    )

    # M = number of conformers in the ensemble, so max-over-M is measured on equal footing.
    conformer_count = None
    for c in candidates:
        n = sum(1 for k in c.get("all_affinities", {}) if "baseline_crystal.pdb" not in k)
        if n:
            conformer_count = n
            break
    if conformer_count is None:
        raise RuntimeError("Could not determine ensemble conformer count from results.json")

    seeds = args.seeds if args.seeds is not None else DEFAULT_SEEDS[:conformer_count]
    if len(seeds) != conformer_count:
        logger.warning(
            f"Using {len(seeds)} seeds but the ensemble had {conformer_count} conformers -- "
            "max-over-M comparison will not be on equal footing unless this is intentional."
        )
    logger.info(f"M = {conformer_count} conformers -> seeds: {seeds}")

    os.makedirs(RIGID_WORK_DIR, exist_ok=True)

    rows = []
    for idx, cand in enumerate(candidates):
        smiles = cand["smiles"]
        cand_id = cand.get("id", f"idx{idx}")
        ensemble_best_affinity = cand.get("ensemble_best_affinity")

        logger.info(f"[{idx + 1}/{len(candidates)}] {cand_id}: {smiles[:60]}")

        scores = []
        for seed in seeds:
            job_name = f"rigid_{cand_id}_seed{seed}"
            result = run_docking(
                pdb_file=baseline_receptor,
                smiles=smiles,
                output_dir=RIGID_WORK_DIR,
                job_name=job_name,
                exhaustiveness=docking_params["exhaustiveness"],
                center_coords=docking_params["center"],
                box_size=docking_params["size"],
                seed=seed,
            )
            if result is None or result.get("affinity") is None:
                logger.warning(f"  ✗ Docking failed for seed {seed}.")
                continue
            affinity = result["affinity"]
            scores.append(affinity)
            logger.info(f"  seed={seed}: affinity={affinity:.3f} kcal/mol")

        if not scores:
            logger.warning(f"  All rigid-control docking runs failed for {cand_id}.")
            summary = {"rigid_max": None, "rigid_mean": None, "rigid_sd": None,
                       "rigid_range": None, "noise_corrected_delta": None}
        else:
            summary = summarize_rigid_scores(scores, ensemble_best_affinity)

        rows.append({
            "id": cand_id,
            "smiles": smiles,
            "seeds": seeds,
            "rigid_scores": scores,
            "ensemble_best_affinity": ensemble_best_affinity,
            **summary,
        })

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "ID", "SMILES", "seeds", "rigid_scores", "rigid_max", "rigid_mean",
            "rigid_sd", "rigid_range", "ensemble_best_affinity", "noise_corrected_delta",
        ])
        for r in rows:
            writer.writerow([
                r["id"], r["smiles"], json.dumps(r["seeds"]), json.dumps(r["rigid_scores"]),
                r["rigid_max"], r["rigid_mean"], r["rigid_sd"], r["rigid_range"],
                r["ensemble_best_affinity"], r["noise_corrected_delta"],
            ])
    logger.info(f"Rigid control results written -> {OUTPUT_CSV}")

    valid_deltas = [r["noise_corrected_delta"] for r in rows if r["noise_corrected_delta"] is not None]
    if valid_deltas:
        logger.info(
            f"noise_corrected_delta across {len(valid_deltas)} candidates: "
            f"mean={statistics.mean(valid_deltas):.3f}, "
            f"median={statistics.median(valid_deltas):.3f}"
        )
        beat_noise = sum(1 for d in valid_deltas if d > 0)
        logger.info(
            f"{beat_noise}/{len(valid_deltas)} candidates beat the noise-matched rigid "
            f"max-over-{conformer_count} baseline."
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())

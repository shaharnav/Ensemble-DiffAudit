"""
Experiment 3, Phase 2: pocket characterization and the volume gate.

No fpocket binary is available in this environment (checked: not on PATH, not in
Homebrew), so this uses a grid-based proxy instead, per the plan's explicit fallback.

Center is read from pocket_config.yaml (the ligand-transfer-based box, not a geometric cavity
finder -- see pocket_config.py) rather than recomputed here, per the rule that every
downstream stage reads center/box from that one file.

Method: a fixed 3D grid (0.5 A spacing) inside a sphere of radius `POCKET_RADIUS` centered on
pocket_config.yaml's center. `cavity_volume_A3` is the total non-clashing volume inside that
sphere -- grid points farther than `CLASH_MIN` (a water-probe-radius cutoff) from every
receptor heavy atom.

**This is a coarser proxy than true SASA-based cavity detection (e.g. fpocket): it does not
do a flood-fill enclosure test, so it measures "how much empty space sits within a fixed
radius of the catalytic triad" rather than the volume of a topologically enclosed pocket.**
An open groove and a genuinely enclosed pocket of the same local empty-space extent would
score similarly. It is used here only for relative comparison across R1/R2/conformers of the
*same* protein with the *same* sphere, which is what Phase 2's gate and Phase 3's
per-conformer tracking need -- not as an absolute physical cavity volume. This limitation is
recorded, not hidden.

Usage:
    ./venv/bin/python pocket_volume.py
"""
import json
import logging
import os
import sys

import numpy as np
import yaml
from Bio.PDB import PDBParser

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

R1_PATH = os.path.join("results", "experiment3", "receptors", "R1_zymogen.pdb")
R2_PATH = os.path.join("results", "experiment3", "receptors", "R2_mature_proxy.pdb")
POCKET_CONFIG_YAML = "pocket_config.yaml"
CHAIN_ID = "A"
POCKET_RADIUS = 10.0  # sphere radius for the cavity-volume estimate
GRID_SPACING = 0.5
CLASH_MIN = 1.4
OUTPUT_JSON = os.path.join("results", "experiment3", "pocket_volume.json")


def get_active_site_centroid() -> tuple[np.ndarray, list]:
    """Reads center/box from pocket_config.yaml -- the ligand-transfer-based
    determination, not a geometric cavity finder. Never recomputed here."""
    with open(POCKET_CONFIG_YAML) as f:
        config = yaml.safe_load(f)
    return np.array(config["center"]), config["box_size_A"]


def compute_cavity_volume(pdb_path: str, centroid: np.ndarray) -> dict:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("r", pdb_path)
    atoms = np.array([
        a.get_coord() for a in structure.get_atoms()
        if a.element != "H"
    ])

    # Only atoms that could plausibly reach into the sphere -- cheap pre-filter.
    nearby = atoms[np.linalg.norm(atoms - centroid, axis=1) <= POCKET_RADIUS + CLASH_MIN]

    axis = np.arange(-POCKET_RADIUS, POCKET_RADIUS + 1e-6, GRID_SPACING)
    n_open = 0
    n_clash = 0
    n_outside_sphere = 0
    for x in axis:
        for y in axis:
            z_offsets = axis
            pts = np.stack([
                np.full_like(z_offsets, centroid[0] + x),
                np.full_like(z_offsets, centroid[1] + y),
                centroid[2] + z_offsets,
            ], axis=1)
            in_sphere = np.linalg.norm(pts - centroid, axis=1) <= POCKET_RADIUS
            if not np.any(in_sphere):
                continue
            pts = pts[in_sphere]
            d = np.linalg.norm(nearby[None, :, :] - pts[:, None, :], axis=2).min(axis=1)
            n_clash += int(np.sum(d < CLASH_MIN))
            n_open += int(np.sum(d >= CLASH_MIN))

    voxel_volume = GRID_SPACING ** 3
    return {
        "cavity_volume_A3": round(n_open * voxel_volume, 1),
        "n_open_voxels": n_open,
        "n_clash_voxels": n_clash,
    }


def main() -> int:
    for path in (R1_PATH, R2_PATH, POCKET_CONFIG_YAML):
        if not os.path.exists(path):
            logger.error(f"{path} not found -- run exp3_receptor_prep.py / pocket_config.py first.")
            return 1

    centroid, box_size = get_active_site_centroid()
    logger.info(f"Active-site centroid (from pocket_config.yaml): {centroid}")
    logger.info(f"Docking box: center {centroid.tolist()}, size {box_size} A")

    results = {}
    for label, path in (("R1_zymogen", R1_PATH), ("R2_mature_proxy", R2_PATH)):
        vol = compute_cavity_volume(path, centroid)
        results[label] = vol
        logger.info(f"{label}: cavity_volume_A3={vol['cavity_volume_A3']} "
                    f"(open voxels={vol['n_open_voxels']}, clash={vol['n_clash_voxels']})")

    r1_vol = results["R1_zymogen"]["cavity_volume_A3"]
    r2_vol = results["R2_mature_proxy"]["cavity_volume_A3"]
    logger.info(f"R1 -> R2 cavity opening: {r1_vol} -> {r2_vol} A3 "
                f"({'+' if r2_vol >= r1_vol else ''}{r2_vol - r1_vol:.1f} A3)")

    if r1_vol < 150:
        logger.warning(f"R1 cavity volume ({r1_vol} A3) is below the ~150 A3 threshold for "
                        f"a drug-sized ligand -- R1 docking scores will likely be dominated "
                        f"by clashes. Proceeding per the plan, but expect and report a "
                        f"large R1->R2 gap.")
    if r2_vol < 150:
        logger.error(f"R2 cavity volume ({r2_vol} A3) is ALSO below ~150 A3 -- "
                      f"reconsider the target. Stopping.")
        return 1

    output = {
        "box_center": centroid.tolist(),
        "box_size_A": box_size,
        "box_source": POCKET_CONFIG_YAML,
        "pocket_radius_A": POCKET_RADIUS,
        "grid_spacing_A": GRID_SPACING,
        "clash_min_A": CLASH_MIN,
        **results,
    }
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)
    logger.info(f"Written -> {OUTPUT_JSON}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

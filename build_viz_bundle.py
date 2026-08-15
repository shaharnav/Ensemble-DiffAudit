"""
Builds results/payload_unpacked/viz_bundle.json for the 4D denoising-trajectory frontend viewer.

Combines:
  - the DiffSBDD denoising trajectories (results/payload_unpacked/valid_trajectories/*.xyz)
  - the final candidate geometry + bonds (results/payload_unpacked/valid_candidates.sdf)
  - a local Vina docking pass (via docking_engine.run_docking) for the 3 candidates whose
    trajectory has been geometry-verified against the SDF (see plan notes)
  - the 6-conformer receptor ensemble (results/payload_unpacked/ensemble_receptors/)

Run with the project venv: venv/bin/python3 build_viz_bundle.py
"""
import glob
import json
import os
import re

import numpy as np
from rdkit import Chem
from scipy.optimize import linear_sum_assignment

from docking_engine import run_docking

UNPACKED = "results/payload_unpacked"
TRAJ_DIR = os.path.join(UNPACKED, "valid_trajectories")
SDF_PATH = os.path.join(UNPACKED, "valid_candidates.sdf")
POCKET_PATH = os.path.join(UNPACKED, "pocket.pdb")
# Prefer the ConforMix variants that ensemble_auditor.py has sequence-aligned to the
# reference crystal frame (see its "Align ConforMix variants" step) over the raw,
# unaligned ensemble_receptors/ copies -- pocket-lining and spatial calculations here
# are only meaningful in the aligned frame. Falls back to the raw copies if the
# auditor hasn't been run yet.
_ALIGNED_RECEPTOR_DIR = os.path.join(UNPACKED, "aligned_receptors")
RECEPTOR_DIR = (
    _ALIGNED_RECEPTOR_DIR
    if glob.glob(os.path.join(_ALIGNED_RECEPTOR_DIR, "conformix_var_*.pdb"))
    else os.path.join(UNPACKED, "ensemble_receptors")
)
METADATA_PATH = os.path.join(UNPACKED, "metadata.json")
OUT_PATH = os.path.join(UNPACKED, "viz_bundle.json")
RECEPTOR_PDB_OUT = os.path.join(UNPACKED, "receptor_breathing.pdb")
DOCK_OUTPUT_DIR = os.path.join(UNPACKED, "traj_docking")
POCKET_LINING_RADIUS = 8.0  # Angstrom, for highlighting pocket-lining residues in the viewer

# Geometry-verified mapping: SDF <OriginalIndex> -> trajectory filename.
# Verified via ICP correspondence + Kabsch RMSD (see icp_correspondence below)
# against every file in valid_trajectories/; true matches score ~0.0005 A, all
# others fail on atom-count alone. The other 4 SDF candidates (of 7 valid) have
# no matching trajectory (likely dropped during the notebook's RDKit
# sanitize/fragment-pick step) and are excluded here.
ORIGIDX_TO_TRAJ = {
    3: "mol_0003.xyz",
    4: "mol_0004.xyz",
    6: "mol_0006.xyz",
}


def parse_sdf(path):
    """Parse valid_candidates.sdf without relying on RDKit's SDF bond perception,
    so we keep the exact atom order/bonds as written (matches the .xyz atom order)."""
    text = open(path).read()
    blocks = [b for b in text.split("$$$$\n") if b.strip()]
    mols = []
    for b in blocks:
        lines = b.split("\n")
        counts_line = lines[3]
        natoms = int(counts_line[0:3])
        nbonds = int(counts_line[3:6])
        elems, coords = [], []
        for l in lines[4 : 4 + natoms]:
            parts = l.split()
            coords.append([float(parts[0]), float(parts[1]), float(parts[2])])
            elems.append(parts[3])
        bonds = []
        for l in lines[4 + natoms : 4 + natoms + nbonds]:
            i, j, order = int(l[0:3]), int(l[3:6]), int(l[6:9])
            bonds.append([i - 1, j - 1, order])  # to 0-indexed
        # SD tag blocks look like: >  <TAG>  (n) \nVALUE\n\n
        props = {}
        for tag in ("OriginalIndex", "QED", "SA_Score"):
            m = re.search(rf"<{tag}>[^\n]*\n([^\n]+)", b)
            if m:
                props[tag] = m.group(1).strip()
        mols.append(
            dict(
                natoms=natoms,
                elements=elems,
                coords=np.array(coords),
                bonds=bonds,
                original_index=int(props["OriginalIndex"]),
                qed=float(props["QED"]),
                sa_score=float(props["SA_Score"]),
            )
        )
    return mols


def read_xyz_trajectory(path):
    lines = open(path).read().strip().split("\n")
    n = int(lines[0])
    frame_size = n + 2
    nframes = len(lines) // frame_size
    frames = np.zeros((nframes, n, 3))
    elements = []
    for f in range(nframes):
        start = f * frame_size
        frame_elems = []
        for a, l in enumerate(lines[start + 2 : start + 2 + n]):
            parts = l.split()
            frame_elems.append(parts[0])
            frames[f, a] = [float(parts[1]), float(parts[2]), float(parts[3])]
        elements.append(frame_elems)
    return frames, elements


def kabsch(mobile, target):
    """Returns (R, t) such that mobile @ R.T + t ~= target (both Nx3, matched order)."""
    mobile_c = mobile - mobile.mean(axis=0)
    target_c = target - target.mean(axis=0)
    H = mobile_c.T @ target_c
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = target.mean(axis=0) - R @ mobile.mean(axis=0)
    return R, t


def icp_correspondence(elems_a, coords_a, elems_b, coords_b, iters=20):
    """The .xyz trajectory and the SDF list the same atoms in different orders (confirmed by
    inspection). Solve atom correspondence via ICP: alternate a per-element Hungarian assignment
    (scipy linear_sum_assignment) with a Kabsch re-fit, until the assignment stabilizes.

    Returns (perm, R, t, rmsd) where perm[j] is the index into `a` for sdf-position j, i.e.
    coords_a[perm] @ R.T + t ~= coords_b.
    """
    n = len(elems_a)
    R, t = np.eye(3), np.zeros(3)
    row = col = None
    rmsd = None
    for _ in range(iters):
        aligned = coords_a @ R.T + t
        cost = np.full((n, n), 1e6)
        for i in range(n):
            for j in range(n):
                if elems_a[i] == elems_b[j]:
                    cost[i, j] = np.linalg.norm(aligned[i] - coords_b[j])
        row, col = linear_sum_assignment(cost)
        R, t = kabsch(coords_a[row], coords_b[col])
        rmsd = np.sqrt(np.mean(np.sum((coords_a[row] @ R.T + t - coords_b[col]) ** 2, axis=1)))
    perm = row[np.argsort(col)]  # perm[j] = trajectory atom index for sdf position j
    return perm, R, t, rmsd


def parse_pdb_atom_lines(path):
    """Return the raw ATOM/HETATM lines (water stripped) plus parsed (chain, resi, coords)
    per line, so the same file can be reused both for geometry (pocket-lining residues) and
    for writing out real PDB text (3Dmol renders cartoon from PDB records directly)."""
    lines, chains, resis, coords = [], [], [], []
    for line in open(path):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname in ("HOH", "WAT"):
                continue
            lines.append(line)
            chains.append(line[21])
            resis.append(int(line[22:26]))
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return lines, chains, resis, np.array(coords)


def pocket_lining_residues(path, pocket_center, radius):
    """Residues (chain, resi) with at least one atom within `radius` of pocket_center,
    in this PDB's own numbering (ConforMix renumbers from 1, unrelated to pocket.pdb)."""
    _, chains, resis, coords = parse_pdb_atom_lines(path)
    dists = np.linalg.norm(coords - np.array(pocket_center), axis=1)
    close = dists <= radius
    return sorted({(chains[i], resis[i]) for i in range(len(resis)) if close[i]})


def write_multimodel_pdb(variant_files, out_path):
    """Concatenate the (index-aligned, same-topology) receptor conformers as MODEL/ENDMDL
    records so 3Dmol can load them with addModelAsFrames and animate between them."""
    with open(out_path, "w") as out:
        for i, vf in enumerate(variant_files, start=1):
            lines, *_ = parse_pdb_atom_lines(vf)
            out.write(f"MODEL     {i:>4}\n")
            out.writelines(lines)
            out.write("ENDMDL\n")
        out.write("END\n")


def _build_baseline_receptor(pdb_id):
    """Strip waters + native ligand from the full reference crystal structure
    (pdbs/<pdb_id>.pdb) to use as the demo-candidate docking target, instead of the
    ~8-residue pocket.pdb crop (too few atoms for a meaningful docking score)."""
    from Bio.PDB import PDBParser, PDBIO, Select

    reference_pdb = os.path.join("pdbs", f"{pdb_id}.pdb")
    if not os.path.exists(reference_pdb):
        print(f"[warn] {reference_pdb} not found; falling back to {POCKET_PATH}")
        return POCKET_PATH

    class _StripHetero(Select):
        def accept_residue(self, residue):
            resname = residue.get_resname().strip().upper()
            return resname not in ("HOH", "WAT", "BEN")

    baseline_path = os.path.join(UNPACKED, f"{pdb_id}_baseline_crystal.pdb")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ref_baseline", reference_pdb)
    io = PDBIO()
    io.set_structure(structure)
    io.save(baseline_path, _StripHetero())
    return baseline_path


def main():
    metadata = json.load(open(METADATA_PATH))
    pocket_center = metadata["pocket_center"]
    pocket_radius = metadata["pocket_radius"]
    box_size = [2 * pocket_radius] * 3
    baseline_receptor = _build_baseline_receptor(metadata["pdb_id"])

    sdf_mols = {m["original_index"]: m for m in parse_sdf(SDF_PATH)}

    os.makedirs(DOCK_OUTPUT_DIR, exist_ok=True)

    candidates = []
    for origidx, traj_file in ORIGIDX_TO_TRAJ.items():
        mol = sdf_mols[origidx]
        frames, elements = read_xyz_trajectory(os.path.join(TRAJ_DIR, traj_file))

        # The .xyz and SDF list the same atoms in different orders. Solve the correspondence
        # (per-element Hungarian assignment) and rigid alignment (Kabsch) from the last frame,
        # then apply both to every frame so atom index j is consistent with the SDF/bonds and
        # coordinates land in the pocket frame throughout the trajectory.
        perm, R, t, rmsd = icp_correspondence(elements[-1], frames[-1], mol["elements"], mol["coords"])
        print(f"[origidx {origidx}] correspondence+alignment RMSD to SDF final frame: {rmsd:.5f} A")

        aligned_frames = frames[:, perm, :] @ R.T + t
        elements = [[frame_elems[i] for i in perm] for frame_elems in elements]

        smiles = Chem.MolToSmiles(Chem.MolFromMolBlock(_molblock_for(mol)))
        print(f"[origidx {origidx}] SMILES: {smiles}")

        dock = run_docking(
            baseline_receptor,
            smiles,
            output_dir=DOCK_OUTPUT_DIR,
            job_name=f"trypsin_traj_{origidx}",
            exhaustiveness=16,
            center_coords=pocket_center,
            box_size=box_size,
        )
        affinity = dock["affinity"] if dock else None
        print(f"[origidx {origidx}] Vina affinity: {affinity}")

        candidates.append(
            dict(
                original_index=origidx,
                smiles=smiles,
                affinity=affinity,
                qed=mol["qed"],
                sa_score=mol["sa_score"],
                alignment_rmsd=round(float(rmsd), 4),
                frames=np.round(aligned_frames, 3).tolist(),
                elements=elements,
                bonds=mol["bonds"],
            )
        )

    scored = [c for c in candidates if c["affinity"] is not None]
    if scored:
        best = min(scored, key=lambda c: c["affinity"])
        for c in candidates:
            c["default"] = c is best
    else:
        candidates[0]["default"] = True
        for c in candidates[1:]:
            c["default"] = False

    # Receptor ensemble: write as a real multi-model PDB (3Dmol renders cartoon straight from
    # PDB records/residue info) instead of a bare coordinate array.
    variant_files = sorted(glob.glob(os.path.join(RECEPTOR_DIR, "conformix_var_*.pdb")))
    write_multimodel_pdb(variant_files, RECEPTOR_PDB_OUT)
    lining = pocket_lining_residues(variant_files[0], pocket_center, POCKET_LINING_RADIUS)
    print(f"Receptor ensemble: {len(variant_files)} variants -> {RECEPTOR_PDB_OUT}; "
          f"{len(lining)} pocket-lining residues within {POCKET_LINING_RADIUS} A")

    bundle = dict(
        meta=dict(
            pdb_id=metadata["pdb_id"],
            pocket_center=pocket_center,
            pocket_radius=pocket_radius,
            receptor_pdb=os.path.basename(RECEPTOR_PDB_OUT),
            pocket_lining_residues=[dict(chain=c, resi=r) for c, r in lining],
        ),
        candidates=candidates,
    )

    with open(OUT_PATH, "w") as f:
        json.dump(bundle, f)
    size_kb = os.path.getsize(OUT_PATH) / 1024
    print(f"Wrote {OUT_PATH} ({size_kb:.1f} KB)")


def _molblock_for(mol):
    """Reconstruct a minimal V2000 molblock from parsed atoms/bonds so RDKit can derive SMILES."""
    lines = ["", "  Built", "", f"{mol['natoms']:>3}{len(mol['bonds']):>3}  0  0  0  0  0  0  0  0999 V2000"]
    for elem, (x, y, z) in zip(mol["elements"], mol["coords"]):
        lines.append(f"{x:>10.4f}{y:>10.4f}{z:>10.4f} {elem:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for i, j, order in mol["bonds"]:
        lines.append(f"{i+1:>3}{j+1:>3}{order:>3}  0")
    lines.append("M  END")
    return "\n".join(lines)


if __name__ == "__main__":
    main()

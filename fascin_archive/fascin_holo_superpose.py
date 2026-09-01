"""
Fascin Step 1b/1c: superpose each of the 11 fascin-inhibitor holo structures
(resolution <=2.5A, UniProt Q16658) onto the 3LLP (chain A) apo frame via CA
atoms, transform each ligand into that frame, report per-pair alignment RMSD,
occupancy/B-factor, and write the final ligand set table + reference-pose PDBs.

Ligand-in-pocket check: confirm each ligand's centroid falls within the same
pocket contact shell established during target screening (3LLP->6I11 access
residues), the fascin analogue of the LasB "distance to catalytic Zn" check --
fascin has no catalytic metal, so proximity to the already-verified pocket
residue set is the site-consistency criterion here.
"""
import gemmi
import numpy as np
import csv
import os

APO_PATH = "results/target_screen/structures/3LLP.pdb"
APO_CHAIN = "A"
STRUCT_DIR = "results/target_screen/structures"
OUT_DIR = "results/fascin_ensemble_rmsd/holo_ligands_aligned"
os.makedirs(OUT_DIR, exist_ok=True)

STRUCTS = {
    "6I0Z": ("GZQ", "B", "pdb"), "6I10": ("GZK", "A", "pdb"), "6I11": ("H0H", "A", "pdb"),
    "6I12": ("H08", "A", "pdb"), "6I13": ("H0Q", "A", "pdb"), "6I14": ("GZN", "A", "pdb"),
    "6I15": ("GZT", "A", "pdb"), "6I16": ("H0B", "A", "pdb"), "6I17": ("GZW", "A", "pdb"),
    "6I18": ("H0N", "A", "pdb"), "9GS6": ("A1IPD", "A", "cif"),
}


def kabsch(mob, ref):
    mob_c, ref_c = mob.mean(axis=0), ref.mean(axis=0)
    H = (mob - mob_c).T @ (ref - ref_c)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = ref_c - R @ mob_c
    return R, t


def get_ca_dict(chain):
    out = {}
    for res in chain:
        atom = next((a for a in res if a.name == "CA"), None)
        if atom is not None:
            out[res.seqid.num] = atom.pos
    return out


apo_st = gemmi.read_structure(APO_PATH)
apo_st.setup_entities()
apo_chain = apo_st[0][APO_CHAIN]
apo_ca = get_ca_dict(apo_chain)

# established pocket residues from the target-screening pass (3LLP(A)->6I11
# access-shell), used as the site-consistency reference for every ligand
from cryptic_pocket_screen import get_chain_residues, align_residues
holo_ref_st = gemmi.read_structure(f"{STRUCT_DIR}/6I11.pdb")
holo_ref_st.setup_entities()
_matched_ref = align_residues(get_chain_residues(apo_chain), get_chain_residues(holo_ref_st[0]["A"]))
_ref_H0H = [a for r in holo_ref_st[0]["A"] if r.name == "H0H" for a in r]
_ref_H0H_coords = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in _ref_H0H])
POCKET_RESNUMS_HOLO_FRAME = set()
for apo_rn, apo_res, holo_rn, holo_res in _matched_ref:
    rc = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in holo_res if a.name != "H"])
    if rc.size and np.min(np.linalg.norm(_ref_H0H_coords[:, None, :] - rc[None, :, :], axis=2)) < 5.0:
        POCKET_RESNUMS_HOLO_FRAME.add(apo_rn)  # store as apo-frame resnum (consistent numbering here)
print(f"Reference pocket (from 6I11/H0H): {len(POCKET_RESNUMS_HOLO_FRAME)} residues -> {sorted(POCKET_RESNUMS_HOLO_FRAME)}")

rows = []
for pdbid, (ligcode, holo_chain_id, ext) in STRUCTS.items():
    st = gemmi.read_structure(f"{STRUCT_DIR}/{pdbid}.{ext}")
    st.setup_entities()
    model = st[0]
    holo_chain = model[holo_chain_id]

    lig_res = next((r for r in holo_chain if r.name == ligcode), None)
    if lig_res is None:
        rows.append({"pdbid": pdbid, "ligand": ligcode, "error": "ligand_not_found"})
        continue
    atoms = [a for a in lig_res if a.altloc in ("", "A", "\x00")]
    lig_atoms = atoms if atoms else list(lig_res)

    apo_res_list = get_chain_residues(apo_chain)
    holo_res_list = get_chain_residues(holo_chain)
    matched = align_residues(apo_res_list, holo_res_list)
    if len(matched) < 30:
        rows.append({"pdbid": pdbid, "ligand": ligcode, "error": "too_few_matched_residues"})
        continue

    apo_ca_arr = np.array([[next(a for a in r[1] if a.name == "CA").pos.x,
                             next(a for a in r[1] if a.name == "CA").pos.y,
                             next(a for a in r[1] if a.name == "CA").pos.z] for r in matched])
    holo_ca_arr = np.array([[next(a for a in r[3] if a.name == "CA").pos.x,
                              next(a for a in r[3] if a.name == "CA").pos.y,
                              next(a for a in r[3] if a.name == "CA").pos.z] for r in matched])
    R, t = kabsch(holo_ca_arr, apo_ca_arr)
    holo_ca_aligned = (R @ holo_ca_arr.T).T + t
    align_rmsd = float(np.sqrt(np.mean(np.sum((holo_ca_aligned - apo_ca_arr) ** 2, axis=1))))

    lig_coords = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lig_atoms])
    lig_aligned = (R @ lig_coords.T).T + t

    # site-consistency: does this ligand's centroid land near the established pocket?
    apo_pocket_ca = np.array([[apo_ca[rn].x, apo_ca[rn].y, apo_ca[rn].z]
                               for rn in POCKET_RESNUMS_HOLO_FRAME if rn in apo_ca])
    centroid = lig_aligned.mean(axis=0)
    min_dist_to_pocket = float(np.min(np.linalg.norm(apo_pocket_ca - centroid, axis=1))) if len(apo_pocket_ca) else None

    occupancies = [a.occ for a in lig_atoms]
    bfactors = [a.b_iso for a in lig_atoms]

    out_pdb = f"{OUT_DIR}/{pdbid}_{ligcode}_ref_pose.pdb"
    with open(out_pdb, "w") as f:
        for i, (atom, xyz) in enumerate(zip(lig_atoms, lig_aligned)):
            elem = atom.element.name.upper()
            f.write(f"HETATM{i+1:5d} {atom.name:<4s} LIG L 900    "
                    f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{atom.occ:6.2f}{atom.b_iso:6.2f}"
                    f"          {elem:>2s}\n")
        f.write("END\n")

    rows.append({
        "pdbid": pdbid, "ligand": ligcode, "n_ca_matched": len(matched),
        "ca_alignment_rmsd_A": round(align_rmsd, 3),
        "min_dist_to_established_pocket_A": round(min_dist_to_pocket, 2) if min_dist_to_pocket is not None else "",
        "n_lig_atoms": len(lig_atoms),
        "mean_occupancy": round(float(np.mean(occupancies)), 2),
        "min_occupancy": round(float(np.min(occupancies)), 2),
        "mean_bfactor": round(float(np.mean(bfactors)), 1),
        "max_bfactor": round(float(np.max(bfactors)), 1),
        "reference_pose_pdb": out_pdb,
    })
    print(f"{pdbid} ({ligcode}): {len(matched)} CA matched, align RMSD={align_rmsd:.3f}A, "
          f"dist-to-pocket={rows[-1]['min_dist_to_established_pocket_A']}A, "
          f"occ[min/mean]={rows[-1]['min_occupancy']}/{rows[-1]['mean_occupancy']}, "
          f"B[mean/max]={rows[-1]['mean_bfactor']}/{rows[-1]['max_bfactor']}")

with open("results/fascin_ensemble_rmsd/holo_ligand_set_final.csv", "w", newline="") as f:
    fieldnames = sorted(set().union(*[r.keys() for r in rows]))
    w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
    w.writeheader()
    w.writerows(rows)
print(f"\n{len(rows)} ligands processed. Written results/fascin_ensemble_rmsd/holo_ligand_set_final.csv")

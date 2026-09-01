"""
Step 2c/2d: superpose each holo LasB structure onto the 1EZM frame via CA
atoms (residue-number matched), transform ligand coordinates into that frame,
and write the reference-pose ligand SDF/PDB plus the final ligand-set table.

Occupancy note: V85 (7OC7) and IEV (7Z68) are excluded -- both modeled at
partial occupancy (0.5/0.38), so "the" reference pose is ambiguous (alt-conf
ligand). 7QH1 (CI8) is kept but flagged: full occupancy, but high B-factor
(mean 50.6, max 105.3) and the largest/most flexible ligand in the set.
"""
import gemmi
import numpy as np
import csv
import os

STRUCTS = {
    "1U4G": "HPI", "3DBK": "RDF", "6F8B": "CXH", "6FZX": "EEK",
    "7QH1": "CI8", "8R1B": "XI5", "8CC4": "U7F",
    "9FQD": "A1IEX", "9FQE": "A1IEU", "9FQX": "A1IJJ", "9FQY": "A1IFN",
    "9FRY": "A1IFU", "9FRZ": "A1IFZ", "9GM4": "A1IM0", "9FS0": "A1IFY",
}
EXCLUDED = {"7OC7": "V85 modeled at 0.5/0.38 occupancy (alt-conf) -- reference pose ambiguous",
            "7Z68": "IEV modeled at 0.5 occupancy (alt-conf) -- reference pose ambiguous"}

REF_PATH = "results/lasb_payload/ensemble_receptors_aligned/1EZM_baseline_crystal.pdb"
ref_st = gemmi.read_structure(REF_PATH)
ref_st.setup_entities()
ref_chain = ref_st[0][0]
ref_ca = {res.seqid.num: res[0].pos for res in ref_chain
          if any(a.name == "CA" for a in res) for a in [next(a for a in res if a.name == "CA")]}

out_dir = "results/lasb_ensemble_rmsd/holo_ligands_aligned"
os.makedirs(out_dir, exist_ok=True)

rows = []
for pdbid, ligcode in STRUCTS.items():
    for ext in ["pdb", "cif"]:
        path = f"results/lasb_ensemble_rmsd/holo_structures/{pdbid}.{ext}"
        if os.path.exists(path):
            break
    st = gemmi.read_structure(path)
    st.setup_entities()
    model = st[0]

    # Multiple copies of the complex can exist in the asymmetric unit --
    # find the specific chain that carries this ligand, and match CA atoms
    # only within that same chain (not blended across copies).
    lig_chain = None
    for chain in model:
        for res in chain:
            if res.name == ligcode:
                lig_chain = chain
                break
        if lig_chain is not None:
            break

    lig_atoms = [a for res in lig_chain for a in res if res.name == ligcode]

    mob_ca, ref_ca_matched = [], []
    for res in lig_chain:
        atom = next((a for a in res if a.name == "CA"), None)
        if atom is not None and res.seqid.num in ref_ca:
            mob_ca.append([atom.pos.x, atom.pos.y, atom.pos.z])
            p = ref_ca[res.seqid.num]
            ref_ca_matched.append([p.x, p.y, p.z])

    mob_ca = np.array(mob_ca)
    ref_ca_matched = np.array(ref_ca_matched)
    n = len(mob_ca)

    # Kabsch superposition: mobile (holo) -> reference (1EZM) frame
    mob_c = mob_ca.mean(axis=0)
    ref_c = ref_ca_matched.mean(axis=0)
    H = (mob_ca - mob_c).T @ (ref_ca_matched - ref_c)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = ref_c - R @ mob_c

    mob_aligned = (R @ mob_ca.T).T + t
    align_rmsd = np.sqrt(np.mean(np.sum((mob_aligned - ref_ca_matched) ** 2, axis=1)))

    # Transform ligand atoms into the 1EZM frame
    lig_coords = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lig_atoms])
    lig_aligned = (R @ lig_coords.T).T + t

    out_pdb = f"{out_dir}/{pdbid}_{ligcode}_ref_pose.pdb"
    with open(out_pdb, "w") as f:
        for i, (atom, xyz) in enumerate(zip(lig_atoms, lig_aligned)):
            elem = atom.element.name.upper()
            f.write(
                f"HETATM{i+1:5d} {atom.name:<4s} {ligcode:<3s} L 900    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{atom.occ:6.2f}{atom.b_iso:6.2f}"
                f"          {elem:>2s}\n"
            )
        f.write("END\n")

    rows.append({
        "pdbid": pdbid, "ligand": ligcode, "n_ca_matched": n,
        "ca_alignment_rmsd_A": round(float(align_rmsd), 3),
        "reference_pose_pdb": out_pdb,
    })
    print(f"{pdbid} ({ligcode}): {n} CA matched, alignment RMSD = {align_rmsd:.3f} A -> {out_pdb}")

for pdbid, reason in EXCLUDED.items():
    print(f"EXCLUDED {pdbid}: {reason}")

with open("results/lasb_ensemble_rmsd/holo_ligand_set_final.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"\n{len(rows)} ligands in final set. Written results/lasb_ensemble_rmsd/holo_ligand_set_final.csv")

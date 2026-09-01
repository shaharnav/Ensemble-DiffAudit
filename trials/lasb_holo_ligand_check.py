"""
Step 2b: check inclusion criteria for the 17 LasB structures with a bound
non-crystallization-additive ligand: catalytic-site binding (proximity to
Zn2+), occupancy, B-factor sanity.
"""
import gemmi
import numpy as np
import csv

STRUCTS = {
    "1U4G": "HPI", "3DBK": "RDF", "6F8B": "CXH", "6FZX": "EEK", "7OC7": "V85",
    "7QH1": "CI8", "7Z68": "IEV", "8R1B": "XI5", "8CC4": "U7F",
    "9FQD": "A1IEX", "9FQE": "A1IEU", "9FQX": "A1IJJ", "9FQY": "A1IFN",
    "9FRY": "A1IFU", "9FRZ": "A1IFZ", "9GM4": "A1IM0", "9FS0": "A1IFY",
}

rows = []
for pdbid, ligcode in STRUCTS.items():
    for ext in ["pdb", "cif"]:
        path = f"results/lasb_ensemble_rmsd/holo_structures/{pdbid}.{ext}"
        import os
        if os.path.exists(path):
            break
    st = gemmi.read_structure(path)
    st.setup_entities()
    model = st[0]

    zn_pos = None
    lig_atoms = []
    for chain in model:
        for res in chain:
            if res.name == "ZN":
                zn_pos = np.array([res[0].pos.x, res[0].pos.y, res[0].pos.z])
            if res.name == ligcode:
                for atom in res:
                    lig_atoms.append(atom)

    if zn_pos is None or not lig_atoms:
        rows.append({"pdbid": pdbid, "ligand": ligcode, "error": "zn_or_ligand_not_found"})
        continue

    lig_coords = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lig_atoms])
    min_dist_to_zn = np.min(np.linalg.norm(lig_coords - zn_pos, axis=1))
    occupancies = [a.occ for a in lig_atoms]
    bfactors = [a.b_iso for a in lig_atoms]

    rows.append({
        "pdbid": pdbid, "ligand": ligcode,
        "min_dist_to_zn_A": round(float(min_dist_to_zn), 2),
        "catalytic_site": min_dist_to_zn < 6.0,
        "mean_occupancy": round(float(np.mean(occupancies)), 2),
        "min_occupancy": round(float(np.min(occupancies)), 2),
        "mean_bfactor": round(float(np.mean(bfactors)), 1),
        "max_bfactor": round(float(np.max(bfactors)), 1),
        "n_atoms": len(lig_atoms),
    })

for r in rows:
    print(r)

with open("results/lasb_ensemble_rmsd/holo_ligand_site_check.csv", "w", newline="") as f:
    fieldnames = sorted(set().union(*[r.keys() for r in rows]))
    w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
    w.writeheader()
    w.writerows(rows)
print("\nWritten results/lasb_ensemble_rmsd/holo_ligand_site_check.csv")

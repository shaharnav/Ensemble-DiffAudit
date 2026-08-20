"""
Step 1: Zn2+/Ca2+ transplant into the 6 ConforMix/Boltz LasB conformers.

Conformers are pure-protein CIFs (zero HETATM) already aligned to the 1EZM
crystal frame. We do NOT re-superpose. Metal coordinates are copied directly
from 1EZM_apo.pdb into each conformer's existing frame, then validated:
distance to the coordinating residues (identified directly from 1EZM, not
assumed) and any steric clash (<1.5 A) to a protein atom.

Zn2+ coordination in 1EZM (from geometric search, <3.0 A): His140 NE2 (2.09 A),
His144 NE2 (2.05 A), Glu164 OE2 (1.84 A) -- the HExxH zinc-binding motif of
thermolysin-family (M4) metalloproteases.
Ca2+ coordination in 1EZM (<3.0 A): Glu172 OE1/OE2, Glu175 OE1/OE2, Asp183 OD1,
Asp136 OD2, Leu185 O.
"""
import gemmi
import numpy as np
import csv

ZN_COORD = np.array([49.202, 37.708, 19.451])
CA_COORD = np.array([51.525, 45.467, 31.069])

# Carboxylate O1/O2 are chemically symmetric (rotamer-arbitrary naming) --
# use whichever named oxygen is actually closest, not a fixed one.
ZN_LIGANDS = [(140, "HIS", ["NE2"]), (144, "HIS", ["NE2"]), (164, "GLU", ["OE1", "OE2"])]
CA_LIGANDS = [(172, "GLU", ["OE1", "OE2"]), (175, "GLU", ["OE1", "OE2"]),
              (183, "ASP", ["OD1", "OD2"]), (136, "ASP", ["OD1", "OD2"]),
              (185, "LEU", ["O"])]

CONFORMERS = ["beta0.0", "beta0.8", "beta1.6", "beta2.4", "beta3.2", "beta4.0"]

def load_chain(pdb_path):
    st = gemmi.read_structure(pdb_path)
    st.setup_entities()
    return st, st[0][0]

def get_atom_pos(chain, resnum, atomname):
    for res in chain:
        if res.seqid.num == resnum:
            for atom in res:
                if atom.name == atomname:
                    return np.array([atom.pos.x, atom.pos.y, atom.pos.z])
    return None

def all_protein_positions(chain):
    coords, labels = [], []
    for res in chain:
        for atom in res:
            coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
            labels.append(f"{res.name}{res.seqid.num}.{atom.name}")
    return np.array(coords), labels

rows = []
for beta in CONFORMERS:
    pdb_path = f"results/lasb_payload/ensemble_receptors_aligned/lasb_conformer_{beta}.pdb"
    st, chain = load_chain(pdb_path)

    def closest_dist(resnum, atomnames, target):
        best_name, best_d = None, None
        for an in atomnames:
            pos = get_atom_pos(chain, resnum, an)
            if pos is None:
                continue
            d = np.linalg.norm(pos - target)
            if best_d is None or d < best_d:
                best_d, best_name = d, an
        return best_name, best_d

    zn_dists = {}
    for resnum, resname, atomnames in ZN_LIGANDS:
        an, d = closest_dist(resnum, atomnames, ZN_COORD)
        zn_dists[f"{resname}{resnum}.{an}"] = d

    ca_dists = {}
    for resnum, resname, atomnames in CA_LIGANDS:
        an, d = closest_dist(resnum, atomnames, CA_COORD)
        ca_dists[f"{resname}{resnum}.{an}"] = d

    coords, labels = all_protein_positions(chain)
    zn_clash_d = np.min(np.linalg.norm(coords - ZN_COORD, axis=1))
    zn_clash_label = labels[np.argmin(np.linalg.norm(coords - ZN_COORD, axis=1))]
    ca_clash_d = np.min(np.linalg.norm(coords - CA_COORD, axis=1))
    ca_clash_label = labels[np.argmin(np.linalg.norm(coords - CA_COORD, axis=1))]

    zn_vals = [v for v in zn_dists.values() if v is not None]
    ca_vals = [v for v in ca_dists.values() if v is not None]
    # plausible coordination window: too short (<1.8 A) = physical overlap/clash,
    # too long (>3.2 A Zn / >3.5 A Ca) = coordination effectively broken
    zn_ok = all(1.8 <= v <= 3.2 for v in zn_vals)
    ca_ok = all(1.8 <= v <= 3.5 for v in ca_vals)
    # separate clash check: any protein atom outside the intended ligand set
    # coming within 1.5 A of the metal
    intended_zn_atoms = {f"{r}{n}" for n, r, a in ZN_LIGANDS}
    intended_ca_atoms = {f"{r}{n}" for n, r, a in CA_LIGANDS}
    clash = (zn_clash_d < 1.5 and zn_clash_label.split(".")[0] not in intended_zn_atoms) or \
            (ca_clash_d < 1.5 and ca_clash_label.split(".")[0] not in intended_ca_atoms)
    accept = zn_ok and ca_ok and not clash

    zn_max = max(zn_vals)
    ca_max = max(ca_vals)
    row = {"conformer": beta, "zn_max_ligand_dist": round(zn_max,2), "zn_closest_protein_atom": zn_clash_label,
           "zn_closest_protein_dist": round(zn_clash_d,2), "ca_max_ligand_dist": round(ca_max,2),
           "ca_closest_protein_atom": ca_clash_label, "ca_closest_protein_dist": round(ca_clash_d,2),
           "accept": accept}
    row.update({f"zn_{k}": round(v,2) for k,v in zn_dists.items()})
    row.update({f"ca_{k}": round(v,2) for k,v in ca_dists.items()})
    rows.append(row)

    print(f"\n=== {beta} ===")
    print("Zn2+ ligand distances:", {k: round(v,2) for k,v in zn_dists.items()})
    print("Ca2+ ligand distances:", {k: round(v,2) for k,v in ca_dists.items()})
    print(f"Zn2+ closest protein atom (any): {zn_clash_label} @ {zn_clash_d:.2f} A")
    print(f"Ca2+ closest protein atom (any): {ca_clash_label} @ {ca_clash_d:.2f} A")
    print(f"ACCEPT: {accept}")

    if accept:
        # write transplanted PDB: conformer protein atoms + Zn/Ca HETATM
        out_path = f"results/lasb_ensemble_rmsd/receptors_raw/{beta}_transplanted.pdb"
        import os
        os.makedirs("results/lasb_ensemble_rmsd/receptors_raw", exist_ok=True)
        doc = st.clone()
        with open(out_path, "w") as f:
            f.write(st.make_pdb_string())
            f.write(f"HETATM 9001  ZN   ZN A 302    {ZN_COORD[0]:8.3f}{ZN_COORD[1]:8.3f}{ZN_COORD[2]:8.3f}  1.00  0.00          ZN\n")
            f.write(f"HETATM 9002  CA   CA A 400    {CA_COORD[0]:8.3f}{CA_COORD[1]:8.3f}{CA_COORD[2]:8.3f}  1.00  0.00          CA\n")
            f.write("END\n")
        print(f"Written {out_path}")

fieldnames = sorted(set().union(*[r.keys() for r in rows]), key=lambda k: (k not in rows[0], k))
import os
os.makedirs("results/lasb_ensemble_rmsd", exist_ok=True)
with open("results/lasb_ensemble_rmsd/metal_transplant_validation.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
    w.writeheader()
    w.writerows(rows)
print("\nWritten results/lasb_ensemble_rmsd/metal_transplant_validation.csv")

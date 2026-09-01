"""
Step 1 (regenerated): validate the 6 Zn/Ca-explicit ConforMix/Boltz conformers.

Unlike the original transplant attempt, Zn2+/Ca2+ here are the model's OWN
predicted positions (chains B/C), not copied in from the crystal -- this is a
self-consistency check, not a transplant. Structures are aligned to the 1EZM
frame via CA superposition (needed for docking box placement; not needed for
the geometry check itself, which is internally self-consistent regardless of
frame).
"""
import gemmi
import numpy as np
import csv
import os

ZN_LIGANDS = [(140, "HIS", ["NE2"]), (144, "HIS", ["NE2"]), (164, "GLU", ["OE1", "OE2"])]
CA_LIGANDS = [(172, "GLU", ["OE1", "OE2"]), (175, "GLU", ["OE1", "OE2"]),
              (183, "ASP", ["OD1", "OD2"]), (136, "ASP", ["OD1", "OD2"]),
              (185, "LEU", ["O"])]
CONFORMERS = ["beta0.0", "beta0.8", "beta1.6", "beta2.4", "beta3.2", "beta4.0"]

REF_PATH = "results/lasb_payload/ensemble_receptors_aligned/1EZM_baseline_crystal.pdb"
ref_st = gemmi.read_structure(REF_PATH)
ref_st.setup_entities()
ref_chain = ref_st[0][0]
ref_ca = {res.seqid.num: res[0].pos for res in ref_chain
          if any(a.name == "CA" for a in res) for a in [next(a for a in res if a.name == "CA")]}

def kabsch(mob, ref):
    mob_c, ref_c = mob.mean(axis=0), ref.mean(axis=0)
    H = (mob - mob_c).T @ (ref - ref_c)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = ref_c - R @ mob_c
    return R, t

def closest_dist(chain, resnum, atomnames, target):
    best_name, best_d = None, None
    for res in chain:
        if res.seqid.num != resnum:
            continue
        for an in atomnames:
            for atom in res:
                if atom.name == an:
                    d = np.linalg.norm(np.array([atom.pos.x, atom.pos.y, atom.pos.z]) - target)
                    if best_d is None or d < best_d:
                        best_d, best_name = d, an
    return best_name, best_d

OUT_DIR = "results/lasb_ensemble_rmsd/receptors_raw"
os.makedirs(OUT_DIR, exist_ok=True)

rows = []
for beta in CONFORMERS:
    cif_path = f"results/lasb_ensemble_rmsd/conformers_zn_regenerated/{beta}.cif"
    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    model = st[0]
    protein_chain = model["A"]
    zn_pos = np.array([model["B"][0][0].pos.x, model["B"][0][0].pos.y, model["B"][0][0].pos.z])
    ca_pos = np.array([model["C"][0][0].pos.x, model["C"][0][0].pos.y, model["C"][0][0].pos.z])

    zn_dists = {}
    for resnum, resname, atomnames in ZN_LIGANDS:
        an, d = closest_dist(protein_chain, resnum, atomnames, zn_pos)
        zn_dists[f"{resname}{resnum}.{an}"] = d
    ca_dists = {}
    for resnum, resname, atomnames in CA_LIGANDS:
        an, d = closest_dist(protein_chain, resnum, atomnames, ca_pos)
        ca_dists[f"{resname}{resnum}.{an}"] = d

    zn_vals = list(zn_dists.values())
    ca_vals = list(ca_dists.values())
    # Same acceptance window as the original transplant validation:
    # <1.8 A = physical overlap/clash, >3.2 A (Zn) / >3.5 A (Ca) = broken coordination
    zn_ok = all(1.8 <= v <= 3.2 for v in zn_vals)
    ca_ok = all(1.8 <= v <= 3.5 for v in ca_vals)
    accept = zn_ok and ca_ok

    row = {"conformer": beta, "accept": accept,
           "zn_max": round(max(zn_vals), 2), "zn_min": round(min(zn_vals), 2),
           "ca_max": round(max(ca_vals), 2), "ca_min": round(min(ca_vals), 2)}
    row.update({f"zn_{k}": round(v, 2) for k, v in zn_dists.items()})
    row.update({f"ca_{k}": round(v, 2) for k, v in ca_dists.items()})
    rows.append(row)
    print(f"{beta}: Zn {[round(v,2) for v in zn_vals]}  Ca {[round(v,2) for v in ca_vals]}  ACCEPT={accept}")

    # Align to 1EZM frame and write a receptor PDB (protein + this model's own Zn/Ca)
    mob_ca, ref_ca_matched = [], []
    for res in protein_chain:
        atom = next((a for a in res if a.name == "CA"), None)
        if atom is not None and res.seqid.num in ref_ca:
            mob_ca.append([atom.pos.x, atom.pos.y, atom.pos.z])
            p = ref_ca[res.seqid.num]
            ref_ca_matched.append([p.x, p.y, p.z])
    R, t = kabsch(np.array(mob_ca), np.array(ref_ca_matched))

    out_pdb = f"{OUT_DIR}/{beta}_zn_regenerated.pdb"
    with open(out_pdb, "w") as f:
        serial = 1
        for res in protein_chain:
            for atom in res:
                xyz = R @ np.array([atom.pos.x, atom.pos.y, atom.pos.z]) + t
                f.write(f"ATOM  {serial:5d} {atom.name:<4s} {res.name:<3s} A{res.seqid.num:4d}    "
                        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{atom.occ:6.2f}{atom.b_iso:6.2f}"
                        f"          {atom.element.name:>2s}\n")
                serial += 1
        for label, pos, resnum in [("ZN", zn_pos, 900), ("CA", ca_pos, 901)]:
            xyz = R @ pos + t
            f.write(f"HETATM{serial:5d} {label:<4s} {label:<3s} A{resnum:4d}    "
                    f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {label:>2s}\n")
            serial += 1
        f.write("END\n")
    print(f"  -> {out_pdb}")

fieldnames = sorted(set().union(*[r.keys() for r in rows]), key=lambda k: (k not in rows[0], k))
with open("results/lasb_ensemble_rmsd/zn_regenerated_validation.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
    w.writeheader()
    w.writerows(rows)

n_accept = sum(1 for r in rows if r["accept"])
print(f"\n{n_accept}/6 conformers ACCEPT")
print("Written results/lasb_ensemble_rmsd/zn_regenerated_validation.csv")

"""
Step 3 prep: prepare each ligand's own source holo structure as a rigid
receptor (Condition C, the self-docking ceiling) through the identical meeko
pipeline used for the apo crystal and ensemble conformers.

Strips the ligand + waters, keeps only protein + Zn2+/Ca2+ (and any other
structural metals present), matching the receptor content used elsewhere in
this experiment.
"""
import csv, os, sys
import gemmi
import numpy as np
sys.path.insert(0, ".")
from casf_pipeline.prep import prepare_receptor_pdbqt

METALS_KEEP = {"ZN", "CA", "MG", "MN", "FE", "CU", "NA", "K"}

with open("results/lasb_ensemble_rmsd/holo_ligand_set_final.csv") as f:
    rows = list(csv.DictReader(f))

OUT_DIR = "results/lasb_ensemble_rmsd/receptors_raw"
os.makedirs(OUT_DIR, exist_ok=True)

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

results = []
for r in rows:
    pdbid, ligcode = r["pdbid"], r["ligand"]
    for ext in ["pdb", "cif"]:
        path = f"results/lasb_ensemble_rmsd/holo_structures/{pdbid}.{ext}"
        if os.path.exists(path):
            break
    st = gemmi.read_structure(path)
    st.setup_entities()
    model = st[0]

    # Same chain-selection logic as ligand prep: the chain carrying this ligand
    lig_chain = None
    for chain in model:
        for res in chain:
            if res.name == ligcode:
                lig_chain = chain
                break
        if lig_chain is not None:
            break

    # Align this receptor's chain onto the 1EZM frame -- the docking box is
    # defined in that frame, and each holo structure has its own arbitrary
    # crystallographic origin, so an unaligned receptor here would put the
    # box nowhere near the actual binding site.
    mob_ca, ref_ca_matched = [], []
    for res in lig_chain:
        atom = next((a for a in res if a.name == "CA"), None)
        if atom is not None and res.seqid.num in ref_ca:
            mob_ca.append([atom.pos.x, atom.pos.y, atom.pos.z])
            p = ref_ca[res.seqid.num]
            ref_ca_matched.append([p.x, p.y, p.z])
    R, t = kabsch(np.array(mob_ca), np.array(ref_ca_matched))

    out_pdb = f"{OUT_DIR}/{pdbid}_holo_receptor.pdb"
    with open(out_pdb, "w") as f:
        serial = 1
        for res in lig_chain:
            is_protein = res.het_flag == "A" or gemmi.find_tabulated_residue(res.name) and gemmi.find_tabulated_residue(res.name).is_amino_acid()
            is_metal = res.name in METALS_KEEP
            if not (is_protein or is_metal):
                continue
            record = "ATOM  " if is_protein else "HETATM"
            for atom in res:
                if atom.altloc not in ("\x00", "A"):
                    continue
                elem = atom.element.name
                xyz = R @ np.array([atom.pos.x, atom.pos.y, atom.pos.z]) + t
                f.write(
                    f"{record}{serial:5d} {atom.name:<4s} {res.name:<3s} {lig_chain.name}{res.seqid.num:4d}    "
                    f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{atom.occ:6.2f}{atom.b_iso:6.2f}"
                    f"          {elem:>2s}\n"
                )
                serial += 1
        f.write("END\n")

    out_pdbqt = f"{OUT_DIR}/{pdbid}_holo_receptor.pdbqt"
    ok = prepare_receptor_pdbqt(out_pdb, out_pdbqt)
    results.append({"pdbid": pdbid, "ligand": ligcode, "receptor_pdb": out_pdb,
                     "receptor_pdbqt": out_pdbqt if ok else "", "prep_ok": ok})
    print(f"{pdbid}: prep {'ok' if ok else 'FAILED'}")

n_ok = sum(1 for r in results if r["prep_ok"])
print(f"\n{n_ok}/{len(results)} holo receptors prepped")

with open("results/lasb_ensemble_rmsd/holo_receptor_prep_log.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader()
    w.writerows(results)
print("Written results/lasb_ensemble_rmsd/holo_receptor_prep_log.csv")

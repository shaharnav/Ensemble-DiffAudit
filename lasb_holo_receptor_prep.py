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
sys.path.insert(0, ".")
from casf_pipeline.prep import prepare_receptor_pdbqt

METALS_KEEP = {"ZN", "CA", "MG", "MN", "FE", "CU", "NA", "K"}

with open("results/lasb_ensemble_rmsd/holo_ligand_set_final.csv") as f:
    rows = list(csv.DictReader(f))

OUT_DIR = "results/lasb_ensemble_rmsd/receptors_raw"
os.makedirs(OUT_DIR, exist_ok=True)

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
                f.write(
                    f"{record}{serial:5d} {atom.name:<4s} {res.name:<3s} {lig_chain.name}{res.seqid.num:4d}    "
                    f"{atom.pos.x:8.3f}{atom.pos.y:8.3f}{atom.pos.z:8.3f}{atom.occ:6.2f}{atom.b_iso:6.2f}"
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

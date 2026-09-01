"""
Step 3 prep: convert each aligned holo-ligand reference pose (already in the
1EZM frame) into an SDF (correct bond orders via CCD template) and a pdbqt
(crystal coordinates preserved, no re-embedding).

Builds the RDKit mol directly from gemmi atoms/coordinates rather than
round-tripping through fixed-column PDB text -- several LasB ligand codes are
5-character extended CCD IDs (e.g. A1IEX) that don't fit the 3-column PDB
resName field, which silently corrupted the coordinate/element columns on
reparse.
"""
import csv, os, sys, warnings
import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Geometry import Point3D
sys.path.insert(0, ".")
from casf_pipeline.extract import get_ccd_smiles
from casf_pipeline.prep import prepare_ligand_pdbqt_from_sdf, ligand_bbox_center_size

STRUCTS = {
    "1U4G": "HPI", "3DBK": "RDF", "6F8B": "CXH", "6FZX": "EEK",
    "7QH1": "CI8", "8R1B": "XI5", "8CC4": "U7F",
    "9FQD": "A1IEX", "9FQE": "A1IEU", "9FQX": "A1IJJ", "9FQY": "A1IFN",
    "9FRY": "A1IFU", "9FRZ": "A1IFZ", "9GM4": "A1IM0", "9FS0": "A1IFY",
}

REF_PATH = "results/lasb_payload/ensemble_receptors_aligned/1EZM_baseline_crystal.pdb"
ref_st = gemmi.read_structure(REF_PATH)
ref_st.setup_entities()
ref_chain = ref_st[0][0]
ref_ca = {res.seqid.num: res[0].pos for res in ref_chain
          if any(a.name == "CA" for a in res) for a in [next(a for a in res if a.name == "CA")]}

def kabsch_align(mob_ca, ref_ca_arr):
    mob_c = mob_ca.mean(axis=0)
    ref_c = ref_ca_arr.mean(axis=0)
    H = (mob_ca - mob_c).T @ (ref_ca_arr - ref_c)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = ref_c - R @ mob_c
    return R, t

def build_mol_from_atoms(ligcode, atoms, coords):
    smiles = get_ccd_smiles(ligcode)
    if not smiles:
        return None, "no_ccd_smiles"
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        return None, "bad_template"

    # Write a fixed-column PDB block using a safe 3-char placeholder resName
    # ("LIG") -- extended 5-char CCD codes (e.g. A1IEX) don't fit the PDB
    # resName field and corrupt fixed-column parsing if used directly. The
    # real ligcode is passed separately to get_ccd_smiles above, not re-read
    # from this file. RDKit's PDB parser infers connectivity from covalent
    # radii, which is more robust for phosphonate/macrocycle chemistry than
    # a simple distance-cutoff connectivity pass.
    heavy = [(atom, xyz) for atom, xyz in zip(atoms, coords) if atom.element.name != "H"]
    lines = []
    for i, (atom, xyz) in enumerate(heavy):
        elem = atom.element.name
        name = atom.name if len(atom.name) <= 4 else atom.name[:4]
        lines.append(
            f"HETATM{i+1:5d} {name:<4s} LIG L 900    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{atom.occ:6.2f}{atom.b_iso:6.2f}"
            f"          {elem:>2s}"
        )
    pdb_block = "\n".join(lines) + "\nEND\n"
    raw_mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=True)
    if raw_mol is None:
        return None, "molfrompdbblock_failed"

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = AllChem.AssignBondOrdersFromTemplate(template, raw_mol)
        Chem.SanitizeMol(mol)
    except Exception as e:
        return None, f"assign_bond_orders_failed: {e}"

    mol.SetProp("_Name", ligcode)
    return mol, None

OUT_SDF_DIR = "results/lasb_ensemble_rmsd/ligands_sdf"
OUT_PDBQT_DIR = "results/lasb_ensemble_rmsd/ligands_pdbqt"
os.makedirs(OUT_SDF_DIR, exist_ok=True)
os.makedirs(OUT_PDBQT_DIR, exist_ok=True)

results = []
for pdbid, ligcode in STRUCTS.items():
    for ext in ["pdb", "cif"]:
        path = f"results/lasb_ensemble_rmsd/holo_structures/{pdbid}.{ext}"
        if os.path.exists(path):
            break
    st = gemmi.read_structure(path)
    st.setup_entities()
    model = st[0]

    lig_chain = None
    for chain in model:
        for res in chain:
            if res.name == ligcode:
                lig_chain = chain
                break
        if lig_chain is not None:
            break

    # Some structures (e.g. CXH/6F8B) have two separate copies of the ligand
    # residue bound to the same protein chain (different resSeq, not an
    # altloc split) -- restrict to the single copy nearest the catalytic Zn.
    lig_residues = [res for res in lig_chain if res.name == ligcode]
    if len(lig_residues) > 1:
        zn_pos = None
        for res in lig_chain:
            if res.name == "ZN":
                zn_pos = np.array([res[0].pos.x, res[0].pos.y, res[0].pos.z])
        def min_dist_to_zn(res):
            return min(np.linalg.norm(np.array([a.pos.x, a.pos.y, a.pos.z]) - zn_pos) for a in res)
        lig_residues = [min(lig_residues, key=min_dist_to_zn)]
    lig_atoms = [a for res in lig_residues for a in res]
    mob_ca, ref_ca_matched = [], []
    for res in lig_chain:
        atom = next((a for a in res if a.name == "CA"), None)
        if atom is not None and res.seqid.num in ref_ca:
            mob_ca.append([atom.pos.x, atom.pos.y, atom.pos.z])
            p = ref_ca[res.seqid.num]
            ref_ca_matched.append([p.x, p.y, p.z])
    mob_ca = np.array(mob_ca)
    ref_ca_matched = np.array(ref_ca_matched)
    R, t = kabsch_align(mob_ca, ref_ca_matched)

    lig_coords = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lig_atoms])
    lig_aligned = (R @ lig_coords.T).T + t

    mol, err = build_mol_from_atoms(ligcode, lig_atoms, lig_aligned)
    row = {"pdbid": pdbid, "ligand": ligcode}
    if mol is None:
        row.update({"sdf_ok": False, "pdbqt_ok": False, "error": err})
        results.append(row)
        print(f"{pdbid} ({ligcode}): FAILED -- {err}")
        continue

    sdf_path = f"{OUT_SDF_DIR}/{pdbid}_{ligcode}.sdf"
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()

    pdbqt_path = f"{OUT_PDBQT_DIR}/{pdbid}_{ligcode}.pdbqt"
    pdbqt_ok = prepare_ligand_pdbqt_from_sdf(sdf_path, pdbqt_path)
    center, size = ligand_bbox_center_size(sdf_path)

    row.update({
        "hac": mol.GetNumHeavyAtoms(), "sdf_ok": True, "pdbqt_ok": pdbqt_ok,
        "sdf_path": sdf_path, "pdbqt_path": pdbqt_path if pdbqt_ok else "",
        "bbox_center": center, "bbox_size": size,
    })
    results.append(row)
    print(f"{pdbid} ({ligcode}): hac={mol.GetNumHeavyAtoms()}, sdf ok, pdbqt {'ok' if pdbqt_ok else 'FAILED'}")

n_ok = sum(1 for r in results if r.get("pdbqt_ok"))
print(f"\n{n_ok}/{len(results)} ligands fully prepped")

with open("results/lasb_ensemble_rmsd/ligand_prep_log.csv", "w", newline="") as f:
    fieldnames = sorted(set().union(*[r.keys() for r in results]))
    w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
    w.writeheader()
    w.writerows(results)
print("Written results/lasb_ensemble_rmsd/ligand_prep_log.csv")

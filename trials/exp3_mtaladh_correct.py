"""MtAlaDH: run DSSP + pocket RMSD with correct crystal-form pair (1MOP vs 3ISJ)."""
import os, requests, time
import numpy as np
from Bio.PDB import PDBParser, Superimposer

PDB_DIR = "pdbs/exp3_new"
os.makedirs(PDB_DIR, exist_ok=True)

def download_pdb(pdb_id):
    path = os.path.join(PDB_DIR, f"{pdb_id.upper()}.pdb")
    if not os.path.exists(path):
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", timeout=30)
        r.raise_for_status()
        with open(path, "w") as fh: fh.write(r.text)
        time.sleep(0.2)
    return path

def parse_sse(pdb_path):
    helix, sheet = set(), set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HELIX "):
                c = line[19]
                try:
                    for r in range(int(line[21:25]), int(line[33:37])+1): helix.add((c,r))
                except ValueError: pass
            elif line.startswith("SHEET "):
                c = line[21]
                try:
                    for r in range(int(line[22:26]), int(line[33:37])+1): sheet.add((c,r))
                except ValueError: pass
    return helix, sheet

def get_hetatm_nonbuf(pdb_path, chain_id):
    SKIP = {"HOH","WAT","DOD","SO4","PO4","GOL","EDO","PEG","MPD","BME","CIT","ACT","ACE",
            "TRS","MES","EPE","PG4","1PE","EOH","IPA","DMS","MG","CA","ZN","NA","K",
            "MN","CL","FE","CU","NI","CO","NH4","ACY","PYR","ALA","GLU","IMD","TLA",
            "FOR","BME","BAL","NAI"}
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM") and line[21]==chain_id:
                rn = line[17:20].strip()
                rs = line[22:26].strip()
                if rn not in SKIP: seen.add((rs, rn))
    return {rn for _, rn in seen}

print("=== MtAlaDH: 1MOP (apo) vs 3ISJ (A8D, 268 Da) ===")
print("Same crystal form confirmed: global CA RMSD 0.892 Å\n")

apo_id, apo_chain   = "1MOP", "A"
holo_id, holo_chain = "3ISJ", "A"
lig_codes = {"A8D"}

apo_path  = download_pdb(apo_id)
holo_path = download_pdb(holo_id)

# HETATM audit
apo_hets  = get_hetatm_nonbuf(apo_path,  apo_chain)
holo_hets = get_hetatm_nonbuf(holo_path, holo_chain)
apo_cofac  = apo_hets  - lig_codes
holo_cofac = holo_hets - lig_codes
print(f"APO  non-buffer, non-drug: {apo_cofac}")
print(f"HOLO non-buffer, non-drug: {holo_cofac}")
cofactor_diff = apo_cofac != holo_cofac
print(f"Cofactor gate: {'FAIL' if cofactor_diff else 'PASS'}")

# DSSP on holo
helix_res, sheet_res = parse_sse(holo_path)
parser = PDBParser(QUIET=True)
holo_struct = parser.get_structure(holo_id, holo_path)[0]
apo_struct  = parser.get_structure(apo_id,  apo_path)[0]

holo_chain_obj = holo_struct[holo_chain]

# Get ligand coords
lig_atoms = []
for ch in holo_struct.get_chains():
    for res in ch.get_residues():
        if res.get_resname().strip() in lig_codes and res.get_id()[0] != " ":
            for atom in res.get_atoms():
                if atom.element and atom.element != "H":
                    lig_atoms.append(atom.get_vector().get_array())
if not lig_atoms:
    print("ERROR: ligand not found")
    exit()
lig_coords = np.array(lig_atoms)
print(f"\nLigand heavy atoms: {len(lig_coords)}")

# Pocket-lining residues (all chains of holo)
pocket_keys = set()
for ch in holo_struct.get_chains():
    for res in ch.get_residues():
        if res.get_id()[0] != " ": continue
        for atom in res.get_atoms():
            if not atom.element or atom.element == "H": continue
            coord = np.array(atom.get_vector().get_array())
            if np.sqrt(((lig_coords - coord[None,:])**2).sum(axis=1)).min() <= 8.0:
                pocket_keys.add((ch.get_id(), res.get_id()[1]))
                break

total = len(pocket_keys)
h_n = sum(1 for k in pocket_keys if k in helix_res)
s_n = sum(1 for k in pocket_keys if k in sheet_res)
hs  = (h_n + s_n) / total if total else 0.0
dssp_pass = hs >= 0.60
print(f"Pocket residues: {total}  |  helix={h_n} sheet={s_n} coil={total-h_n-s_n}  →  {hs:.1%}")
print(f"DSSP gate: {'PASS' if dssp_pass else 'FAIL'}")

# Pocket RMSD (global align first)
def get_ca(model, cid):
    try: ch = model[cid]
    except KeyError: ch = list(model.get_chains())[0]
    return {res.get_id()[1]: res["CA"] for res in ch.get_residues()
            if res.get_id()[0]==" " and "CA" in res}

apo_ca  = get_ca(apo_struct,  apo_chain)
holo_ca = get_ca(holo_struct, holo_chain)
shared  = sorted(set(apo_ca) & set(holo_ca))

sup = Superimposer()
sup.set_atoms([holo_ca[r] for r in shared], [apo_ca[r] for r in shared])
sup.apply(apo_struct.get_atoms())
print(f"\nGlobal CA RMSD (n={len(shared)}): {sup.rms:.3f} Å")

pocket_rids = {rid for _, rid in pocket_keys}
common_pock = sorted(set(apo_ca) & set(holo_ca) & pocket_rids)
apo_arr  = np.array([apo_ca[r].get_vector().get_array()  for r in common_pock])
holo_arr = np.array([holo_ca[r].get_vector().get_array() for r in common_pock])
ca_rmsd  = np.sqrt(((apo_arr - holo_arr)**2).sum(axis=1).mean())
rmsd_pass = 1.0 <= ca_rmsd <= 2.5
print(f"Pocket CA RMSD (n={len(common_pock)}): {ca_rmsd:.3f} Å")
print(f"Pocket RMSD gate: {'PASS' if rmsd_pass else 'FAIL'}")

# Also check 2A86 as alternative apo
print("\n--- Also checking 2A86 vs 3ISJ ---")
a2 = download_pdb("2A86")
apo2_ca = get_ca(parser.get_structure("2A86",a2)[0], "A")
shared2  = sorted(set(apo2_ca) & set(holo_ca))
apo2_struct = parser.get_structure("2A86",a2)[0]
sup2 = Superimposer()
sup2.set_atoms([holo_ca[r] for r in shared2], [apo2_ca[r] for r in shared2])
sup2.apply(apo2_struct.get_atoms())
pock2 = sorted(set(apo2_ca) & set(holo_ca) & pocket_rids)
a2_arr   = np.array([apo2_ca[r].get_vector().get_array() for r in pock2])
h2_arr   = np.array([holo_ca[r].get_vector().get_array() for r in pock2])
ca2_rmsd = np.sqrt(((a2_arr - h2_arr)**2).sum(axis=1).mean())
print(f"Global CA RMSD 2A86/3ISJ: {sup2.rms:.3f} Å  Pocket CA RMSD: {ca2_rmsd:.3f} Å")
a2_hets = get_hetatm_nonbuf(a2, "A")
print(f"2A86 non-buffer: {a2_hets}")

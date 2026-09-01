"""DSSP gate for 4 PocketMiner pairs. Downloads holo PDB, measures helix/sheet fraction
over pocket-lining residues (8 Å from holo ligand heavy atoms, any-atom distance).
Uses HELIX/SHEET records from PDB files."""
import os, requests, time, re
import numpy as np
from Bio.PDB import PDBParser

PDB_DIR = "pdbs/exp3_pocketminer"
os.makedirs(PDB_DIR, exist_ok=True)

BUFFER_CODES = {
    "SO4","PO4","GOL","EDO","PEG","MPD","BME","DTT","CIT","ACT","ACE",
    "TRS","MES","HEP","TAR","TLA","FMT","IMD","EPE","PGE","PG4","1PE",
    "15P","EOH","ETF","IOD","CL","NA","K","MG","CA","ZN","MN","FE","CU",
    "NI","CO","CD","HG","PT","AU","BR","I","F","IPA","DMS","GLC","FRU",
    "GAL","MAN","XYL","RIB","FUC","NAG","NDG","BMA","MMA","SEP","PTR",
    "TPO","CSO","HOH","WAT","DOD","ACY","NH4","PYR","ALA","GLU"
}

def download_pdb(pdb_id, pdb_dir=PDB_DIR):
    path = os.path.join(pdb_dir, f"{pdb_id.upper()}.pdb")
    if not os.path.exists(path):
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        with open(path, "w") as f:
            f.write(r.text)
        time.sleep(0.2)
    return path

def parse_sse(pdb_path):
    """Return sets of (chain, res_id) for helix and sheet residues."""
    helix_res = set()
    sheet_res = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HELIX "):
                c = line[19]
                try:
                    start, end = int(line[21:25].strip()), int(line[33:37].strip())
                    for r in range(start, end + 1):
                        helix_res.add((c, r))
                except ValueError:
                    pass
            elif line.startswith("SHEET "):
                c = line[21]
                try:
                    start, end = int(line[22:26].strip()), int(line[33:37].strip())
                    for r in range(start, end + 1):
                        sheet_res.add((c, r))
                except ValueError:
                    pass
    return helix_res, sheet_res

def dssp_gate(apo_pdb_id, holo_pdb_id, holo_chain, ligand_codes,
              label, motion_type, pocket_ca_rmsd, ca_aa_ratio):
    print(f"\n{'─'*60}")
    print(f"{label}  ({apo_pdb_id}/{holo_pdb_id}, chain {holo_chain})")
    print(f"  Motion: {motion_type}  |  Existing pocket CA RMSD: {pocket_ca_rmsd:.3f} Å  |  CA/AA: {ca_aa_ratio:.3f}")

    holo_path = download_pdb(holo_pdb_id)
    helix_res, sheet_res = parse_sse(holo_path)

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(holo_pdb_id, holo_path)[0]

    try:
        chain = struct[holo_chain]
    except KeyError:
        # Try first chain
        chain = list(struct.get_chains())[0]
        print(f"  WARNING: chain {holo_chain} not found, using chain {chain.id}")
        holo_chain = chain.id

    # Get ligand heavy atoms
    lig_atoms = []
    for res in chain.get_residues():
        if res.get_resname().strip() in ligand_codes and res.get_id()[0] != " ":
            for atom in res.get_atoms():
                if atom.element and atom.element != "H":
                    lig_atoms.append(atom.get_vector().get_array())
    if not lig_atoms:
        # Try all chains
        for ch in struct.get_chains():
            for res in ch.get_residues():
                if res.get_resname().strip() in ligand_codes and res.get_id()[0] != " ":
                    for atom in res.get_atoms():
                        if atom.element and atom.element != "H":
                            lig_atoms.append(atom.get_vector().get_array())
            if lig_atoms:
                break
    if not lig_atoms:
        print(f"  ERROR: ligand {ligand_codes} not found")
        return None

    lig_coords = np.array(lig_atoms)
    print(f"  Ligand heavy atoms: {len(lig_coords)}")

    # Pocket-lining residues: protein heavy atoms within 8 Å
    pocket_res = set()
    for ch in struct.get_chains():
        for res in ch.get_residues():
            if res.get_id()[0] != " ":  # hetflag — skip
                continue
            rid = res.get_id()[1]
            cid = ch.get_id()
            for atom in res.get_atoms():
                if atom.element == "H" or not atom.element:
                    continue
                coord = np.array(atom.get_vector().get_array())
                diffs = lig_coords - coord[None, :]
                if np.sqrt((diffs**2).sum(axis=1)).min() <= 8.0:
                    pocket_res.add((cid, rid))
                    break

    print(f"  Pocket-lining residues (≤8 Å): {len(pocket_res)}")

    helix_n = sum(1 for key in pocket_res if key in helix_res)
    sheet_n = sum(1 for key in pocket_res if key in sheet_res)
    coil_n  = len(pocket_res) - helix_n - sheet_n
    total   = len(pocket_res)
    hs_frac = (helix_n + sheet_n) / total if total else 0.0

    print(f"  SSE: helix={helix_n}, sheet={sheet_n}, coil={coil_n}  →  {hs_frac:.1%}")
    gate = "PASS" if hs_frac >= 0.60 else "FAIL"
    print(f"  DSSP gate (≥60%): {gate}")

    rmsd_gate = "PASS" if 1.0 <= pocket_ca_rmsd <= 2.5 else "FAIL"
    print(f"  Pocket RMSD gate (1.0–2.5 Å): {rmsd_gate}")

    return {
        "label": label,
        "apo": apo_pdb_id,
        "holo": holo_pdb_id,
        "motion_type": motion_type,
        "pocket_ca_rmsd": pocket_ca_rmsd,
        "ca_aa_ratio": ca_aa_ratio,
        "n_pocket_res": total,
        "helix_n": helix_n,
        "sheet_n": sheet_n,
        "coil_n": coil_n,
        "hs_frac": hs_frac,
        "dssp_pass": hs_frac >= 0.60,
        "rmsd_pass": 1.0 <= pocket_ca_rmsd <= 2.5,
    }


# PocketMiner pairs — apo/holo from Exp2 screen CSV, chains from CSV
pairs = [
    # (apo, holo, chain, [lig_codes], label, motion, pocket_ca_rmsd, ca_aa_ratio)
    ("5NIA", "5NI6", "A", ["DJ3"],
     "LTA4H (Leukotriene A-4 hydrolase) — H. sapiens", "ssm, loop", 2.336, 0.793),
    ("3NX1", "3NX2", "A", ["FER"],
     "FDC (Ferulic acid decarboxylase) — Enterobacter sp.", "ssm", 1.937, 0.764),
    ("2OHG", "2OHV", "A", ["NHL"],
     "GluR (Glutamate racemase) — S. pyogenes", "loop", 1.662, 0.907),
    ("1EZM", "3DBK", "A", ["RDF"],
     "LasB (Elastase) — P. aeruginosa", "id", 1.569, 0.895),
]

results = []
for args in pairs:
    r = dssp_gate(*args)
    if r:
        results.append(r)

print("\n\n" + "="*70)
print("SUMMARY TABLE — PocketMiner pairs")
print("="*70)
print(f"{'Pair':<12} {'Protein':<32} {'Motion':<12} {'CA RMSD':>8} {'DSSP%':>6} {'DSSP':>5} {'RMSD':>5}")
print("-"*85)
for r in results:
    pair = f"{r['apo']}/{r['holo']}"
    prot = r["label"].split("(")[1].split(")")[0] if "(" in r["label"] else r["label"][:30]
    flag = " *LOOP*" if "loop" in r["motion_type"] else ""
    print(f"{pair:<12} {prot:<32} {r['motion_type']:<12} {r['pocket_ca_rmsd']:>8.3f} "
          f"{r['hs_frac']:>6.1%} {'✓' if r['dssp_pass'] else '✗':>5} "
          f"{'✓' if r['rmsd_pass'] else '✗':>5}{flag}")

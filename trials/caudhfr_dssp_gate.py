"""DSSP gate for CauDHFR 8A0Z using biopython + PDB HELIX/SHEET records."""
import numpy as np
from Bio.PDB import PDBParser

pdb_path = "pdbs/exp3_caudhfr/8A0Z.pdb"

# Parse HELIX and SHEET records directly from PDB file
helix_residues = set()  # (chain, res_id)
sheet_residues = set()

with open(pdb_path) as f:
    for line in f:
        if line.startswith("HELIX "):
            chain = line[19]
            start = int(line[21:25].strip())
            end = int(line[33:37].strip())
            for r in range(start, end + 1):
                helix_residues.add((chain, r))
        elif line.startswith("SHEET "):
            chain = line[21]
            start = int(line[22:26].strip())
            end = int(line[33:37].strip())
            for r in range(start, end + 1):
                sheet_residues.add((chain, r))

print(f"HELIX residue positions (chain, res_id): {len(helix_residues)}")
print(f"SHEET residue positions (chain, res_id): {len(sheet_residues)}")

# Parse structure
parser = PDBParser(QUIET=True)
struct = parser.get_structure("8A0Z", pdb_path)
model = struct[0]

# Use chain B
chain = model["B"]

# Get CP6 heavy atoms
cp6_atoms = [a for a in chain.get_atoms()
             if a.get_parent().get_resname() == "CP6" and a.element != "H"]
print(f"\nCP6 heavy atoms in chain B: {len(cp6_atoms)}")

if not cp6_atoms:
    chain = model["A"]
    cp6_atoms = [a for a in chain.get_atoms()
                 if a.get_parent().get_resname() == "CP6" and a.element != "H"]
    print(f"CP6 heavy atoms in chain A: {len(cp6_atoms)}")

cp6_coords = np.array([a.get_vector().get_array() for a in cp6_atoms])

# Pocket-lining residues: any protein heavy atom within 8 Å of any CP6 heavy atom
pocket_res = set()
for res in chain.get_residues():
    if res.get_resname() == "CP6":
        continue
    res_id = res.get_id()[1]
    res_chain = res.get_parent().get_id()
    # Check if it's a standard amino acid residue
    if res.get_id()[0] != " ":  # hetflag
        continue
    for atom in res.get_atoms():
        if atom.element == "H":
            continue
        coord = np.array(atom.get_vector().get_array())
        diffs = cp6_coords - coord[None, :]
        dists = np.sqrt((diffs**2).sum(axis=1))
        if dists.min() <= 8.0:
            pocket_res.add((res_chain, res_id, res.get_resname()))
            break

print(f"\nPocket-lining residues (any heavy atom ≤8 Å from any CP6): {len(pocket_res)}")

# Classify each pocket residue
helix_count = 0
sheet_count = 0
coil_count = 0
print("\nPocket residues and SSE:")
for (c, rid, rname) in sorted(pocket_res, key=lambda x: x[1]):
    key = (c, rid)
    if key in helix_residues:
        sse = "helix"
        helix_count += 1
    elif key in sheet_residues:
        sse = "sheet"
        sheet_count += 1
    else:
        sse = "coil"
        coil_count += 1
    print(f"  {c} {rid:4d} {rname}: {sse}")

total = len(pocket_res)
hs = helix_count + sheet_count
fraction = hs / total if total > 0 else 0.0
print(f"\nHelix: {helix_count}, Sheet: {sheet_count}, Coil: {coil_count}")
print(f"Helix+sheet fraction: {hs}/{total} = {fraction:.1%}")
print(f"DSSP gate (floor ≥60%): {'PASS' if fraction >= 0.60 else 'FAIL'}")

"""Pocket CA RMSD gate: 8A0N vs 8A0Z for CauDHFR, with global alignment first."""
import numpy as np
from Bio.PDB import PDBParser, Superimposer

parser = PDBParser(QUIET=True)
s_apo = parser.get_structure("8A0N", "pdbs/exp3_caudhfr/8A0N.pdb")[0]
s_holo = parser.get_structure("8A0Z", "pdbs/exp3_caudhfr/8A0Z.pdb")[0]

# Pocket-lining residue IDs (chain B, from DSSP gate)
pocket_res_ids = {
    8, 9, 10, 11, 12, 19, 23, 24, 25, 27, 29, 30, 31, 32, 33, 34, 35, 36, 37,
    39, 40, 52, 53, 54, 56, 57, 58, 59, 60, 61, 62, 68, 118, 119, 120, 121,
    122, 126, 139, 140, 141, 142, 143, 196
}

def get_chain_ca(model, chain_id):
    """All CA atoms from protein residues in chain, sorted by res_id."""
    chain = model[chain_id]
    atoms = []
    for res in sorted(chain.get_residues(), key=lambda r: r.get_id()[1]):
        if res.get_id()[0] == " " and "CA" in res:
            atoms.append(res["CA"])
    return atoms

# Global alignment: align all chain B CAs from 8A0N onto 8A0Z
apo_ca_all = get_chain_ca(s_apo, "B")
holo_ca_all = get_chain_ca(s_holo, "B")

# Match on shared res_ids only
apo_ids = {a.get_parent().get_id()[1]: a for a in apo_ca_all}
holo_ids = {a.get_parent().get_id()[1]: a for a in holo_ca_all}
shared_ids = sorted(set(apo_ids.keys()) & set(holo_ids.keys()))
print(f"Shared CA residues for alignment: {len(shared_ids)}")

apo_align = [apo_ids[r] for r in shared_ids]
holo_align = [holo_ids[r] for r in shared_ids]

sup = Superimposer()
sup.set_atoms(holo_align, apo_align)  # fixed=holo, moving=apo
sup.apply(s_apo.get_atoms())  # transform all apo atoms
print(f"Global alignment RMSD (all-CA): {sup.rms:.3f} Å")

# Now compute pocket CA RMSD post-alignment
pocket_apo = [(r, apo_ids[r].get_vector().get_array()) for r in shared_ids if r in pocket_res_ids]
pocket_holo = [(r, holo_ids[r].get_vector().get_array()) for r in shared_ids if r in pocket_res_ids]
pocket_apo_d = {r: c for r, c in pocket_apo}
pocket_holo_d = {r: c for r, c in pocket_holo}

common_pocket = sorted(set(pocket_apo_d.keys()) & set(pocket_holo_d.keys()))
print(f"Pocket residues with CA in both (post-alignment): {len(common_pocket)}")

diffs_ca = np.array([pocket_apo_d[r] - pocket_holo_d[r] for r in common_pocket])
ca_rmsd = np.sqrt((diffs_ca**2).sum(axis=1).mean())
print(f"\nPocket CA RMSD (8A0N vs 8A0Z, post-alignment): {ca_rmsd:.3f} Å")

# All-atom RMSD over pocket residues (post-alignment)
def get_pocket_allatom(model, chain_id, res_ids):
    chain = model[chain_id]
    result = {}
    for res in chain.get_residues():
        rid = res.get_id()[1]
        if rid in res_ids and res.get_id()[0] == " ":
            atom_map = {a.get_name(): a.get_vector().get_array()
                        for a in res.get_atoms() if a.element != "H"}
            result[rid] = atom_map
    return result

apo_aa = get_pocket_allatom(s_apo, "B", pocket_res_ids)
holo_aa = get_pocket_allatom(s_holo, "B", pocket_res_ids)

all_sq = []
for rid in common_pocket:
    apo_atoms = apo_aa.get(rid, {})
    holo_atoms = holo_aa.get(rid, {})
    for aname in set(apo_atoms.keys()) & set(holo_atoms.keys()):
        d = np.array(apo_atoms[aname]) - np.array(holo_atoms[aname])
        all_sq.append((d**2).sum())

aa_rmsd = np.sqrt(np.mean(all_sq))
print(f"Pocket all-atom RMSD (8A0N vs 8A0Z, post-alignment): {aa_rmsd:.3f} Å")
print(f"CA/all-atom ratio: {ca_rmsd/aa_rmsd:.3f}")

required = (1.0, 2.5)
status = "PASS" if required[0] <= ca_rmsd <= required[1] else "FAIL"
print(f"\nPocket RMSD gate ({required[0]}–{required[1]} Å): {status}")

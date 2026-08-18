"""
Feature extraction, native poses only. Groups A-F per Stage 1 spec.

Group A (Vina terms) is reimplemented directly from pdbqt geometry using Vina's
published functional forms (Trott & Olson 2010) -- the official vina4dv fork
that emits per-term breakdowns is unbuildable here (see Phase 7), and the
`vina` python package needs a Boost build this environment doesn't have ready.
This yields UNWEIGHTED per-term sums (not the final linear combination) so the
downstream model can learn its own weighting -- that's the point of delta
learning.

H-bond donor/acceptor typing (groups D, E) is geometric: acceptor = AD4 type in
{NA,OA,SA,...}; donor = a polar heavy atom with an attached HD hydrogen within
1.3A. This is a simplification relative to a full chemistry-aware perception
but is standard practice in scoring-function reimplementations.
"""
import numpy as np
from scipy.spatial import cKDTree
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import freesasa

from .pdbqt_atoms import parse_pdbqt_atoms, atoms_to_arrays, find_donor_heavy_atoms, VDW_RADII, DEFAULT_RADIUS

CUTOFF = 8.0

# ---------------- Group A: Vina terms (unweighted, summed over cross pairs) ----------------

def vina_terms(rec_coords, rec_radii, rec_hydrophobic, rec_acc, rec_don,
                lig_coords, lig_radii, lig_hydrophobic, lig_acc, lig_don):
    tree = cKDTree(rec_coords)
    pairs = tree.query_ball_point(lig_coords, r=CUTOFF)

    gauss1 = gauss2 = repulsion = hydrophobic = hbond = 0.0
    for li, rec_idx in enumerate(pairs):
        if not rec_idx:
            continue
        rec_idx = np.array(rec_idx)
        dxyz = rec_coords[rec_idx] - lig_coords[li]
        dist = np.linalg.norm(dxyz, axis=1)
        d = dist - (rec_radii[rec_idx] + lig_radii[li])  # surface distance

        gauss1 += np.sum(np.exp(-(d / 0.5) ** 2))
        gauss2 += np.sum(np.exp(-((d - 3.0) / 2.0) ** 2))
        repulsion += np.sum(np.where(d < 0, d ** 2, 0.0))

        if lig_hydrophobic[li]:
            hp = rec_hydrophobic[rec_idx]
            dd = d[hp]
            hydrophobic += np.sum(np.clip(1.5 - dd, 0.0, 1.0))

        # hbond: ligand donor-receptor acceptor, or ligand acceptor-receptor donor
        if lig_don[li]:
            acc = rec_acc[rec_idx]
            dd = d[acc]
            hbond += np.sum(np.clip(dd / (-0.7), 0.0, 1.0) * (dd < 0))
        if lig_acc[li]:
            don = rec_don[rec_idx]
            dd = d[don]
            hbond += np.sum(np.clip(dd / (-0.7), 0.0, 1.0) * (dd < 0))

    return {"vina_gauss1": gauss1, "vina_gauss2": gauss2, "vina_repulsion": repulsion,
            "vina_hydrophobic": hydrophobic, "vina_hbond": hbond}

# ---------------- Group C: distance-binned contacts by ligand heavy-atom element ----------------

ELEMENT_GROUPS = {"C": "C", "A": "C", "N": "N", "NA": "N", "NS": "N",
                   "O": "O", "OA": "O", "OS": "O", "S": "S", "SA": "S"}
def element_group(adtype):
    return ELEMENT_GROUPS.get(adtype, "X")

BINS = [(0, 2), (2, 4), (4, 6), (6, 8)]

def contact_counts(rec_coords, rec_types, lig_coords, lig_types):
    tree = cKDTree(rec_coords)
    pairs = tree.query_ball_point(lig_coords, r=8.0)
    counts = {f"contact_{eg}_{lo}_{hi}": 0 for eg in ("C","N","O","S","X") for lo,hi in BINS}
    for li, rec_idx in enumerate(pairs):
        if not rec_idx:
            continue
        eg = element_group(lig_types[li])
        rec_idx = np.array(rec_idx)
        dist = np.linalg.norm(rec_coords[rec_idx] - lig_coords[li], axis=1)
        for lo, hi in BINS:
            n = int(np.sum((dist >= lo) & (dist < hi)))
            counts[f"contact_{eg}_{lo}_{hi}"] += n
    return counts

# ---------------- Groups D & E: H-bonds + satisfaction ----------------

def hbond_and_satisfaction(rec_coords, rec_acc, rec_don,
                            lig_coords, lig_acc, lig_don, lig_types):
    """Geometric H-bond: heavy-atom donor...acceptor distance < 3.5A (no explicit
    donor-H-acceptor angle since we don't reconstruct H positions for the receptor)."""
    lig_don_idx = np.where(lig_don)[0]
    lig_acc_idx = np.where(lig_acc)[0]
    rec_don_idx = np.where(rec_don)[0]
    rec_acc_idx = np.where(rec_acc)[0]

    n_satisfied_donor = 0
    n_satisfied_acceptor = 0
    HBOND_DIST = 3.5

    if len(lig_don_idx) and len(rec_acc_idx):
        tree = cKDTree(rec_coords[rec_acc_idx])
        for li in lig_don_idx:
            d, _ = tree.query(lig_coords[li], k=1)
            if d < HBOND_DIST:
                n_satisfied_donor += 1
    if len(lig_acc_idx) and len(rec_don_idx):
        tree = cKDTree(rec_coords[rec_don_idx])
        for li in lig_acc_idx:
            d, _ = tree.query(lig_coords[li], k=1)
            if d < HBOND_DIST:
                n_satisfied_acceptor += 1

    n_donors = len(lig_don_idx)
    n_acceptors = len(lig_acc_idx)
    n_unsatisfied_donor = n_donors - n_satisfied_donor
    n_unsatisfied_acceptor = n_acceptors - n_satisfied_acceptor
    total_polar = n_donors + n_acceptors
    satisfied = n_satisfied_donor + n_satisfied_acceptor
    polar_satisfaction_fraction = satisfied / total_polar if total_polar else 0.0

    return {
        "hbond_count": satisfied,
        "n_donors": n_donors, "n_acceptors": n_acceptors,
        "buried_unsatisfied_donors": n_unsatisfied_donor,
        "buried_unsatisfied_acceptors": n_unsatisfied_acceptor,
        "polar_satisfaction_fraction": polar_satisfaction_fraction,
    }

# ---------------- Group B: buried SASA (FreeSASA) ----------------

def buried_sasa(ligand_pdbqt_path, complex_pdb_path=None, rec_coords=None, rec_atoms=None,
                 lig_coords=None, lig_atoms=None):
    """Free ligand SASA minus bound ligand SASA, total + by element.
    Uses freesasa directly on coordinate arrays (Shrake-Rupley), constructing a
    minimal PDB-like structure in memory via freesasa's Python API."""
    def make_structure(atoms_subset, coords_subset):
        structure = freesasa.Structure(options={"skip-unknown": False})
        radii = []
        for i, (a, c) in enumerate(zip(atoms_subset, coords_subset)):
            elem = a["adtype"][0] if a["adtype"][0] in "CNOSH" else "C"
            structure.addAtom(f" {elem:<3}", "LIG", str(i % 9999), "L", c[0], c[1], c[2])
            radii.append(VDW_RADII.get(a["adtype"], DEFAULT_RADIUS))
        structure.setRadii(radii)
        return structure

    # free ligand (ligand atoms only)
    free_struct = make_structure(lig_atoms, lig_coords)
    free_result = freesasa.calc(free_struct)
    free_total = free_result.totalArea()

    # bound: ligand + nearby receptor atoms (within 8A, sufficient for burial calc)
    tree = cKDTree(rec_coords) if len(rec_coords) else None
    if tree is not None:
        nearby_idx = set()
        for lc in lig_coords:
            idx = tree.query_ball_point(lc, r=10.0)
            nearby_idx.update(idx)
        nearby_idx = sorted(nearby_idx)
    else:
        nearby_idx = []

    combined_atoms = lig_atoms + [rec_atoms[i] for i in nearby_idx]
    combined_coords = np.vstack([lig_coords, rec_coords[nearby_idx]]) if nearby_idx else lig_coords
    bound_struct = make_structure(combined_atoms, combined_coords)
    bound_result = freesasa.calc(bound_struct)

    # per-atom SASA for just the ligand portion (first len(lig_atoms) atoms)
    lig_bound_total = sum(bound_result.atomArea(i) for i in range(len(lig_atoms)))

    buried_total = free_total - lig_bound_total

    # by element
    by_elem_free = {}
    by_elem_bound = {}
    for i, a in enumerate(lig_atoms):
        eg = element_group(a["adtype"])
        by_elem_free[eg] = by_elem_free.get(eg, 0.0) + free_result.atomArea(i)
        by_elem_bound[eg] = by_elem_bound.get(eg, 0.0) + bound_result.atomArea(i)

    out = {"buried_sasa_total": buried_total}
    for eg in ("C", "N", "O", "S", "X"):
        out[f"buried_sasa_{eg}"] = by_elem_free.get(eg, 0.0) - by_elem_bound.get(eg, 0.0)
    return out

# ---------------- Group F: ligand descriptors (RDKit) ----------------

def ligand_descriptors(sdf_path):
    mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
    if mol is None:
        return None
    molH = Chem.AddHs(mol, addCoords=True) if mol.GetNumAtoms() == mol.GetNumHeavyAtoms() else mol
    return {
        "hac": mol.GetNumHeavyAtoms(),
        "mw": Descriptors.MolWt(mol),
        "clogp": Descriptors.MolLogP(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "ring_count": rdMolDescriptors.CalcNumRings(mol),
        "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(mol),
    }

# ---------------- Orchestration ----------------

def extract_all_features(receptor_pdbqt, ligand_pdbqt, ligand_sdf):
    rec_atoms = parse_pdbqt_atoms(receptor_pdbqt)
    lig_atoms = parse_pdbqt_atoms(ligand_pdbqt)
    if not rec_atoms or not lig_atoms:
        return None

    rec_coords, rec_radii, rec_types, rec_hp, rec_acc, rec_don_heavy = atoms_to_arrays(rec_atoms)
    lig_coords, lig_radii, lig_types, lig_hp, lig_acc, lig_don_heavy = atoms_to_arrays(lig_atoms)
    rec_don = find_donor_heavy_atoms(rec_atoms, rec_coords) & rec_don_heavy
    lig_don = find_donor_heavy_atoms(lig_atoms, lig_coords) & lig_don_heavy

    feats = {}
    feats.update(vina_terms(rec_coords, rec_radii, rec_hp, rec_acc, rec_don,
                             lig_coords, lig_radii, lig_hp, lig_acc, lig_don))
    feats.update(contact_counts(rec_coords, rec_types, lig_coords, lig_types))
    feats.update(hbond_and_satisfaction(rec_coords, rec_acc, rec_don,
                                         lig_coords, lig_acc, lig_don, lig_types))
    feats.update(buried_sasa(ligand_pdbqt, rec_coords=rec_coords, rec_atoms=rec_atoms,
                              lig_coords=lig_coords, lig_atoms=lig_atoms))
    lig_desc = ligand_descriptors(ligand_sdf)
    if lig_desc is None:
        return None
    feats.update(lig_desc)
    return feats

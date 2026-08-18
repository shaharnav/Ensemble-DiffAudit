"""Minimal PDBQT atom parser + AD4/Vina atom-typing helpers, shared by features.py."""
import numpy as np

# AutoDock4/Vina vdW radii (Angstrom), standard table
VDW_RADII = {
    "C": 2.00, "A": 2.00, "N": 1.75, "NA": 1.75, "NS": 1.75,
    "OA": 1.60, "OS": 1.60, "O": 1.60, "SA": 2.00, "S": 2.00,
    "HD": 1.00, "H": 1.00, "F": 1.54, "Cl": 2.04, "Br": 2.16, "I": 2.36,
    "P": 2.10, "Zn": 1.48, "Mg": 1.30, "Ca": 1.97, "Mn": 1.30, "Fe": 1.30,
    "Cu": 1.30, "Na": 1.36, "K": 1.76,
}
DEFAULT_RADIUS = 1.90

HYDROPHOBIC_TYPES = {"C", "A", "F", "Cl", "Br", "I"}
ACCEPTOR_TYPES = {"NA", "OA", "SA", "NS", "OS"}
DONOR_HEAVY_TYPES = {"N", "NA", "NS", "OA", "OS", "O", "SA", "S"}
POLAR_HEAVY_TYPES = ACCEPTOR_TYPES | DONOR_HEAVY_TYPES

def parse_pdbqt_atoms(path):
    """Returns list of dicts: name, resname, x,y,z, adtype."""
    atoms = []
    with open(path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                name = line[12:16].strip()
                resname = line[17:20].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                adtype = line[77:79].strip() if len(line) >= 79 else line.split()[-1]
            except (ValueError, IndexError):
                continue
            atoms.append({"name": name, "resname": resname, "xyz": (x, y, z), "adtype": adtype})
    return atoms

def atoms_to_arrays(atoms):
    """Returns (coords: Nx3 array, radii: N array, types: list, is_hydrophobic: bool array,
    is_acceptor: bool array, is_donor_heavy: bool array)."""
    coords = np.array([a["xyz"] for a in atoms], dtype=float)
    radii = np.array([VDW_RADII.get(a["adtype"], DEFAULT_RADIUS) for a in atoms])
    types = [a["adtype"] for a in atoms]
    is_hydrophobic = np.array([t in HYDROPHOBIC_TYPES for t in types])
    is_acceptor = np.array([t in ACCEPTOR_TYPES for t in types])
    is_donor_heavy = np.array([t in DONOR_HEAVY_TYPES for t in types])
    return coords, radii, types, is_hydrophobic, is_acceptor, is_donor_heavy

def find_donor_heavy_atoms(atoms, coords):
    """A heavy atom (N/O/S family) is a real donor if it has an attached HD within
    bonding distance (~1.2A). Returns boolean array over `atoms`."""
    hd_idx = [i for i, a in enumerate(atoms) if a["adtype"] == "HD"]
    is_donor = np.zeros(len(atoms), dtype=bool)
    if not hd_idx:
        return is_donor
    hd_coords = coords[hd_idx]
    for i, a in enumerate(atoms):
        if a["adtype"] not in DONOR_HEAVY_TYPES:
            continue
        d = np.linalg.norm(hd_coords - coords[i], axis=1)
        if np.any(d < 1.3):
            is_donor[i] = True
    return is_donor

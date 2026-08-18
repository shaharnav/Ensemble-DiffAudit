import os, sys, subprocess
from rdkit import Chem
from meeko import MoleculePreparation

MK_PREPARE_RECEPTOR = os.path.abspath("./venv/bin/mk_prepare_receptor.py")

def prepare_receptor_pdbqt(protein_pdb, out_pdbqt):
    cmd = [
        sys.executable, MK_PREPARE_RECEPTOR,
        "--read_pdb", protein_pdb,
        "-o", out_pdbqt[:-6] if out_pdbqt.endswith(".pdbqt") else out_pdbqt,
        "-p", out_pdbqt,
        "--allow_bad_res",
        "--default_altloc", "A",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode == 0 and os.path.exists(out_pdbqt)

def prepare_ligand_pdbqt_from_sdf(sdf_path, out_pdbqt):
    """Preserves the crystal-pose 3D coordinates (no re-embedding)."""
    mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
    if mol is None:
        return False
    try:
        molH = Chem.AddHs(mol, addCoords=True)
        preparator = MoleculePreparation()
        preparator.prepare(molH, conformer_id=-1)
        pdbqt_string = preparator.write_pdbqt_string()
        with open(out_pdbqt, "w") as f:
            f.write(pdbqt_string)
        return True
    except Exception:
        return False

def ligand_bbox_center_size(sdf_path, padding=8.0):
    mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
    if mol is None:
        return None, None
    conf = mol.GetConformer()
    coords = conf.GetPositions()
    mins, maxs = coords.min(axis=0), coords.max(axis=0)
    center = ((mins + maxs) / 2).tolist()
    size = (maxs - mins + padding).tolist()
    size = [max(s, 20.0) for s in size]  # floor at 20A to give the grid enough room
    return center, size

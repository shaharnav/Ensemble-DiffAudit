import os
import subprocess
import logging
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from Bio.PDB import PDBParser

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

VINA_PATH = os.path.abspath("./bin/vina_1.2.7_mac_aarch64")
MK_PREPARE_RECEPTOR = os.path.abspath("./venv/bin/mk_prepare_receptor.py")

def prepare_receptor(pdb_file, output_pdbqt):
    """
    Prepares a PDB file for docking using Meeko's mk_prepare_receptor (maintained by the
    Vina/AutoDock authors), which handles polar-H merging, correct AD4 typing (backbone
    amide N as donor, only true ring carbons as aromatic 'A', proper acceptor assignment)
    and metal ions natively.
    """
    try:
        import sys
        cmd = [
            sys.executable, MK_PREPARE_RECEPTOR,
            "--read_pdb", pdb_file,
            "-o", output_pdbqt[:-6] if output_pdbqt.endswith(".pdbqt") else output_pdbqt,
            "-p", output_pdbqt,
            "--allow_bad_res",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 or not os.path.exists(output_pdbqt):
            logging.error(f"mk_prepare_receptor failed: {result.stderr}")
            return False

        logging.info("Receptor successfully prepared with Meeko (mk_prepare_receptor).")
        return True
    except Exception as e:
        logging.error(f"Error preparing receptor: {e}")
        return False
def get_center_and_size(pdb_file, target_residue=None):
    """
    Calculates center and size. 
    1. Prioritizes a specific target residue (e.g., "ASN253") if provided.
    2. Falls back to catalytic metals.
    3. Falls back to cavity search.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_file)
    
    # 1. TARGETED RESIDUE SEARCH
    if target_residue:
        # Parse target residue gracefully
        res_name_target = "".join(filter(str.isalpha, target_residue)).upper()
        res_num_str = "".join(filter(str.isdigit, target_residue))
        try:
            res_num_target = int(res_num_str) if res_num_str else None
            target_coords = []
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        res_name = residue.get_resname().strip().upper()
                        res_num = residue.get_id()[1]
                        
                        if res_name == res_name_target and (res_num_target is None or res_num == res_num_target):
                            for atom in residue:
                                target_coords.append(atom.get_coord())
            
            if target_coords:
                logging.info(f"Target residue {target_residue} found. Centering docking box here.")
                coords = np.array(target_coords)
                center = np.mean(coords, axis=0)
                # 15.0A is a standard high-precision box size for a known active site
                return center, [15.0, 15.0, 15.0]
            else:
                logging.warning(f"Target residue {target_residue} not found in PDB. Falling back to Metals/Cavity.")
        except ValueError:
            logging.warning(f"Could not parse target residue format '{target_residue}'. Use format like 'ASN253'.")

    # 2. METAL SEARCH FALLBACK
    metal_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip().upper() in ["ZN", "MG", "MN", "FE", "CO", "CA"]:
                    for atom in residue:
                        metal_coords.append(atom.get_coord())
    
    if metal_coords:
        logging.info("Metal detected - Snapping box to metal center.")
        coords = np.array(metal_coords)
        center = np.mean(coords, axis=0)
        return center, [22.5, 22.5, 22.5]

    logging.info("No targets or metals found - Initializing Cavity Search...")

    # 3. CAVITY SEARCH FALLBACK (Simplified for brevity)
    all_coords = []
    for model in structure:
        for atom in model.get_atoms():
            all_coords.append(atom.get_coord())
            
    if not all_coords:
        return None, None
        
    all_coords_np = np.array(all_coords)
    center = np.mean(all_coords_np, axis=0)
    logging.info("Using geometric center of the protein.")
    return center, [25.0, 25.0, 25.0]

def prepare_ligand(smiles, output_pdbqt):
    """
    Converts SMILES to PDBQT using RDKit and Meeko.
    (Meeko is maintained by the Forli lab, creators of Vina, so it handles ligands perfectly).
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logging.error("Invalid SMILES string.")
            return False
            
        mol = Chem.AddHs(mol)
        embed_status = AllChem.EmbedMolecule(mol, randomSeed=42)
        if embed_status != 0:
            logging.error("RDKit conformer embedding failed for ligand.")
            return False
        AllChem.MMFFOptimizeMolecule(mol)

        preparator = MoleculePreparation()
        preparator.prepare(mol)
        pdbqt_string = preparator.write_pdbqt_string()
        
        with open(output_pdbqt, 'w') as f:
            f.write(pdbqt_string)
            
        return True
    except Exception as e:
        logging.error(f"Error preparing ligand: {e}")
        return False

def run_docking(pdb_file, smiles, output_dir="./results", job_name="job", exhaustiveness=32, target_residue=None, center_coords=None, box_size=None, seed=42, num_modes=9):
    """
    Runs Vina docking.

    Args:
        center_coords: Optional explicit (x, y, z) pocket center. If provided
                       together with box_size, bypasses get_center_and_size().
        box_size:      Optional explicit [sx, sy, sz] search box dimensions.
        seed:          Vina random seed, for reproducibility across runs.
        num_modes:     Number of binding modes Vina should generate.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    receptor_pdbqt = os.path.join(output_dir, f"{job_name}_receptor.pdbqt")
    ligand_pdbqt = os.path.join(output_dir, f"{job_name}_ligand.pdbqt")
    output_docked = os.path.join(output_dir, f"{job_name}_out.pdbqt")
    params_sidecar = os.path.join(output_dir, f"{job_name}_params.json")

    logging.info("Preparing receptor with Meeko...")
    if not prepare_receptor(pdb_file, receptor_pdbqt):
        return None

    logging.info("Preparing ligand with Meeko...")
    if not prepare_ligand(smiles, ligand_pdbqt):
        return None
        
    # Use explicit center/size if provided, otherwise auto-detect
    if center_coords is not None and box_size is not None:
        center = center_coords
        size = box_size
        logging.info(f"Using explicit pocket center: {center}, box: {size}")
    else:
        center, size = get_center_and_size(pdb_file, target_residue=target_residue)
        if center is None:
            logging.error("Could not calculate center/size.")
            return None
        
    logging.info(f"Center: {center}, Size: {size}")

    # Parameters that must match for a cached result to be reused. Any change here
    # (box, exhaustiveness, seed, num_modes) or a newer receptor file invalidates the cache.
    import json
    run_params = {
        "center": [round(float(c), 4) for c in center],
        "size": [round(float(s), 4) for s in size],
        "exhaustiveness": exhaustiveness,
        "seed": seed,
        "num_modes": num_modes,
        "smiles": smiles,
        "receptor_mtime": os.path.getmtime(pdb_file),
    }

    # ── Smart Vina Caching ──
    if os.path.exists(output_docked) and os.path.getsize(output_docked) > 0 and os.path.exists(params_sidecar):
        has_affinity = False
        has_atoms = False
        cached_affinity = None
        try:
            with open(params_sidecar, 'r') as f:
                cached_params = json.load(f)

            if cached_params == run_params:
                with open(output_docked, 'r') as f:
                    for line in f:
                        if line.startswith("REMARK VINA RESULT:"):
                            has_affinity = True
                            parts = line.split()
                            if len(parts) >= 4:
                                cached_affinity = float(parts[3])
                        elif line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ENDMDL"):
                            has_atoms = True

                        if has_affinity and has_atoms:
                            break

                if has_affinity and has_atoms and cached_affinity is not None:
                    logging.info(f"Using cached Vina output for {job_name} (Affinity: {cached_affinity})")
                    return {
                        "affinity": cached_affinity,
                        "docked_file": output_docked,
                        "stdout": "CACHED - Skipped Vina Subprocess"
                    }
                else:
                    logging.warning(f"Found corrupted/incomplete PDBQT cache for {job_name}. Re-running Vina.")
            else:
                logging.info(f"Docking parameters changed for {job_name}. Re-running Vina.")
        except Exception as e:
            logging.warning(f"Error parsing cache for {job_name}: {e}. Re-running Vina.")

    # Run Vina
    cmd = [
        VINA_PATH,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--cpu", "4",
        "--exhaustiveness", str(exhaustiveness),
        "--seed", str(seed),
        "--num_modes", str(num_modes),
        "--out", output_docked
    ]

    logging.info(f"Running Vina: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"Vina failed: {result.stderr}")
            return None

        logging.info("Vina completed successfully.")

        # Parse the best-mode affinity from the output PDBQT's REMARK VINA RESULT line
        # rather than stdout (stdout table rows are fragile to match/re-order).
        best_affinity = None
        with open(output_docked, 'r') as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    if len(parts) >= 4:
                        best_affinity = float(parts[3])
                    break

        with open(params_sidecar, 'w') as f:
            json.dump(run_params, f)

        return {
            "affinity": best_affinity,
            "docked_file": output_docked,
            "stdout": result.stdout
        }

    except Exception as e:
        logging.error(f"Error executing Vina: {e}")
        return None

if __name__ == "__main__":
    # CONTROL TEST: Experimental Structure + Flexible Ligand
    # This will prove the engine works when the biological pocket is actually open.
    
    test_pdb = "./pdbs/P00918_exp.pdb"  # You already have this file
    test_smiles = "CC(=O)Nc1nnc(s1)S(=O)(=O)N" # Acetazolamide (has rotatable bonds)
    
    if os.path.exists(test_pdb):
        print("Running CAII Control Test...")
        # We target the Zinc (ZN) atom specifically.
        res = run_docking(
            test_pdb, 
            test_smiles, 
            job_name="test_acetazolamide", 
            exhaustiveness=32,
            target_residue="ZN"
        )
        print("Control Result:", res)
        
        # Run the analyzer to prove H-Bonds work
        from analyzer import analyze_docking
        if res and res.get("docked_file"):
            analysis = analyze_docking(test_pdb, res["docked_file"])
            print(f"H-Bonds Detected: {analysis['h_bond_count']}")
            print(f"Metal Bonds Detected: {analysis['metal_bond_count']}")
    else:
        print(f"Test PDB {test_pdb} not found. Please ensure it is in the pdbs/ folder.")
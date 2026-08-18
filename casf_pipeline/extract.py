"""
Extract (protein, ligand) from a raw RCSB PDB file with no PDBbind pre-processing.

Ligand identification heuristic (documented, not "official" PDBbind curation --
this is the reproduction-not-replication gap the user flagged):
  - among HETATM groups, exclude a blocklist of crystallization additives/ions/
    solvents/buffers and water
  - pick the largest remaining group by heavy-atom count
  - assign correct bond orders via RDKit's AssignBondOrdersFromTemplate against
    the RCSB Chemical Component Dictionary (CCD) ideal SMILES for that residue
    code (fetched once per unique ligand code, cached)
"""
import os, json, subprocess, warnings
from rdkit import Chem
from rdkit.Chem import AllChem

BLOCKLIST = {
    "HOH","GOL","EDO","PEG","PG4","1PE","PGE","DMS","SO4","PO4","ACT","TRS",
    "MPD","BME","DTT","IPA","MES","EPE","CIT","IMD","FMT","ACY","NA","CL","K",
    "MG","CA","ZN","MN","FE","CU","NI","CO","CD","HG","BR","IOD","NH4","SR",
    "BA","AZI","NO3","GSH","BOG","LDA","P6G","1PG","MRD","CO3","SIN","OCT",
    "OXL","TAM","PLM","MYR","PLC","LMT","PLP","FLC","EOH","MOH","ACE","1CU",
    "SCN","CAC","NO2","GD","TL","CS","LI","RB","EU","YB","LU","W",
    # glycans / sugars commonly attached as glycosylation, not the scored ligand
    "NAG","BMA","MAN","FUC","GAL","NDG","BGC","GLC","XYP","SIA",
    # covalently modified amino acids (appear as HETATM but are part of the protein chain)
    "LLP","TYS","PTR","SEP","TPO","MSE","CSO","CSD","KCX","CGU","MLY",
}

CCD_CACHE_DIR = "results/casf2016/ccd_cache"
os.makedirs(CCD_CACHE_DIR, exist_ok=True)

def get_ccd_smiles(resname):
    """Fetch ideal isomeric SMILES for a CCD component code. Cached to disk."""
    cache_f = os.path.join(CCD_CACHE_DIR, f"{resname}.json")
    if os.path.exists(cache_f):
        with open(cache_f) as f:
            d = json.load(f)
        return d.get("smiles")
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{resname}"
    r = subprocess.run(["curl","-sL","--max-time","15",url], capture_output=True, text=True)
    smiles = None
    try:
        d = json.loads(r.stdout)
        descriptors = d.get("pdbx_chem_comp_descriptor", [])
        canonical = [x for x in descriptors if x.get("type") == "SMILES_CANONICAL"]
        # prefer OpenEye's canonical isomeric SMILES (deterministic tie-break), else CACTVS, else any SMILES
        for program in ("OpenEye OEToolkits", "CACTVS"):
            match = [x["descriptor"] for x in canonical if x.get("program") == program]
            if match:
                smiles = match[0]; break
        if smiles is None:
            any_smiles = [x["descriptor"] for x in descriptors if x.get("type","").startswith("SMILES")]
            smiles = any_smiles[0] if any_smiles else None
    except Exception:
        pass
    with open(cache_f, "w") as f:
        json.dump({"smiles": smiles}, f)
    return smiles

def parse_pdb_hetatm_groups(pdb_path):
    """Group HETATM lines by (chain, resnum, icode, resname)."""
    groups = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            chain = line[21]
            resnum = line[22:26].strip()
            icode = line[26]
            key = (chain, resnum, icode, resname)
            groups.setdefault(key, []).append(line)
    return groups

def pick_ligand_group(groups):
    """Largest non-blocklisted HETATM group by atom count."""
    candidates = [(k, v) for k, v in groups.items() if k[3] not in BLOCKLIST]
    if not candidates:
        return None
    candidates.sort(key=lambda kv: len(kv[1]), reverse=True)
    return candidates[0]

def extract_ligand_mol(pdb_path, out_sdf_path):
    """
    Returns (resname, heavy_atom_count) on success, None on failure.
    Writes an SDF with correct bond orders + crystal 3D coordinates.
    """
    groups = parse_pdb_hetatm_groups(pdb_path)
    picked = pick_ligand_group(groups)
    if picked is None:
        return None
    (chain, resnum, icode, resname), lines = picked

    pdb_block = "".join(lines) + "END\n"
    raw_mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=True)
    if raw_mol is None:
        return None

    smiles = get_ccd_smiles(resname)
    if not smiles:
        return None
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        return None

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = AllChem.AssignBondOrdersFromTemplate(template, raw_mol)
        Chem.SanitizeMol(mol)
    except Exception:
        return None

    mol.SetProp("_Name", resname)
    writer = Chem.SDWriter(out_sdf_path)
    writer.write(mol)
    writer.close()
    return resname, mol.GetNumHeavyAtoms()

def extract_protein(pdb_path, out_pdb_path, ligand_key):
    """Write a protein-only PDB: standard ATOM records + retained metal ions
    (blocklist metals are still real receptor cofactors, e.g. catalytic Zn/Ca --
    only excluded from LIGAND consideration, not stripped from the receptor)."""
    chain, resnum, icode, resname = ligand_key
    # limited to metals meeko's mk_prepare_receptor has covalent radii for
    metal_hetatm_codes = {"ZN","MG","MN","FE","CA","CU","NA","K"}
    with open(pdb_path) as f, open(out_pdb_path, "w") as out:
        for line in f:
            if line.startswith("ATOM"):
                out.write(line)
            elif line.startswith("HETATM"):
                rn = line[17:20].strip()
                if rn in metal_hetatm_codes:
                    out.write(line)
        out.write("END\n")

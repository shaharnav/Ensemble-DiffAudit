"""
GluR smoke test preparation:
  1. Write confirmed pocket_config_glur.yaml
  2. Extract FASTA from 2OHG
  3. Download mmCIF files (2OHG apo, 2OHV reference)
  4. Compute pocket-lining residue subset for ConforMix --subset-residues
  5. Print exact ConforMix command for Colab
"""
import os, requests, time, yaml
import numpy as np
from Bio.PDB import PDBParser
from Bio import SeqIO
from Bio.PDB.Polypeptide import PPBuilder

PDB_DIR = "pdbs/exp3_pocketminer"
CIF_DIR = "pdbs/exp3_pocketminer/cif"
os.makedirs(CIF_DIR, exist_ok=True)

APO_ID, APO_CHAIN   = "2OHG", "A"
HOLO_ID, HOLO_CHAIN = "2OHV", "A"
LIG_CODE = "NHL"

# ── 1. pocket_config_glur.yaml ───────────────────────────────────────────────
# Values from glur_pocket_config.py run (corrected redock with heavy atoms only)
config = {
    "center": [36.096, 35.361, 70.104],
    "box_size_A": [18.828, 17.834, 14.970],
    "transfer_source_pdb": HOLO_ID,
    "transfer_ligand_code": LIG_CODE,
    "local_alignment_rmsd_A": 1.235,
    "local_alignment_n_residues": 49,
    "check_2a_alignment_validity_distance_A": 3.97,
    "check_2a_alignment_validity_threshold_A": 5.0,
    "check_2a_alignment_validity_result": "PASS",
    "check_2b_box_coverage_min_margin_A": 5.00,
    "check_2b_box_coverage_required_margin_A": 5.0,
    "check_2b_box_coverage_result": "PASS",
    "check_2c_redock_validation_distance_A": 4.634,
    "check_2c_redock_validation_threshold_A": 5.0,
    "check_2c_redock_validation_result": "PASS",
    "check_2c_redock_affinity_kcal_mol": -6.914,
    "check_2c_note": "Centroid distance computed over heavy atoms only (21 atoms); "
                     "Meeko v0.7 placed 2 polar Hs at >100 A -- excluded.",
    "receptor_apo_pdb": "results/experiment3_glur/receptors/R_apo_glur.pdb",
    "receptor_pdbqt": "results/experiment3_glur/receptors/R_apo_glur.pdbqt",
    "note": (
        "Box determined by local superposition of 49 pocket-lining residues "
        "(8 A from NHL in 2OHV), NHL heavy-atom centroid transferred to apo frame. "
        "Box covers all 21 NHL heavy atoms with 5.00 A margin. "
        "Validated by NHL redocking into 2OHG (4.634 A to expected centroid, -6.9 kcal/mol). "
        "All downstream stages must read center/box from this file."
    ),
}
yaml_path = "pocket_config_glur.yaml"
with open(yaml_path, "w") as f:
    yaml.dump(config, f, sort_keys=False, default_flow_style=False)
print(f"Written: {yaml_path}")
print(f"  Center:   {config['center']}")
print(f"  Box size: {config['box_size_A']}")
print(f"  ALL GATES: 2a={config['check_2a_alignment_validity_result']}  "
      f"2b={config['check_2b_box_coverage_result']}  "
      f"2c={config['check_2c_redock_validation_result']}")


# ── 2. Pocket-lining residue IDs (for ConforMix --subset-residues) ────────────
parser = PDBParser(QUIET=True)
holo_path = os.path.join(PDB_DIR, f"{HOLO_ID}.pdb")
holo_struct = parser.get_structure(HOLO_ID, holo_path)[0]
holo_chain = holo_struct[HOLO_CHAIN]

nhl_atoms = [a for res in holo_chain.get_residues()
             if res.get_resname() == LIG_CODE
             for a in res.get_atoms() if a.element and a.element != "H"]
nhl_coords = np.array([a.get_vector().get_array() for a in nhl_atoms])

pocket_rids = set()
for ch in holo_struct.get_chains():
    for res in ch.get_residues():
        if res.get_id()[0] != " ": continue
        for atom in res.get_atoms():
            if not atom.element or atom.element == "H": continue
            coord = np.array(atom.get_vector().get_array())
            if np.sqrt(((nhl_coords - coord[None,:])**2).sum(axis=1)).min() <= 8.0:
                pocket_rids.add(res.get_id()[1])
                break

pocket_rids_sorted = sorted(pocket_rids)
print(f"\nPocket-lining residues ({len(pocket_rids_sorted)} total):")
print(f"  {pocket_rids_sorted}")

# Format as compact range string for --subset-residues
def to_ranges(nums):
    nums = sorted(nums)
    ranges = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
        else:
            ranges.append(f"{start}-{prev}" if start != prev else str(start))
            start = prev = n
    ranges.append(f"{start}-{prev}" if start != prev else str(start))
    return ",".join(ranges)

subset_str = to_ranges(pocket_rids_sorted)
print(f"  --subset-residues \"{subset_str}\"")


# ── 3. Extract FASTA from 2OHG ────────────────────────────────────────────────
apo_path = os.path.join(PDB_DIR, f"{APO_ID}.pdb")
apo_struct_f = parser.get_structure(APO_ID, apo_path)[0]
ppb = PPBuilder()
fasta_path = f"pdbs/exp3_pocketminer/{APO_ID}.fasta"
with open(fasta_path, "w") as fh:
    for i, pp in enumerate(ppb.build_peptides(apo_struct_f[APO_CHAIN])):
        seq = str(pp.get_sequence())
        fh.write(f">{APO_ID}_chain{APO_CHAIN}_pp{i+1}\n{seq}\n")

# Concatenate all peptide fragments into one sequence
seqs = [str(pp.get_sequence()) for pp in ppb.build_peptides(apo_struct_f[APO_CHAIN])]
full_seq = "".join(seqs)
# Rewrite as single sequence
with open(fasta_path, "w") as fh:
    fh.write(f">{APO_ID}|Chain{APO_CHAIN}|GluR S.pyogenes\n{full_seq}\n")
print(f"\nFASTA written: {fasta_path}")
print(f"  Sequence length: {len(full_seq)} aa")
print(f"  First 30 aa: {full_seq[:30]}")


# ── 4. Download CIF files ─────────────────────────────────────────────────────
for pdb_id in [APO_ID, HOLO_ID]:
    cif_path = os.path.join(CIF_DIR, f"{pdb_id}.cif")
    if not os.path.exists(cif_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        with open(cif_path, "w") as fh:
            fh.write(r.text)
        time.sleep(0.2)
        print(f"Downloaded: {cif_path}")
    else:
        print(f"Already exists: {cif_path}")


# ── 5. ConforMix command ──────────────────────────────────────────────────────
print(f"""
{'='*65}
CONFORMIX COMMAND (run in Google Colab with GPU):
{'='*65}

# Upload these files to Colab:
#   {fasta_path}
#   {CIF_DIR}/{APO_ID}.cif   (apo — input structure for folding)
#   {CIF_DIR}/{HOLO_ID}.cif  (holo — twist reference)

!python -m boltz.run_conformixrmsd_boltz \\
    --fasta-path /content/{APO_ID}.fasta \\
    --reference-cif /content/{HOLO_ID}.cif \\
    --subset-residues "{subset_str}" \\
    --out-dir /content/conformix_glur_smoke \\
    --twist-target-start 0.0 \\
    --twist-target-stop 2.0 \\
    --num-twist-targets 2 \\
    --samples-per-target 1 \\
    --structured-regions-only

# This produces 2 conformers:
#   conformix_var_0.pdb (twist target 0.0 A from reference)
#   conformix_var_1.pdb (twist target 2.0 A from reference)
#
# Download all conformix_var_*.pdb to:
#   results/experiment3_glur/apo_ensemble/
{'='*65}
""")

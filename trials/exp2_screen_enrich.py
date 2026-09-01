"""
Enrich Experiment 2 target_screen.csv qualifying pairs with:
  - has_cofactor_difference (from existing cofactor_difference column, verified against HETATM)
  - ligand_mw_holo (from RCSB CCD)
  - source_organism (from RCSB entry)
Then filter and rank.
"""
import csv, re, time
import requests

SCREEN_CSV = "target_screen.csv"

def rcsb_chemcomp(comp_id):
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            d = r.json()
            mw = d.get("chem_comp", {}).get("formula_weight")
            name = d.get("chem_comp", {}).get("name", "")
            return float(mw) if mw else None, name
    except Exception:
        pass
    return None, ""

def rcsb_source_organism(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            d = r.json()
            # Try polymer entities for organism
            orgs = d.get("rcsb_entry_info", {}).get("source_organism_common_names", [])
            if orgs:
                return orgs[0]
    except Exception:
        pass
    # Fallback: polymer entity
    url2 = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
    try:
        r = requests.get(url2, timeout=10)
        if r.status_code == 200:
            d = r.json()
            taxa = d.get("rcsb_entity_source_organism", [])
            if taxa:
                return taxa[0].get("scientific_name", "")
    except Exception:
        pass
    return ""

def parse_ligand_code(field):
    """Extract primary ligand CCD code from field like '2xCHD', 'AHK:403', 'DJ3'."""
    field = field.strip().split(",")[0].strip()  # take first if comma-separated
    field = re.sub(r"^\d+x", "", field)          # remove count prefix like '2x'
    field = field.split(":")[0]                  # remove residue suffix like ':403'
    return field

# Read qualifying pairs
qualifying = []
with open(SCREEN_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row.get("qualifies") == "True":
            qualifying.append(row)

print(f"Qualifying pairs: {len(qualifying)}\n")

results = []
for row in qualifying:
    apo = row["apo_pdb"]
    holo = row["holo_pdb"]
    ligand_field = row["ligand_field"]
    lig_code = parse_ligand_code(ligand_field)
    pocket_ca = float(row["pocket_ca_rmsd"])
    ratio = float(row["ca_to_allatom_ratio"])
    motion = row["motion_type"]
    cofactor_diff = row["cofactor_difference"] == "True"  # pre-existing column

    print(f"  {apo}/{holo}  ligand={lig_code}  ...", end=" ", flush=True)
    mw, lig_name = rcsb_chemcomp(lig_code)
    organism = rcsb_source_organism(holo)
    drug_like = (mw is not None) and (150 <= mw <= 600)
    mw_str = f"{mw:.1f}" if mw else "N/A"
    print(f"MW={mw_str}  organism={organism}")
    time.sleep(0.1)

    results.append({
        "apo_pdb": apo,
        "holo_pdb": holo,
        "pocket_ca_rmsd": pocket_ca,
        "ca_to_allatom_ratio": ratio,
        "motion_type": motion,
        "source_organism": organism,
        "has_cofactor_difference": cofactor_diff,
        "ligand_code": lig_code,
        "ligand_name": lig_name[:60] if lig_name else "",
        "ligand_mw_holo": mw,
        "drug_like_mw": drug_like,
    })

print("\n\n=== Filtered list (cofactor_diff=False, drug-like MW), ranked by pocket CA RMSD ===")
filtered = [r for r in results if not r["has_cofactor_difference"] and r["drug_like_mw"]]
filtered.sort(key=lambda x: x["pocket_ca_rmsd"], reverse=True)

print(f"\n{'Pair':<15} {'CA RMSD':>8} {'CA/AA':>6} {'Motion':<12} {'MW':>7} {'Organism':<30} {'Ligand'}")
print("-" * 110)
for r in filtered:
    pair = f"{r['apo_pdb']}/{r['holo_pdb']}"
    print(f"{pair:<15} {r['pocket_ca_rmsd']:>8.3f} {r['ca_to_allatom_ratio']:>6.3f} "
          f"{r['motion_type']:<12} {r['ligand_mw_holo']:>7.1f} "
          f"{r['source_organism']:<30} {r['ligand_code']} — {r['ligand_name'][:40]}")

print(f"\nTotal passing pairs: {len(filtered)}")

print("\n\n=== Excluded (cofactor_diff or non-drug-like MW) ===")
excluded = [r for r in results if r["has_cofactor_difference"] or not r["drug_like_mw"]]
for r in excluded:
    reason = []
    if r["has_cofactor_difference"]: reason.append("cofactor_diff")
    if not r["drug_like_mw"]: reason.append(f"MW={r['ligand_mw_holo']}")
    pair = f"{r['apo_pdb']}/{r['holo_pdb']}"
    print(f"  {pair:<15} {r['ligand_code']}  EXCLUDED: {', '.join(reason)}")

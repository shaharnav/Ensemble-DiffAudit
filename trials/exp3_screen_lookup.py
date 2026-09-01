"""Look up protein identity for 4 PocketMiner pairs + holo gates for 3 new candidates."""
import requests, json, time

def rcsb_entry(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
    r = requests.get(url, timeout=15)
    return r.json() if r.status_code == 200 else {}

def rcsb_polymer_entity(pdb_id, entity_id=1):
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.upper()}/{entity_id}"
    r = requests.get(url, timeout=15)
    return r.json() if r.status_code == 200 else {}

def rcsb_chemcomp(comp_id):
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id.upper()}"
    r = requests.get(url, timeout=10)
    if r.status_code == 200:
        d = r.json().get("chem_comp", {})
        return d.get("formula_weight"), d.get("name", ""), d.get("type", "")
    return None, "", ""

def uniprot_search_pdbs(uniprot_id):
    """Get all PDB IDs mapped to a UniProt accession."""
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{uniprot_id}"
    r = requests.get(url, timeout=20)
    if r.status_code == 200:
        data = r.json().get(uniprot_id, {}).get("PDB", {})
        return list(data.keys())
    return []

def get_struct_title_organism(pdb_id):
    d = rcsb_entry(pdb_id)
    title = d.get("struct", {}).get("title", "")
    e = rcsb_polymer_entity(pdb_id, 1)
    org = ""
    taxa = e.get("rcsb_entity_source_organism", [])
    if taxa:
        org = taxa[0].get("scientific_name", "")
    desc = e.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
    return title, desc, org

def hetatm_ligands(pdb_id, drug_only=True):
    """Return non-water, non-standard HETATM residue names from RCSB nonpolymer entities."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
    d = requests.get(url, timeout=15).json() if True else {}
    # Use chem_comp ligands from nonpolymer entities
    url2 = f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id.upper()}/1"
    # Actually iterate entity IDs
    results = []
    for eid in range(1, 10):
        url3 = f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id.upper()}/{eid}"
        r = requests.get(url3, timeout=10)
        if r.status_code != 200:
            break
        data = r.json()
        comp_id = data.get("pdbx_entity_nonpoly", {}).get("comp_id", "")
        if comp_id and comp_id not in ("HOH", "WAT", "DOD"):
            mw, name, ctype = rcsb_chemcomp(comp_id)
            results.append({"comp_id": comp_id, "mw": mw, "name": name, "type": ctype})
        time.sleep(0.05)
    return results


# ── Part 1: PocketMiner pair protein identities ──────────────────────────────
pairs = [
    ("5NIA", "5NI6", "DJ3", "ssm,loop", "H. sapiens"),
    ("3NX1", "3NX2", "FER", "ssm",      "Enterobacter sp."),
    ("2OHG", "2OHV", "NHL", "loop",     "S. pyogenes"),
    ("1EZM", "3DBK", "RDF", "id",       "P. aeruginosa"),
]

print("=" * 70)
print("PART 1 — PocketMiner pair protein identities")
print("=" * 70)
for apo, holo, lig, motion, org in pairs:
    title, desc, orgsci = get_struct_title_organism(holo)
    print(f"\n{apo}/{holo} ({lig}):")
    print(f"  Title:    {title}")
    print(f"  Protein:  {desc}")
    print(f"  Organism: {orgsci}")
    time.sleep(0.15)


# ── Part 2: Holo gate for 3 new candidates ───────────────────────────────────
print("\n\n" + "=" * 70)
print("PART 2 — Holo gate for new candidates")
print("=" * 70)

BUFFER_CODES = {
    "SO4","PO4","GOL","EDO","PEG","MPD","BME","DTT","CIT","ACT","ACE",
    "TRS","MES","HEP","TAR","TLA","FMT","IMD","EPE","PGE","PG4","1PE",
    "15P","EOH","ETF","IOD","CL","NA","K","MG","CA","ZN","MN","FE","CU",
    "NI","CO","CD","HG","PT","AU","BR","I","F","IPA","DMS","GLC","FRU",
    "GAL","MAN","XYL","RIB","FUC","NAG","NDG","BMA","MMA","SEP","PTR",
    "TPO","CSO","HOH","WAT","DOD"
}
COFACTOR_CODES = {
    "NAD","NADH","NAP","NADP","FAD","FMN","ATP","ADP","AMP","GTP","GDP",
    "TDP","TTP","SAM","SAH","COA","HEM","HEC","CLA","BCL","SF4","FES",
    "FEO","F3S","F4S","NDP","ANP","AGS","APC","MGT","GNP","POP","FOZ"
}

# LdG6PD — search by UniProt Q9NEB1 (Ld G6PD) or keyword
print("\n-- LdG6PD (Leishmania donovani G6PD) --")
print("Checking known series 7ZHT-7ZHZ + broader search...")
ldg6pd_pdbs = ["7ZHT","7ZHU","7ZHV","7ZHW","7ZHX","7ZHY","7ZHZ"]
# Also try full-text search via RCSB
search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
query = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {"type": "terminal", "service": "full_text",
             "parameters": {"value": "Leishmania donovani glucose-6-phosphate dehydrogenase"}},
        ]
    },
    "return_type": "entry",
    "request_options": {"results_verbosity": "compact", "paginate": {"start": 0, "rows": 50}}
}
r = requests.post(search_url, json=query, timeout=20)
if r.status_code == 200:
    hits = [h["identifier"] for h in r.json().get("result_set", [])]
    print(f"  RCSB full-text hits: {hits}")
    ldg6pd_pdbs = list(set(ldg6pd_pdbs + hits))
else:
    print(f"  RCSB search failed ({r.status_code}), using known series only")

for pdb in sorted(ldg6pd_pdbs):
    ligs = hetatm_ligands(pdb)
    drug_ligs = [l for l in ligs if l["comp_id"] not in BUFFER_CODES and
                 l["comp_id"] not in COFACTOR_CODES and
                 l["mw"] and 150 <= l["mw"] <= 600]
    tag = "DRUG-LIKE" if drug_ligs else ""
    lig_str = "; ".join(f"{l['comp_id']} {l['mw']:.0f} Da" for l in ligs[:6])
    print(f"  {pdb}: {lig_str[:80]}  {tag}")
    time.sleep(0.1)

# LdNMT — Leishmania donovani NMT
print("\n-- LdNMT (Leishmania donovani N-myristoyltransferase) --")
query2 = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {"type": "terminal", "service": "full_text",
             "parameters": {"value": "Leishmania donovani N-myristoyltransferase"}},
        ]
    },
    "return_type": "entry",
    "request_options": {"results_verbosity": "compact", "paginate": {"start": 0, "rows": 50}}
}
r2 = requests.post(search_url, json=query2, timeout=20)
ldnmt_pdbs = []
if r2.status_code == 200:
    ldnmt_pdbs = [h["identifier"] for h in r2.json().get("result_set", [])]
    print(f"  RCSB full-text hits: {ldnmt_pdbs}")
for pdb in sorted(ldnmt_pdbs):
    ligs = hetatm_ligands(pdb)
    drug_ligs = [l for l in ligs if l["comp_id"] not in BUFFER_CODES and
                 l["comp_id"] not in COFACTOR_CODES and
                 l["mw"] and 150 <= l["mw"] <= 600]
    tag = "DRUG-LIKE" if drug_ligs else ""
    lig_str = "; ".join(f"{l['comp_id']} {l['mw']:.0f} Da" for l in ligs[:6])
    print(f"  {pdb}: {lig_str[:80]}  {tag}")
    time.sleep(0.1)

# MtAlaDH — Mycobacterium tuberculosis alanine dehydrogenase
print("\n-- MtAlaDH (M. tuberculosis alanine dehydrogenase) --")
query3 = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {"type": "terminal", "service": "full_text",
             "parameters": {"value": "Mycobacterium tuberculosis alanine dehydrogenase"}},
        ]
    },
    "return_type": "entry",
    "request_options": {"results_verbosity": "compact", "paginate": {"start": 0, "rows": 50}}
}
r3 = requests.post(search_url, json=query3, timeout=20)
mtaladh_pdbs = []
if r3.status_code == 200:
    mtaladh_pdbs = [h["identifier"] for h in r3.json().get("result_set", [])]
    print(f"  RCSB full-text hits: {mtaladh_pdbs}")
for pdb in sorted(mtaladh_pdbs):
    ligs = hetatm_ligands(pdb)
    drug_ligs = [l for l in ligs if l["comp_id"] not in BUFFER_CODES and
                 l["comp_id"] not in COFACTOR_CODES and
                 l["mw"] and 150 <= l["mw"] <= 600]
    tag = "DRUG-LIKE" if drug_ligs else ""
    lig_str = "; ".join(f"{l['comp_id']} {l['mw']:.0f} Da" for l in ligs[:6])
    print(f"  {pdb}: {lig_str[:80]}  {tag}")
    time.sleep(0.1)

print("\nDone.")

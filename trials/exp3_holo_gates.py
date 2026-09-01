"""Holo gate for LdG6PD, LdNMT, MtAlaDH — uses PDBe ligand API (more reliable)."""
import requests, time

BUFFER_CODES = {
    "SO4","PO4","GOL","EDO","PEG","MPD","BME","DTT","CIT","ACT","ACE",
    "TRS","MES","HEP","TAR","TLA","FMT","IMD","EPE","PGE","PG4","1PE",
    "15P","EOH","ETF","IOD","CL","NA","K","MG","CA","ZN","MN","FE","CU",
    "NI","CO","CD","HG","PT","AU","BR","I","F","IPA","DMS","GLC","FRU",
    "GAL","MAN","XYL","RIB","FUC","NAG","NDG","BMA","MMA","SEP","PTR",
    "TPO","CSO","HOH","WAT","DOD","ACY","NH4","PYR","ALA","GLU","G6P",
    "6PG","6PGL","PGS","MOH","BOG","OGA","TLA","FOR","MLI","MPO"
}
COFACTOR_CODES = {
    "NAD","NADH","NAP","NADP","FAD","FMN","ATP","ADP","AMP","GTP","GDP",
    "TDP","TTP","SAM","SAH","COA","HEM","HEC","CLA","BCL","SF4","FES",
    "FEO","F3S","F4S","NDP","ANP","AGS","APC","MGT","GNP","POP","FOZ",
    "MYR"  # myristic acid — NMT substrate
}

def pdbe_ligands(pdb_id):
    """Get all ligands for a PDB entry via PDBe API."""
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb_id.lower()}"
    r = requests.get(url, timeout=15)
    if r.status_code != 200:
        return []
    data = r.json().get(pdb_id.lower(), [])
    seen = set()
    results = []
    for lig in data:
        chem_comp_id = lig.get("chem_comp_id", "")
        if chem_comp_id and chem_comp_id not in seen:
            seen.add(chem_comp_id)
            results.append({
                "comp_id": chem_comp_id,
                "name": lig.get("chem_comp_name", ""),
                "mw": None  # will fill below if needed
            })
    return results

def rcsb_chemcomp(comp_id):
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
    r = requests.get(url, timeout=10)
    if r.status_code == 200:
        d = r.json().get("chem_comp", {})
        return d.get("formula_weight"), d.get("name", "")
    return None, ""

def screen_drug_like(ligs):
    drug = []
    for l in ligs:
        cid = l["comp_id"]
        if cid in BUFFER_CODES or cid in COFACTOR_CODES:
            continue
        if cid in ("HOH","WAT","DOD"):
            continue
        mw = l.get("mw")
        if mw is None:
            mw, name = rcsb_chemcomp(cid)
            if name and not l.get("name"):
                l["name"] = name
            l["mw"] = mw
            time.sleep(0.05)
        if mw and 150 <= mw <= 600:
            drug.append(l)
    return drug

def rcsb_search(query_text):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    payload = {
        "query": {"type": "terminal", "service": "full_text",
                  "parameters": {"value": query_text}},
        "return_type": "entry",
        "request_options": {"results_verbosity": "compact",
                            "paginate": {"start": 0, "rows": 60}}
    }
    r = requests.post(url, json=payload, timeout=20)
    if r.status_code != 200:
        return []
    out = []
    for item in r.json().get("result_set", []):
        if isinstance(item, str):
            out.append(item)
        elif isinstance(item, dict):
            pid = item.get("identifier") or item.get("entry_id") or ""
            if pid:
                out.append(pid)
    return out

def check_holo(label, pdbs, substrate_note=""):
    print(f"\n{'='*65}")
    print(f"{label}")
    if substrate_note:
        print(f"  NOTE: {substrate_note}")
    print(f"  Total structures: {len(pdbs)}")
    print(f"{'='*65}")
    drug_structs = []
    for pdb in sorted(pdbs):
        ligs = pdbe_ligands(pdb)
        if not ligs:
            print(f"  {pdb}: no heteroatoms")
            time.sleep(0.05)
            continue
        drug_ligs = screen_drug_like(ligs)
        all_codes = [l["comp_id"] for l in ligs]
        if drug_ligs:
            d_str = "; ".join(f"{l['comp_id']} {l['mw']:.0f} Da" for l in drug_ligs)
            print(f"  {pdb}: DRUG-LIKE >> {d_str}")
            drug_structs.append((pdb, drug_ligs))
        else:
            print(f"  {pdb}: {', '.join(all_codes)}")
        time.sleep(0.08)
    if drug_structs:
        print(f"\n  HOLO GATE: PASS ({len(drug_structs)} inhibitor-bound structure(s))")
        for pdb, dl in drug_structs:
            for l in dl:
                print(f"    {pdb}: {l['comp_id']} {l['mw']:.0f} Da — {l['name'][:60]}")
    else:
        print(f"\n  HOLO GATE: FAIL")
    return drug_structs


# ── LdG6PD ──────────────────────────────────────────────────────────────────
known_ldg6pd = ["7ZHT","7ZHU","7ZHV","7ZHW","7ZHX","7ZHY","7ZHZ","9VP5","9VP6","9VP7"]
hits = rcsb_search("Leishmania donovani glucose-6-phosphate dehydrogenase")
ldg6pd_pdbs = sorted(set(known_ldg6pd + hits))
check_holo("LdG6PD — Ld G6PD", ldg6pd_pdbs,
           "substrates (G6P, 6PG) and cofactor (NADP) excluded from drug-like filter")

# ── LdNMT ───────────────────────────────────────────────────────────────────
hits2a = rcsb_search("Leishmania donovani N-myristoyltransferase")
hits2b = rcsb_search("Leishmania N-myristoyltransferase inhibitor")
hits2c = rcsb_search("Leishmania NMT DDD")
ldnmt_pdbs = sorted(set(hits2a + hits2b + hits2c))
check_holo("LdNMT — Ld NMT", ldnmt_pdbs,
           "myristic acid (MYR, the acyl substrate) excluded as cofactor")

# ── MtAlaDH ─────────────────────────────────────────────────────────────────
hits3 = rcsb_search("Mycobacterium tuberculosis alanine dehydrogenase")
hits3b = rcsb_search("M. tuberculosis AlaDH inhibitor")
mtaladh_pdbs = sorted(set(hits3 + hits3b))
check_holo("MtAlaDH — Mtb alanine dehydrogenase", mtaladh_pdbs,
           "pyruvate (PYR) and alanine (ALA) excluded as substrates")

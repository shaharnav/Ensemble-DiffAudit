"""
Better pair selection for LdNMT and MtAlaDH.
LdNMT: find structure with MYA only (no inhibitor) to match 6QD9 cofactor environment.
MtAlaDH: find pairs from same crystal form (3ISJ series); check which apo-like
         structure aligns with 3ISJ before running gates.
"""
import os, requests, time
import numpy as np
from Bio.PDB import PDBParser, Superimposer

PDB_DIR = "pdbs/exp3_new"
os.makedirs(PDB_DIR, exist_ok=True)

def download_pdb(pdb_id):
    path = os.path.join(PDB_DIR, f"{pdb_id.upper()}.pdb")
    if not os.path.exists(path):
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", timeout=30)
        r.raise_for_status()
        with open(path, "w") as fh: fh.write(r.text)
        time.sleep(0.2)
    return path

def pdbe_ligands(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb_id.lower()}"
    r = requests.get(url, timeout=15)
    if r.status_code != 200: return []
    data = r.json().get(pdb_id.lower(), [])
    seen = set()
    out = []
    for lig in data:
        cid = lig.get("chem_comp_id","")
        if cid and cid not in seen:
            seen.add(cid)
            out.append({"comp_id": cid, "chain": lig.get("chain_id",""), "name": lig.get("chem_comp_name","")})
    return out

def global_ca_rmsd(pdb1, chain1, pdb2, chain2, max_res=400):
    parser = PDBParser(QUIET=True)
    s1 = parser.get_structure("s1", pdb1)[0]
    s2 = parser.get_structure("s2", pdb2)[0]
    def get_ca(model, cid):
        try: ch = model[cid]
        except KeyError: ch = list(model.get_chains())[0]
        return {res.get_id()[1]: res["CA"] for res in ch.get_residues()
                if res.get_id()[0] == " " and "CA" in res}
    ca1 = get_ca(s1, chain1)
    ca2 = get_ca(s2, chain2)
    shared = sorted(set(ca1) & set(ca2))[:max_res]
    if len(shared) < 30:
        return 999.0, len(shared)
    sup = Superimposer()
    sup.set_atoms([ca2[r] for r in shared], [ca1[r] for r in shared])
    return sup.rms, len(shared)

SKIP = {"HOH","WAT","DOD","SO4","PO4","GOL","EDO","PEG","MPD","BME","CIT","ACT","ACE",
        "TRS","MES","EPE","PG4","1PE","EOH","IPA","DMS","MG","CA","ZN","NA","K",
        "MN","CL","FE","CU","NI","CO","NH4","ACY","PYR","ALA","GLU","IMD"}

print("=" * 65)
print("LdNMT — searching for apo with MYA (myristoyl-CoA analog) only")
print("=" * 65)

# Try structures that may have MYA only (no inhibitor)
nmt_candidates = ["2WUU","4CGL","4CGM","4CGN","4CGO","4CYN","4CYO","4CYP","4CYQ",
                  "5A27","5A28","5G20","5G21","6EU5","6EWF",
                  "6QD9","6QDA","6QDB","6QDC","6QDD","6QDE","6QDF","6QDG","6QDH",
                  "6TW7","6TW8","9RLU","9RLV","9RLW","9RLX","9RLY"]
DRUG_MW_MIN = 150
DRUG_MW_MAX = 600
COFACTORS_NMT = {"NHW","MYA","MYS","COA","AMP","ADP","SCS","MYT"}  # NMT CoA-related

def rcsb_mw(comp_id):
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
    r = requests.get(url, timeout=10)
    if r.status_code == 200:
        return r.json().get("chem_comp", {}).get("formula_weight")
    return None

for pdb in nmt_candidates[:15]:  # spot-check first 15
    ligs = pdbe_ligands(pdb)
    drug_ligs = []
    coa_ligs  = []
    for l in ligs:
        cid = l["comp_id"]
        if cid in SKIP: continue
        if cid in COFACTORS_NMT:
            coa_ligs.append(cid)
        else:
            mw = rcsb_mw(cid)
            if mw and DRUG_MW_MIN <= mw <= DRUG_MW_MAX:
                drug_ligs.append(f"{cid}({mw:.0f})")
        time.sleep(0.05)
    tag = ""
    if coa_ligs and not drug_ligs:
        tag = " << APO-LIKE (CoA only)"
    elif not coa_ligs and not drug_ligs:
        tag = " << APO-LIKE (empty)"
    print(f"  {pdb}: CoA={coa_ligs}  drug={drug_ligs}{tag}")

print("\n\nChecking global alignment for promising apo candidates vs 6QD9 (holo):")
holo_path = download_pdb("6QD9")
for pdb_id in ["2WUU"]:
    apo_path = download_pdb(pdb_id)
    rmsd, n = global_ca_rmsd(apo_path, "A", holo_path, "A")
    print(f"  {pdb_id} vs 6QD9: global CA RMSD = {rmsd:.3f} Å  (n={n})")


print("\n\n" + "=" * 65)
print("MtAlaDH — finding pairs within same crystal form as 3ISJ")
print("=" * 65)

# 3ISJ is one drug-bound structure. Need an apo-like from the same series.
# Check 3IOB, 3IUB, 3LE8, etc. - these are all in the drug-bound series but
# maybe there's a structure with just the NAD cofactor and no drug.
# Also need to check what NAD-related species each has.

ALADH_DRUG_SERIES = ["3COW","3COY","3COZ","3IMC","3IME","3IMG",
                     "3IOB","3IOC","3IOD","3IOE","3ISJ","3IUB","3IUE",
                     "3IVC","3IVG","3IVX","3LE8",
                     "4DDH","4DDK","4DDM","4DE5","4EF6","4FZJ","4G5F","4G5Y",
                     "4LMP","4MQ6","4MUE","4MUF","4MUG","4MUH","4MUI",
                     "4MUJ","4MUK","4MUL","4MUN"]
NAD_CODES = {"NAD","NADH","NAP","NAI","NDP"}
SUBSTRATE_CODES = {"PYR","ALA","2OP","LAC","PDO"}  # pyruvate/lactate/alanine

print("Scanning for apo-like (NAD only, no drug):")
for pdb in ALADH_DRUG_SERIES[:20]:
    ligs = pdbe_ligands(pdb)
    nad = [l["comp_id"] for l in ligs if l["comp_id"] in NAD_CODES]
    sub = [l["comp_id"] for l in ligs if l["comp_id"] in SUBSTRATE_CODES]
    drug = []
    for l in ligs:
        cid = l["comp_id"]
        if cid in SKIP or cid in NAD_CODES or cid in SUBSTRATE_CODES: continue
        mw = rcsb_mw(cid)
        if mw and DRUG_MW_MIN <= mw <= DRUG_MW_MAX:
            drug.append(f"{cid}({mw:.0f})")
        time.sleep(0.04)
    tag = " << DRUG" if drug else (" << APO-like (NAD only)" if nad and not sub else "")
    print(f"  {pdb}: NAD={nad} sub={sub} drug={drug}{tag}")

print("\nChecking alignment of apo candidates (1MOP, 2VOE, 2VHY) vs 3ISJ:")
holo_3isj = download_pdb("3ISJ")
for apo_id, chain in [("1MOP","A"),("2VOE","A"),("2VHY","A"),("2VHX","A"),("2A86","A")]:
    path = download_pdb(apo_id)
    for hchain in ["A","B"]:
        rmsd, n = global_ca_rmsd(path, chain, holo_3isj, hchain)
        if rmsd < 50:
            print(f"  {apo_id}:{chain} vs 3ISJ:{hchain} → {rmsd:.3f} Å (n={n})")
            break

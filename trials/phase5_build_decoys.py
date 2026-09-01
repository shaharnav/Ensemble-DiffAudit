"""
Phase 5 Control B — Decoy Builder (PubChem source)
====================================================
Scans PubChem CID ranges in batches to find property-matched decoys for each
of the 20 Phase 4 candidates. All property filtering is done locally with RDKit.

Matching criteria (per pre-registration):
  MW      ±25 Da
  cLogP   ±0.5
  HBD     ±1
  HBA     ±1
  RotBonds ±2
  Tanimoto < 0.3 (Morgan r=2, 2048 bits) vs all 20 candidates

Source: PubChem PUG REST (https://pubchem.ncbi.nlm.nih.gov/rest/pug)

Writes:
  results/phase5_decoys.sdf
  results/phase5_decoy_matching.csv
"""

import csv, json, time, requests
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, DataStructs

PUBCHEM_PROP_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}"
    "/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,"
    "RotatableBondCount,CanonicalSMILES,IsomericSMILES/JSON"
)

# ── Candidate properties ─────────────────────────────────────────────────
with open('results/lasb_clean_results.json') as f:
    phase4 = json.load(f)

def get_props(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
    return {
        'mw':   Descriptors.ExactMolWt(mol),
        'logp': Descriptors.MolLogP(mol),
        'hbd':  rdMolDescriptors.CalcNumHBD(mol),
        'hba':  rdMolDescriptors.CalcNumHBA(mol),
        'rotb': rdMolDescriptors.CalcNumRotatableBonds(mol),
        'fp': fp, 'mol': mol,
    }

cand_props = [get_props(c['smiles']) for c in phase4]
cand_fps   = [p['fp'] for p in cand_props if p]

def max_tanimoto(fp):
    return max(DataStructs.BulkTanimotoSimilarity(fp, cand_fps))

def fetch_pubchem_batch(cid_start, cid_end):
    cids = ",".join(str(c) for c in range(cid_start, cid_end + 1))
    url = PUBCHEM_PROP_URL.format(cids=cids)
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 404:
            return []
        r.raise_for_status()
        return r.json().get('PropertyTable', {}).get('Properties', [])
    except Exception as e:
        print(f"  PubChem error [{cid_start}-{cid_end}]: {e}")
        return []

TANIMOTO_THRESHOLD = 0.3
BATCH = 500   # CIDs per request

# We need one decoy per candidate; build a pool of candidates still needing matches
unmatched = list(range(len(phase4)))
decoys    = [None] * len(phase4)
used_cids = set()

print(f"{'#':<4} {'cand_MW':>8} {'cand_logP':>9} {'HBD':>4} {'HBA':>4} {'Rb':>4}  "
      f"| {'dec_MW':>7} {'dec_logP':>8} {'Tan':>6}  PubChem_CID")
print("-"*95)

# Scan CID ranges; early CIDs are densely populated with small drug-like molecules
cid = 1
MAX_CID = 200_000

while unmatched and cid <= MAX_CID:
    batch = fetch_pubchem_batch(cid, min(cid + BATCH - 1, MAX_CID))
    time.sleep(0.22)  # PubChem rate limit ~5 req/s

    for entry in batch:
        pcid = entry.get('CID')
        if pcid in used_cids:
            continue
        smi = entry.get('SMILES') or entry.get('IsomericSMILES') or entry.get('CanonicalSMILES', '')
        if not smi:
            continue
        dp = get_props(smi)
        if dp is None:
            continue

        # Try to match each still-unmatched candidate
        for i in list(unmatched):
            cp = cand_props[i]
            if cp is None:
                unmatched.remove(i)
                continue
            if abs(dp['mw']   - cp['mw'])   > 25:  continue
            if abs(dp['logp'] - cp['logp']) > 0.5: continue
            if abs(dp['hbd']  - cp['hbd'])  > 1:   continue
            if abs(dp['hba']  - cp['hba'])  > 1:   continue
            if abs(dp['rotb'] - cp['rotb']) > 2:   continue
            tan = max_tanimoto(dp['fp'])
            if tan >= TANIMOTO_THRESHOLD:           continue

            # Match found
            cp2 = cand_props[i]
            decoys[i] = {
                'cand_idx': i, 'cand_smiles': phase4[i]['smiles'],
                'pubchem_cid': pcid, 'decoy_smiles': smi,
                'cand_mw': cp2['mw'], 'dec_mw': dp['mw'],
                'cand_logp': cp2['logp'], 'dec_logp': dp['logp'],
                'cand_hbd': cp2['hbd'], 'dec_hbd': dp['hbd'],
                'cand_hba': cp2['hba'], 'dec_hba': dp['hba'],
                'cand_rotb': cp2['rotb'], 'dec_rotb': dp['rotb'],
                'max_tanimoto': tan, 'mol': dp['mol'],
            }
            used_cids.add(pcid)
            unmatched.remove(i)
            print(f"{i:<4} {cp2['mw']:>8.1f} {cp2['logp']:>9.2f} {cp2['hbd']:>4} "
                  f"{cp2['hba']:>4} {cp2['rotb']:>4}  | "
                  f"{dp['mw']:>7.1f} {dp['logp']:>8.2f} {tan:>6.3f}  CID{pcid}")
            break  # one decoy per CID

    cid += BATCH

# Report failures
for i in unmatched:
    cp = cand_props[i]
    print(f"{i:<4} {cp['mw']:>8.1f} {cp['logp']:>9.2f} {cp['hbd']:>4} "
          f"{cp['hba']:>4} {cp['rotb']:>4}  | NO MATCH FOUND (scanned CIDs 1–{MAX_CID})")

found = [d for d in decoys if d is not None]
print(f"\nMatched {len(found)}/20 candidates")

# Write SDF
writer = Chem.SDWriter('results/phase5_decoys.sdf')
for d in found:
    mol = d['mol']
    mol.SetProp('_Name', f"CID{d['pubchem_cid']}")
    mol.SetProp('cand_idx', str(d['cand_idx']))
    writer.write(mol)
writer.close()
print(f"Wrote results/phase5_decoys.sdf")

# Write matching table
fields = ['cand_idx','pubchem_cid','cand_smiles','decoy_smiles',
          'cand_mw','dec_mw','cand_logp','dec_logp',
          'cand_hbd','dec_hbd','cand_hba','dec_hba',
          'cand_rotb','dec_rotb','max_tanimoto']
with open('results/phase5_decoy_matching.csv','w',newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
    w.writeheader()
    w.writerows(found)
print(f"Wrote results/phase5_decoy_matching.csv")

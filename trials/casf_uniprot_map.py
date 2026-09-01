"""Fetch UniProt accession for each CASF complex's protein, for family-grouping
(leave-one-protein-family-out CV, since PDBbind's own Pfam clusters aren't
available without registration)."""
import csv, json, subprocess, time
from multiprocessing import Pool

def fetch_uniprot(pdbid):
    query = '{ entry(entry_id: "%s") { polymer_entities { rcsb_polymer_entity_container_identifiers { uniprot_ids } } } }' % pdbid.upper()
    r = subprocess.run(
        ["curl","-s","--max-time","15","-X","POST","https://data.rcsb.org/graphql",
         "-H","Content-Type: application/json","-d", json.dumps({"query": query})],
        capture_output=True, text=True
    )
    try:
        d = json.loads(r.stdout)
        entities = d["data"]["entry"]["polymer_entities"]
        uniprots = set()
        for e in entities:
            ids = e["rcsb_polymer_entity_container_identifiers"].get("uniprot_ids") or []
            uniprots.update(ids)
        return pdbid, ";".join(sorted(uniprots)) if uniprots else None
    except Exception:
        return pdbid, None

if __name__ == "__main__":
    with open("results/casf2016/casf2016_core_pdbid_pkd.csv") as f:
        ids = [r["pdbid"] for r in csv.DictReader(f)]

    t0 = time.time()
    with Pool(processes=12) as pool:
        results = pool.map(fetch_uniprot, ids)
    print(f"Fetched {sum(1 for _,u in results if u)}/{len(ids)} UniProt mappings in {time.time()-t0:.1f}s")

    with open("results/casf2016/uniprot_map.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["pdbid","uniprot_ids"]); w.writerows(results)

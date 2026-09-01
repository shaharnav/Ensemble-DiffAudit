"""
Download CASF-2016 core set structures from RCSB (no registration).
PDB ID / pKd list sourced from GIGN's test2016.csv (github.com/guaguabujianle/GIGN),
which reproduces the CASF-2016 core set composition. This is a reproduction of
CASF-2016, not the official PDBbind-processed structures.
"""
import csv, os, time, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

OUT = "results/casf2016/structures"
os.makedirs(OUT, exist_ok=True)

with open("results/casf2016/casf2016_core_pdbid_pkd.csv") as f:
    ids = [r["pdbid"] for r in csv.DictReader(f)]

def fetch(pdbid):
    dest = os.path.join(OUT, f"{pdbid}.pdb")
    if os.path.exists(dest) and os.path.getsize(dest) > 500:
        return pdbid, "cached"
    url = f"https://files.rcsb.org/download/{pdbid.upper()}.pdb"
    r = subprocess.run(["curl","-sL","--max-time","20","-o",dest,url], capture_output=True)
    if r.returncode == 0 and os.path.exists(dest) and os.path.getsize(dest) > 500:
        return pdbid, "ok"
    return pdbid, f"FAIL: rc={r.returncode}"

t0 = time.time()
results = []
with ThreadPoolExecutor(max_workers=12) as pool:
    futs = {pool.submit(fetch, i): i for i in ids}
    for fut in as_completed(futs):
        pdbid, status = fut.result()
        results.append((pdbid, status))
        if not status.startswith(("ok","cached")):
            print(pdbid, status)

elapsed = time.time()-t0
ok = sum(1 for _,s in results if s in ("ok","cached"))
print(f"Downloaded {ok}/{len(ids)} in {elapsed:.1f}s")
with open("results/casf2016/download_log.csv","w",newline="") as f:
    w = csv.writer(f); w.writerow(["pdbid","status"]); w.writerows(results)

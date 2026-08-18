"""Featurize all Stage-0-successful complexes with the full A-F feature set."""
import csv, time
from multiprocessing import Pool
from casf_pipeline.features import extract_all_features

def process(pdbid):
    try:
        f = extract_all_features(
            f"results/casf2016/pdbqt/{pdbid}_receptor.pdbqt",
            f"results/casf2016/pdbqt/{pdbid}_ligand.pdbqt",
            f"results/casf2016/ligands/{pdbid}_ligand.sdf",
        )
        if f is None:
            return pdbid, None
        return pdbid, f
    except Exception as e:
        return pdbid, None

if __name__ == "__main__":
    with open("results/casf2016/stage0_scores.csv") as f:
        rows = {r["pdbid"]: r for r in csv.DictReader(f) if r["status"] == "ok"}
    ok_ids = list(rows.keys())

    t0 = time.time()
    with Pool(processes=8) as pool:
        results = pool.map(process, ok_ids)
    elapsed = time.time() - t0

    ok = {pdbid: f for pdbid, f in results if f is not None}
    print(f"Featurized {len(ok)}/{len(ok_ids)} in {elapsed:.1f}s")

    all_cols = sorted(set().union(*[f.keys() for f in ok.values()]) - {"hac"})
    with open("results/casf2016/stage1_features.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pdbid","pKd","vina","vinardo","hac"] + all_cols)
        w.writeheader()
        for pdbid, feat in ok.items():
            row = {"pdbid": pdbid, "pKd": rows[pdbid]["pKd"], "vina": rows[pdbid]["vina"],
                   "vinardo": rows[pdbid]["vinardo"]}
            row.update(feat)  # feat's own 'hac' becomes the single hac column
            w.writerow(row)
    print("Written results/casf2016/stage1_features.csv")

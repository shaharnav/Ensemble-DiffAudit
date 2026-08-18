"""
Phase 7 Step 3 (partial) -- Vinardo rescoring of all existing poses.
No new docking: reads existing *_out.pdbqt (best pose) + *_receptor.pdbqt + *_params.json
already on disk from Phase 4 (ensemble_audit), Control A, Control B, and Phase 6.
Uses the local vina 1.2.7 binary's native --scoring vinardo --score_only mode.
"""
import os, glob, json, subprocess, csv, re

VINA = os.path.abspath("./bin/vina_1.2.7_mac_aarch64")
TMP = "/tmp/p7_vinardo_pose.pdbqt"

DIRS = {
    "phase4_ensemble": "results/ensemble_audit",
    "control_a": "results/phase5_control_a",
    "control_b": "results/phase5_control_b",
    "phase6_actives": "results/phase6_docking",
}

def extract_best_pose(out_pdbqt, dest):
    with open(out_pdbqt) as f:
        lines = f.readlines()
    keep, in_model1, seen = [], False, False
    for ln in lines:
        if ln.startswith("MODEL"):
            if not seen:
                in_model1 = True
                seen = True
            else:
                break
            continue
        if ln.startswith("ENDMDL"):
            if in_model1:
                break
            continue
        if in_model1:
            keep.append(ln)
    if not keep:
        # no MODEL tags at all -- already a single-pose file
        keep = [l for l in lines if not l.startswith(("MODEL","ENDMDL"))]
    with open(dest, "w") as f:
        f.writelines(keep)

def vinardo_score(receptor_pdbqt, ligand_pdbqt, center, size):
    cmd = [VINA, "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--scoring", "vinardo", "--score_only"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.\d+)", r.stdout)
    if m:
        return float(m.group(1))
    return None

rows = []
for label, d in DIRS.items():
    out_files = sorted(glob.glob(os.path.join(d, "*_out.pdbqt")))
    print(f"{label}: {len(out_files)} pose files")
    for outf in out_files:
        base = outf[:-len("_out.pdbqt")]
        recf = base + "_receptor.pdbqt"
        paramsf = base + "_params.json"
        if not (os.path.exists(recf) and os.path.exists(paramsf)):
            rows.append({"label": label, "job": os.path.basename(base), "vinardo": None, "error": "missing_receptor_or_params"})
            continue
        with open(paramsf) as f:
            p = json.load(f)
        extract_best_pose(outf, TMP)
        score = vinardo_score(recf, TMP, p["center"], p["size"])
        rows.append({"label": label, "job": os.path.basename(base), "vinardo": score,
                     "error": None if score is not None else "vina_score_only_failed"})

fails = sum(1 for r in rows if r["vinardo"] is None)
print(f"Done. {len(rows)-fails}/{len(rows)} succeeded, {fails} failed.")

with open("results/phase7_vinardo_scores.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["label","job","vinardo","error"])
    w.writeheader()
    w.writerows(rows)
print("Written results/phase7_vinardo_scores.csv")

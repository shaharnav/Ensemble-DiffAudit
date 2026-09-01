"""
Step 1d: prepare all receptors (the 4 validated Zn-regenerated conformers,
the 1EZM crystal, and the zinc-stripped 1EZM control) through the identical
meeko pipeline, confirm metal presence/atom typing in each output pdbqt, and
write the receptor provenance table.
"""
import csv, os, sys
sys.path.insert(0, ".")
from casf_pipeline.prep import prepare_receptor_pdbqt

RAW_DIR = "results/lasb_ensemble_rmsd/receptors_raw"

ACCEPTED_CONFORMERS = ["beta0.0", "beta3.2"]
REJECTED_CONFORMERS = ["beta0.8", "beta1.6", "beta2.4", "beta4.0"]
REJECT_REASONS = {
    "beta0.8": "Zn/Ca validation FAIL -- Glu164-Zn 1.40 A (too short)",
    "beta1.6": "Zn/Ca validation FAIL -- His144-Zn 1.79 A (0.01 A under floor)",
    "beta2.4": "structural clash unrelated to Zn site -- Arg208(NH2)-His224(CE1) 1.73 A",
    "beta4.0": "structural clash unrelated to Zn site -- Arg274(NH2)-Asn278(OD1) 1.56 A",
}

receptors = [
    {"receptor_id": "1EZM_crystal",
     "source_pdb": "results/lasb_payload/ensemble_receptors_aligned/1EZM_baseline_crystal.pdb",
     "role": "Condition B (apo crystal)"},
    {"receptor_id": "1EZM_zinc_stripped",
     "source_pdb": f"{RAW_DIR}/1EZM_zinc_stripped.pdb",
     "role": "zinc-strip control (isolates pure volume effect)"},
]
for beta in ACCEPTED_CONFORMERS:
    receptors.append({"receptor_id": beta,
                       "source_pdb": f"{RAW_DIR}/{beta}_zn_regenerated.pdb",
                       "role": "Condition A (ensemble)"})

rows = []
for r in receptors:
    out_pdbqt = f"{RAW_DIR}/{r['receptor_id']}_prepped.pdbqt"
    ok = prepare_receptor_pdbqt(r["source_pdb"], out_pdbqt)
    n_atoms = zn_present = ca_present = 0
    if ok and os.path.exists(out_pdbqt):
        with open(out_pdbqt) as f:
            lines = [l for l in f if l.startswith(("ATOM", "HETATM"))]
        n_atoms = len(lines)
        zn_present = any(l.rstrip().endswith("Zn") for l in lines)
        ca_present = any(l.rstrip().endswith("Ca") for l in lines)
    rows.append({**r, "prep_ok": ok, "pdbqt": out_pdbqt if ok else "",
                 "n_atoms": n_atoms, "zn_present": zn_present, "ca_present": ca_present})
    print(f"{r['receptor_id']}: prep {'ok' if ok else 'FAILED'}, "
          f"n_atoms={n_atoms}, Zn={zn_present}, Ca={ca_present}")

for beta in REJECTED_CONFORMERS:
    rows.append({"receptor_id": beta, "source_pdb": f"{RAW_DIR}/{beta}_zn_regenerated.pdb",
                 "role": f"EXCLUDED -- {REJECT_REASONS[beta]}",
                 "prep_ok": "", "pdbqt": "", "n_atoms": "", "zn_present": "", "ca_present": ""})

with open("results/lasb_ensemble_rmsd/receptor_provenance_table.csv", "w", newline="") as f:
    fieldnames = ["receptor_id", "role", "source_pdb", "prep_ok", "pdbqt", "n_atoms", "zn_present", "ca_present"]
    w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
    w.writeheader()
    w.writerows(rows)
print("\nWritten results/lasb_ensemble_rmsd/receptor_provenance_table.csv")

"""
Fascin Step 3 (continued): prep the 4 conformers that survived the geometry
gates (beta4.0, beta4.8, beta5.6, beta6.0) through meeko into pdbqt receptors,
plus the apo crystal (3LLP) as Condition B, matching the LasB pipeline.
No metal cofactor exists for fascin, so no transplant step is needed here.
"""
import gemmi
import subprocess
import os

CONF_DIR = "results/fascin_ensemble_rmsd/conformers_raw"
OUT_DIR = "results/fascin_ensemble_rmsd/receptors_raw"
os.makedirs(OUT_DIR, exist_ok=True)

SURVIVORS = {
    "beta4.0": "4.0",
    "beta4.8": "4.800000000000001",
    "beta5.6": "5.6000000000000005",
    "beta6.0": "6.0",
}

results = []
for label, filekey in SURVIVORS.items():
    cif_path = f"{CONF_DIR}/beta{filekey}.cif"
    pdb_path = f"{OUT_DIR}/{label}.pdb"
    pdbqt_path = f"{OUT_DIR}/{label}.pdbqt"

    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    st.write_pdb(pdb_path)

    r = subprocess.run(
        [os.path.abspath("venv/bin/python3"), os.path.abspath("venv/bin/mk_prepare_receptor.py"),
         "--read_pdb", os.path.abspath(pdb_path),
         "-o", os.path.abspath(pdbqt_path.replace(".pdbqt", "")), "-p"],
        capture_output=True, text=True
    )
    ok = os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 0
    results.append({"label": label, "ok": ok, "stderr_tail": r.stderr[-300:] if not ok else ""})
    print(f"{label}: prep {'OK' if ok else 'FAILED'}" + (f"\n  {r.stderr[-300:]}" if not ok else ""))

# also prep the apo crystal (Condition B) -- strip crystallization additives
# (Br-, HEPES, glycerol, K+, sulfate, water) first; none are biologically
# relevant here, but meeko's bond-perception chokes on K+ (no covalent radius
# implemented). Correction to Step 1a's summary: 3LLP is NOT literally free of
# all HETATM records -- it has these standard buffer components, just no real
# (non-additive) ligand, which is what "apo" meant there.
apo_src = "results/target_screen/structures/3LLP.pdb"
apo_st = gemmi.read_structure(apo_src)
apo_st.setup_entities()
apo_st.remove_ligands_and_waters()
apo_pdb = f"{OUT_DIR}/3LLP_stripped.pdb"
apo_st.write_pdb(apo_pdb)
apo_pdbqt = f"{OUT_DIR}/3LLP_apo"
r = subprocess.run(
    [os.path.abspath("venv/bin/python3"), os.path.abspath("venv/bin/mk_prepare_receptor.py"),
     "--read_pdb", os.path.abspath(apo_pdb),
     "-o", os.path.abspath(apo_pdbqt), "-p"],
    capture_output=True, text=True
)
ok = os.path.exists(apo_pdbqt + ".pdbqt")
print(f"3LLP (apo, Condition B): prep {'OK' if ok else 'FAILED'}")
if not ok:
    print(r.stderr[-500:])

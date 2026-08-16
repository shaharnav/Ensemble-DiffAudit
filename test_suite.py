import logging
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from fetcher import fetch_alphafold_structure
from docking_engine import run_docking
from ensemble_auditor import aggregate_candidate_affinities

logging.basicConfig(level=logging.INFO)

# Test: crystal beats every conformer -> delta must be negative, not clamped to zero.
print("--- TESTING aggregate_candidate_affinities (crystal beats all conformers) ---")
synthetic_affinities = {
    "3PTB_baseline_crystal.pdb": -6.684,
    "conformix_var_0.pdb": -6.247,
    "conformix_var_1.pdb": -6.100,
    "conformix_var_2.pdb": -5.954,
}
metrics = aggregate_candidate_affinities(synthetic_affinities)
assert metrics["crystal_affinity"] == -6.684
assert metrics["ensemble_best_affinity"] == -6.247
assert metrics["delta_ensemble_vs_crystal"] < 0, (
    f"Expected negative delta (crystal outperformed ensemble), got "
    f"{metrics['delta_ensemble_vs_crystal']}"
)
print(f"Success: delta_ensemble_vs_crystal = {metrics['delta_ensemble_vs_crystal']:.3f} (negative, as expected)")

# Test Metal Priority with P00918 (Carbonic Anhydrase)
print("--- TESTING P00918 ---")
pdb_path1 = fetch_alphafold_structure("P00918")
res1 = run_docking(pdb_path1, "CC", job_name="test_p00918", exhaustiveness=1)
if res1:
    print("Success P00918")

# Test Cavity Finding with P29274 (Caffeine / A2A)
print("\n--- TESTING P29274 ---")
pdb_path2 = fetch_alphafold_structure("P29274")
res2 = run_docking(pdb_path2, "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", job_name="test_p29274", exhaustiveness=1)
if res2:
    print("Success P29274")

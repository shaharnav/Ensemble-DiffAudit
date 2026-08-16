import logging
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from fetcher import fetch_alphafold_structure
from docking_engine import run_docking
from ensemble_auditor import aggregate_candidate_affinities
from rigid_control import summarize_rigid_scores
from chem_metrics import structural_alerts
from rdkit import Chem

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

# Test: if the rigid-control seeds land exactly on the ensemble's own conformer scores
# (i.e. zero real conformational change), noise_corrected_delta should be ~0.
print("--- TESTING summarize_rigid_scores (zero conformational change -> zero delta) ---")
conformer_scores = [-6.247, -6.135, -6.173, -5.954, -6.238, -6.156]
ensemble_best_affinity = min(conformer_scores)
rigid_summary = summarize_rigid_scores(conformer_scores, ensemble_best_affinity)
assert abs(rigid_summary["noise_corrected_delta"]) < 1e-9, (
    f"Expected noise_corrected_delta ~0 when rigid scores == conformer scores, "
    f"got {rigid_summary['noise_corrected_delta']}"
)
print(f"Success: noise_corrected_delta = {rigid_summary['noise_corrected_delta']:.6f} (~0, as expected)")

# Regression fixtures: known structural-alert matches from the current N=7 candidate set.
print("--- TESTING structural_alerts (known fixtures) ---")
known_fixtures = {
    "Cmpd-0006": ("CC[C@H]1SCC[C@H](F)C(=O)SC(C)(C)CC=C=C1[C@H](O)CCO", {"allene"}),
    "Cmpd-0000": ("O=CN1C(=O)[C@H](CP(=O)(O)O)S[C@H]1CCC[C@H](C[C@H]1CC=CO1)[C@H]1CO1", {"epoxide", "phosphonic_acid"}),
    "Cmpd-0003": ("C[C@@H](F)[C@@H](F)[C@@H]1O[C@@H]([C@@H]2C[C@H](C)[C@@H]3O[C@@H]3[C@@H]2S)C(F)(F)OC[C@]1(C)O", {"epoxide", "free_thiol"}),
    "Cmpd-0007": ("Cc1nc(N)sc1C[C@H](C)CNC(=O)[C@H]1C=CC=CC1=O", {"michael_acceptor"}),
    "Cmpd-0008": ("C[C@H](O)[C@]1(C(=O)O)CC2=C(C[C@H](Cl)CC2)[C@@H]1N", {"alkyl_halide"}),
}
for cmpd_id, (smiles, expected_alerts) in known_fixtures.items():
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"{cmpd_id}: failed to parse SMILES"
    found = set(structural_alerts(mol))
    assert found == expected_alerts, f"{cmpd_id}: expected {expected_alerts}, got {found}"
print(f"Success: all {len(known_fixtures)} known structural-alert fixtures matched")

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

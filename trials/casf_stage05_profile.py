"""Stage 0.5 -- profile full feature set on 50 complexes, broken down by group,
extrapolate to the full available set. Time box: 15 min."""
import csv, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from casf_pipeline.pdbqt_atoms import parse_pdbqt_atoms, atoms_to_arrays, find_donor_heavy_atoms
from casf_pipeline import features as F

with open("results/casf2016/stage0_scores.csv") as f:
    ok_ids = [r["pdbid"] for r in csv.DictReader(f) if r["status"] == "ok"][:50]

print(f"Profiling {len(ok_ids)} complexes...")

group_times = {"parse": 0.0, "A_vina_terms": 0.0, "C_contacts": 0.0, "D_E_hbond": 0.0, "B_sasa": 0.0, "F_descriptors": 0.0}
n_ok = 0
t_total0 = time.time()

for pdbid in ok_ids:
    rec_pdbqt = f"results/casf2016/pdbqt/{pdbid}_receptor.pdbqt"
    lig_pdbqt = f"results/casf2016/pdbqt/{pdbid}_ligand.pdbqt"
    lig_sdf = f"results/casf2016/ligands/{pdbid}_ligand.sdf"
    try:
        t0 = time.time()
        rec_atoms = parse_pdbqt_atoms(rec_pdbqt)
        lig_atoms = parse_pdbqt_atoms(lig_pdbqt)
        rec_coords, rec_radii, rec_types, rec_hp, rec_acc, rec_don_heavy = atoms_to_arrays(rec_atoms)
        lig_coords, lig_radii, lig_types, lig_hp, lig_acc, lig_don_heavy = atoms_to_arrays(lig_atoms)
        rec_don = find_donor_heavy_atoms(rec_atoms, rec_coords) & rec_don_heavy
        lig_don = find_donor_heavy_atoms(lig_atoms, lig_coords) & lig_don_heavy
        group_times["parse"] += time.time() - t0

        t0 = time.time()
        F.vina_terms(rec_coords, rec_radii, rec_hp, rec_acc, rec_don, lig_coords, lig_radii, lig_hp, lig_acc, lig_don)
        group_times["A_vina_terms"] += time.time() - t0

        t0 = time.time()
        F.contact_counts(rec_coords, rec_types, lig_coords, lig_types)
        group_times["C_contacts"] += time.time() - t0

        t0 = time.time()
        F.hbond_and_satisfaction(rec_coords, rec_acc, rec_don, lig_coords, lig_acc, lig_don, lig_types)
        group_times["D_E_hbond"] += time.time() - t0

        t0 = time.time()
        F.buried_sasa(lig_pdbqt, rec_coords=rec_coords, rec_atoms=rec_atoms, lig_coords=lig_coords, lig_atoms=lig_atoms)
        group_times["B_sasa"] += time.time() - t0

        t0 = time.time()
        F.ligand_descriptors(lig_sdf)
        group_times["F_descriptors"] += time.time() - t0

        n_ok += 1
    except Exception as e:
        print(f"  {pdbid} FAILED: {e}")

elapsed = time.time() - t_total0
print(f"\n{n_ok}/{len(ok_ids)} succeeded in {elapsed:.2f}s total ({elapsed/n_ok*1000:.1f}ms/complex)")
print("\nBreakdown by group (ms/complex):")
for g, t in group_times.items():
    print(f"  {g:20s} {t/n_ok*1000:8.2f} ms/complex  ({t/elapsed*100:5.1f}% of total)")

full_n = 244
projected = elapsed / n_ok * full_n
print(f"\nProjected time for full {full_n}-complex set: {projected:.1f}s ({projected/60:.2f} min)")
print("GATE: 2-hour box.", "PASS - well under budget." if projected < 7200 else "FAIL - exceeds budget, stop for review.")

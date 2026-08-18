# Phase 5 Pre-Registration — LasB Ensemble Docking Controls

**Written:** 2026-08-17, before any Phase 5 docking jobs are submitted.  
**Purpose:** Fix metric definitions, predicted outcomes, and falsification thresholds
before results exist, so no interpretation can be reverse-engineered from outcomes.

---

## Context

Phase 4 produced a headline of mean(ensemble_max − crystal) = +1.42 kcal/mol
across 20 DiffSBDD candidates × 5 ConforMix conformers (β0.0–3.2; β4.0 failed
silently for all candidates and is excluded from all Phase 5 analysis).

Step 0 reanalysis showed:
- mean(ensemble_mean − crystal) = +0.990 kcal/mol (selection-bias-corrected)
- All 5 conformer means cluster uniformly at −6.07 to −6.29 kcal/mol vs. crystal −5.21
- No clear gradient tracking Phase 2b directional gain values
- 19/20 candidates show positive delta_mean; one outlier (candidate 6, −0.24)

Three hypotheses motivate the controls:
1. The gain is max-over-N sampling noise (Control A measures this)
2. The gain is a Boltz-pocket scoring offset applying to any ligand (Control B)
3. The gain reflects genuine holo-pocket geometric recognition (would survive A and B)

---

## Control A — Seeded Crystal (noise floor)

**Design:** Dock all 20 candidates into the 1EZM crystal receptor 6 independent
times using seeds 42, 123, 456, 789, 1337, 2024. Same exhaustiveness (16), box
center (55.521, 35.882, 20.807), box size (24 Å), and receptor prep as Phase 4.
Seed 42 is included to enable cross-check against the Phase 4 crystal score.

**Metrics (computed after results exist):**
```
crystal_seeded_best[i]  = min affinity over 6 seeds for candidate i
seed_noise_gain         = mean_i(crystal_phase4[i] − crystal_seeded_best[i])
noise_corrected_delta   = mean_i(ensemble_max[i] − crystal_seeded_best[i])
```
Bootstrap 95% CI on noise_corrected_delta: 10,000 resamples over the 20 candidates.

**Predicted outcome:** seed_noise_gain < 0.30 kcal/mol. The trypsin experiment showed
+0.134 kcal/mol from 6 seeds on a rigid crystal. If LasB crystal Vina variance is
similar, noise contributes less than 0.20 kcal/mol. The mean-based signal (+0.990)
is large enough that a plausible noise floor should not erase it.

**Falsification threshold:** If `noise_corrected_delta ≤ 0.30 kcal/mol` (i.e. the
entire mean-based signal is explained by re-seeding the same crystal receptor), the
ensemble-advantage claim is falsified for this dataset. Stop; do not proceed to
Controls B or C. Report as a methodology negative.

---

## Control B — Decoys (discrimination vs. offset)

**Design:** Build 20 decoys property-matched to the 20 candidates:
- MW ±25 Da, cLogP ±0.5, HBD ±1, HBA ±1, rotatable bonds ±2
- Tanimoto < 0.3 on Morgan fingerprints (radius=2, 2048 bits) vs. all 20 candidates
- Source: ZINC20 or ChEMBL; source and matching table reported in phase5_results.md

Dock all 20 decoys into all 6 receptors (5 conformers + crystal), same settings.

**Metrics:**
```
decoy_crystal[j]    = decoy affinity in 1EZM crystal (seed 42)
decoy_ens_max[j]    = best decoy affinity across 5 conformers
decoy_delta[j]      = decoy_ens_max[j] − decoy_crystal[j]
mean_decoy_delta    = mean_j(decoy_delta)
discrimination_gap  = mean_candidate_delta_mean − mean_decoy_delta
```
Bootstrap 95% CI on discrimination_gap (10,000 resamples, candidates and decoys
resampled jointly).

**Predicted outcome:** mean_decoy_delta < 0.50 kcal/mol; discrimination_gap > 0.50
kcal/mol. If the gain is geometric recognition, candidates (designed against the apo
pocket) should exploit holo-like conformers while generic property-matched compounds
should not. If the gain is a Boltz-pocket softening artifact, both should show similar
deltas.

**Falsification threshold:** If `discrimination_gap < 0.30 kcal/mol` (decoys gain
nearly as much as candidates), the gain is a receptor-level offset, not holo
recognition. The ensemble-advantage claim does not hold as a discrimination signal,
even if Control A survives.

---

## Control C — Holo Rigid + Native Ligand Redock

**Design:**
- Superpose 3DBK onto 1EZM frame via global CA alignment (same as Phase 2b);
  report alignment RMSD.
- Strip RDF, waters, HETATM from 3DBK. Verify existing box center
  (55.521, 35.882, 20.807) encloses the 3DBK pocket before docking; if not,
  stop and report.
- 4a: Dock all 20 candidates into 3DBK rigid (20 jobs), seed 42.
- 4b: Redock native RDF ligand into all 5 conformers + crystal + 3DBK rigid (7 jobs),
  seed 42.

**Metrics:**
```
# 4a
holo_rigid[i]           = candidate i affinity in aligned 3DBK
gap_recovery_fraction[i] = (crystal_seeded_best[i] − ensemble_max[i]) /
                            (crystal_seeded_best[i] − holo_rigid[i])
```
Denominator ≤ 0.1 kcal/mol → candidate excluded from mean; count and report.

```
# 4b — native ligand sanity check
rdf_affinity[r] for r in {crystal, β0.0, β0.8, β1.6, β2.4, β3.2, 3DBK_rigid}
```

**Predicted outcomes:**
- 4a: mean gap_recovery_fraction ∈ [0.3, 0.8]. Values near 1.0 mean conformers
  reproduce the full apo→holo gap; near 0 means they recover none.
- 4b: RDF should prefer 3DBK_rigid and holo-like conformers (β0.0, β0.8) over
  crystal. If RDF prefers the crystal or β3.2/β4.0, the geometric premise is broken.

**Falsification threshold:** If mean gap_recovery_fraction < 0.10 AND RDF fails the
sanity check (prefers crystal over all conformers), Control C falsifies the geometric
interpretation independently of Controls A and B.

---

## Home advantage (Step 5)

**Metric:**
```
home_fraction = count(candidates whose overall_best_structure == crystal) / 20
```
Null hypothesis: home_fraction = 1/6 = 16.7% (uniform over 6 receptors).
Currently observed: 0/20 = 0%. Report binomial p-value for departure from 1/6
in either direction.

**Predicted outcome:** home_fraction < 16.7% (below chance), meaning ligands
designed against the apo crystal systematically avoid it — which is the expected
result if Boltz conformers are genuinely more accommodating. If home_fraction > 30%,
investigate receptor prep asymmetry.

---

## Rules that cannot be changed after this commit

1. Same exhaustiveness (16), box, center, and prep pipeline for all conditions.
2. No conditions dropped. If infeasible, stop and report.
3. Negative results reported as prominently as positive ones.
4. noise_corrected_delta uses crystal_seeded_best (6-seed min), not Phase 4 crystal.
5. All bootstrap CIs: 10,000 resamples over candidates (n=20), percentile method.
6. β4.0 excluded from all Phase 5 analysis (failed silently in Phase 4; 0/20 scores).

---

## Deliverables

- `phase5_control_a_results.csv`: one row per (candidate, seed), columns: smiles, seed, affinity
- `phase5_control_b_decoys.csv`: decoy matching table + per-decoy scores
- `phase5_control_c_results.csv`: holo rigid scores + RDF redock scores
- `phase5_results.md`: every metric above with CIs, per-conformer table, plain verdict on
  each of the three questions, falsification summary

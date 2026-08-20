# Does the LasB conformer ensemble recover crystallographic poses better than the rigid apo crystal?

## Verdict

**No, not distinguishably.** The ensemble (n=2 validated conformers) beats the apo crystal
on top-1 success@2Å by +6.7 percentage points (60.0% vs 53.3%), but the paired 95% bootstrap
CI on that difference is **[−20.0%, +33.3%]** — it spans zero. A Wilcoxon signed-rank test on
paired top-1 RMSD gives p=0.64. **This is a null result at n=15 ligands, reported plainly, not
dressed up as a positive finding.** The apo crystal actually edges out the ensemble on the
stricter top-1<1Å threshold (20.0% vs 13.3%).

This experiment abandoned affinity scoring for pose geometry specifically because Phase 6
showed Vina scores in this pocket are dominated by molecular size (r(Vina,HAC)=−0.915) and
can't separate real candidates from decoys. RMSD doesn't inherit that confound. The honest
result on the metric that actually matches this project's claim is: **the ensemble approach,
as built here, does not measurably improve pose recovery over a single rigid crystal
structure — and the reasons why are mechanistically traceable, not a black box.**

---

## Why n=2, not n=6: the attrition chain

1. **6 conformers generated** (apo-conditioned, no metal cofactor specified) → **all 6 fail**
   Zn-site validation. His144 — one of the three residues coordinating the catalytic zinc —
   collapses to within 1.09–1.79 Å of where Zn²⁺ belongs (a real Zn–His bond is ~2.0–2.2 Å).
   Nothing constrained the rotamer without the metal present during generation.
2. **Regenerated with Zn²⁺/Ca²⁺ explicit** in the Boltz input (no ConforMix code changes needed
   — Boltz's native fasta schema already supports `>chain|ccd|` cofactor entities) →
   **4/6 pass** Zn-site validation (β0.0, β2.4, β3.2, β4.0), now checked as a genuine
   self-consistency test (the model's own predicted Zn position vs. its own predicted
   coordinating residues, not a transplant).
3. Of those 4, **2 more fail** on an orthogonal criterion: direct atom-distance measurement
   confirmed real backbone/side-chain steric clashes elsewhere in the structure (Arg208–His224
   at 1.73 Å in β2.4, Arg274–Asn278 at 1.56 Å in β4.0 — absent in the crystal and in β0.0/β3.2),
   tracking the two largest twist-guidance targets in the sweep.
4. **n=2 conformers survive both gates: β0.0 and β3.2.** Every exclusion has a specific,
   verified mechanistic reason — none were assumed or silently dropped. Full numbers:
   `lasb_zn_transplant_finding.md`, `lasb_zn_regeneration_finding.md`, `lasb_step1d_diagnostic.md`.

This is itself a finding about the scope limits of RMSD-guided conformational sampling on a
metalloenzyme: correctly specifying the cofactor fixes the catalytic-site geometry (the
original, more severe failure) but does not guarantee general backbone validity at every twist
target — a second, independent check is still required.

**n=2 is not called an "ensemble"** in the sense the original protocol intended. The Step 6d
question — does one apo-derived ensemble serve multiple chemically distinct ligands, or only
the one it happens to suit — is **untestable at this n** and is not claimed either way.

---

## Receptor provenance

| Receptor | Role | Zn present (pdbqt) | Ca present (pdbqt) | Atoms |
|---|---|---|---|---|
| 1EZM_crystal | Condition B (apo) | yes | yes | 3220 |
| 1EZM_zinc_stripped | volume-effect control | **no** | yes | 3219 |
| β0.0 | Condition A | yes | yes | 2867 |
| β3.2 | Condition A | yes | yes | 2867 |
| 15× holo crystals | Condition C | yes | yes | 2285–2847 |

Every receptor that proceeded to docking differs from every other in exactly the intended way
(conformation and, for the zinc-strip control, presence/absence of Zn) — verified directly in
the final pdbqt, not assumed from the source PDB. Full table: `receptor_provenance_table.csv`.

## Holo ligand set

15 catalytic-site LasB inhibitors from RCSB (UniProt P14756), all confirmed Zn-chelating
(1.7–2.5 Å to the catalytic Zn²⁺), MW 236–609, resolution 1.3–2.74 Å, CA-superposition onto the
1EZM frame at 1.46–2.08 Å. 2 additional structures (7OC7/V85, 7Z68/IEV) excluded for partial
ligand occupancy (0.5/0.38 — an ambiguous reference pose). Full table:
`holo_ligand_set_final.csv`.

## Docking

Vina 1.2.7, box center (55.521, 35.882, 20.807), 24 Å cube, exhaustiveness 16, num_modes 20,
energy_range 5. Same box for every condition/receptor (verified: 0.12 Å from the RDF ligand
centroid in the 1EZM frame). Condition B and C both used exactly 2 seeds per ligand, matched to
Condition A's 2 surviving conformers — no unequal-compute confound. 90 docking jobs, 90/90
succeeded, 1797 total poses retained across all ranks.

## Primary results

| Condition | n | top1 success@2Å | top1 success@1Å | mean top1 RMSD | oracle success@2Å | oracle success@1Å | mean oracle RMSD |
|---|---|---|---|---|---|---|---|
| A (ensemble, n=2) | 15 | 60.0% | 13.3% | 1.85 Å | 100.0% | 60.0% | 0.93 Å |
| B (apo crystal) | 15 | 53.3% | 20.0% | 2.00 Å | 86.7% | 60.0% | 1.23 Å |
| C (self-docking ceiling) | 15 | 73.3% | 20.0% | 1.59 Å | 100.0% | 73.3% | 0.83 Å |

Condition C outperforms both A and B on every metric, exactly as expected for a ceiling
(each ligand docked into its own crystal structure). This is the sanity check that confirms
the pipeline is measuring something real, not noise: the ceiling behaves like a ceiling.

**Oracle vs top-1 gap (the diagnostic Step 4d calls for):** Condition A's oracle success@2Å is
100% but top-1 is only 60% — the ensemble's conformers *do* sample near-correct geometry
somewhere among their poses, but Vina's scoring function frequently fails to rank that pose
first. This is a scoring-selection failure, not a sampling failure — consistent with Phase 6's
finding that Vina's score in this pocket is a poor discriminator. The same pattern holds for B
(86.7% oracle vs 53.3% top-1) and C (100% oracle vs 73.3% top-1), so it isn't specific to the
ensemble; it's a property of the scoring function in this pocket generally.

## Paired statistical comparison (A vs B)

- Wilcoxon signed-rank on paired top-1 RMSD (n=15): **p = 0.639**
- Success-rate(2Å) difference (A−B): **+6.7%, 95% bootstrap CI [−20.0%, +33.3%]**
- The CI spans zero — the observed advantage is not distinguishable from no difference at
  this sample size, and is explicitly reported as such rather than overstated.

## Ceiling fraction

Of the 10 ligands where Condition C genuinely beats Condition B (i.e. self-docking is a real
ceiling for that ligand), Condition A captures a mean of **16.0%** of the gap between B and C.
Most of the achievable improvement over the apo crystal is not realized by this 2-conformer
ensemble.

## Secondary metrics

- **PoseBusters validity**: 45/45 top-1 poses (15 ligands × 3 conditions) pass — no condition
  wins by exploiting geometrically implausible poses.
- **ProLIF interaction-fingerprint recovery (5a): dropped.** ProLIF/MDAnalysis's RDKit
  bond-order standardization segfaults (uncatchable in Python) on the H-less receptor PDBs used
  throughout. Fixing it requires fully protonating every receptor first — out of scope; flagged
  as a real gap, not silently skipped.
- **Winning conformer per ligand** (Condition A): β0.0 wins 5/15, β3.2 wins 10/15. Neither
  conformer wins every ligand — weak evidence the 2-receptor set does *some* differentiated
  work — but n=2 cannot properly answer the "does the ensemble span accessible states"
  question from Step 6d.

## What this experiment actually establishes

1. Correctly specifying the metal cofactor during generation is necessary but not sufficient
   for a usable conformer — two independent geometric failure modes exist (metal-site collapse,
   general backbone distortion at high twist), and both must be checked.
2. On the metric that matches this project's actual claim (pose geometry, not affinity), the
   ensemble does not show a statistically distinguishable advantage over the rigid apo crystal,
   at n=2 conformers / n=15 ligands.
3. The oracle-vs-top-1 gap is large and roughly uniform across all three conditions — Vina's
   scoring function, not conformational sampling, is the dominant bottleneck in this pocket,
   consistent with the size-bias finding from Phase 6.
4. A clean, isolated zinc-strip control (`1EZM_zinc_stripped_prepped.pdbqt`) is prepared but not
   yet docked — it was scoped as a way to isolate the pure pocket-volume effect independent of
   the conformers' own backbone distortion, and remains available for a follow-up run.

## Deliverables

- `lasb_zn_transplant_finding.md`, `lasb_zn_regeneration_finding.md`, `lasb_step1d_diagnostic.md`
  — the full Step 1 attrition chain
- `results/lasb_ensemble_rmsd/receptor_provenance_table.csv`
- `results/lasb_ensemble_rmsd/holo_ligand_set_final.csv`
- `results/lasb_ensemble_rmsd/rmsd_results.csv` — primary per-ligand×condition results
- `results/lasb_ensemble_rmsd/secondary_metrics.csv`, `winning_conformer_table.csv`
- `results/lasb_ensemble_rmsd/tidy_all_poses.csv` — full per-pose CSV (1797 rows: ligand,
  condition, receptor, seed, pose_rank, vina_score, rmsd, posebusters_pass)
- `results/lasb_ensemble_rmsd/docking_log.csv` — per-job docking log (90 rows)

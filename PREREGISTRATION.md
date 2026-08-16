# Pre-registration — Experiment 2 Pilot: Apo/Holo Cross-Docking

This document is committed before any ConforMix generation or docking run for this
experiment. Any number reported for this pilot must postdate the commit hash of this file.

## 1. Screening criteria (fixed before the Phase 1b screen ran)

A candidate apo/holo pair qualifies only if **all** hold, computed by `target_screen.py`
against `apo_holo_geometry.compute_apo_holo_geometry`:

1. **No cofactor or construct difference.** Same organism and construct, no mutations, no
   bound species other than the drug-like ligand under study. A multi-component annotated
   ligand field (e.g. "NDP,TOP") counts as a cofactor-difference risk on its own, since it's
   ambiguous which component is the drug and which might be a second bound species.
   Differing cofactor occupancy is **disqualifying, not documentable**.
2. **Backbone-driven motion.** CA-to-all-atom pocket RMSD ratio >= 0.5. Sidechain-dominated
   pairs (low ratio) are out of scope for a backbone sampler.
3. **Displacement within the sampler's reach.** Pocket CA RMSD in **1.0-2.5 A** -- clear of
   the ~0.2 A noise floor seen with trypsin, inside the `--twist-target-stop` 2.0 A shell.
4. **Drug-like holo ligand** (MW 150-600; not a metal ion or crystallization additive).
5. **Tractable structure.** Single chain or clean single domain, manually confirmed.

## 2. Rejection record

Three candidates were rejected before generating or docking anything, for structural reasons
recorded in `targets.yaml`'s `target_swap_history`:

- **DHFR (2W9T/2W9S)** -- rejected. Holo carries NADPH in addition to the inhibitor; the
  Met20 loop motion this pair is known for closes over NADPH's nicotinamide ring, so it's
  cofactor-driven, not inhibitor-driven. A null would be ambiguous between H2 (mode collapse)
  and "no apo-conditioned method could reach this," which defeats the pilot's purpose.
  Separately, its pocket CA RMSD (2.64 A) exceeds the 2.0 A sampling shell.
- **Ricin A chain (1RTC/1BR6)** -- rejected. Cofactor-free and otherwise clean, but the
  mechanism is a Tyr80 sidechain swing (pocket CA RMSD 0.38 A vs. all-atom 1.05 A, ratio
  0.36) -- backbone sampling is the wrong tool for this motion type, same ambiguity problem
  via a different mechanism.
- **LpqN (6E5D/6E5F)**, the plan's suggested backup -- excluded automatically by the Phase 1b
  screen rather than hand-checked: pocket CA RMSD 0.56 A (below the 1.0 A floor) and ratio
  0.49 (sidechain-gated, ricin's problem again).

## 3. Selected target

**CRBP1 (human cellular retinol-binding protein 1)**, apo **5H9A**, holo **6E5L**. Selected by
`target_screen.py` as the top-ranked qualifying pair (highest CA-to-allatom ratio inside the
band) out of 43 PocketMiner apo/holo pairs -- see `target_screen.csv` for the full screen.
Verified against RCSB:

- 5H9A: 140-residue single-chain construct (chain A), 139/140 residues modeled, no gaps. Only
  heteroatom is BTB (Bis-Tris buffer, a crystallization additive, not a cofactor).
- 6E5L: identical 140-residue sequence (verified, no mutations), 140/140 residues modeled, no
  gaps, chain A. Contains abnormal-cannabidiol (ligand code `HVD`, MW 314.5, drug-like) and no
  other heteroatoms.
- Mechanism: portal-loop opening around residue 76 (max CA displacement 6.22 A there; every
  other pocket residue moves < 1.1 A except residue 77 at 2.04 A) -- a localized backbone
  motion (global CA RMSD 1.15 A vs. pocket CA RMSD 1.57 A, ratio 0.947, the highest of any
  qualifying pair).
- Chain A used from each structure.

**This target stays in the final ten-target panel and in the final writeup regardless of
outcome.** A null result does not remove CRBP1 from the panel; a positive result does not
exempt it from the panel's full statistical treatment. This line is what makes this a pilot
rather than a search over targets for a favorable result.

## 4. Co-primary metrics

1. **`directional_gain`** (Phase 2b), reported on **both** CA and all-atom pocket RMSD:
   `apo_holo_pocket_rmsd - best_conformer_holo_rmsd`, tested against an isotropic-displacement
   null (not against zero), one-sided p-value, for each metric.
2. **`gap_recovery_fraction`** (Phase 4): `(apo_rigid_score - apo_ensemble_score) / (apo_rigid_score - holo_rigid_score)`,
   per ligand, excluding ligands with |gap| < 0.5 kcal/mol.

Neither metric is designated primary over the other; both are reported regardless of sign or
significance. Reporting both CA and all-atom `directional_gain` matters because different
motion types hide their signal in different places -- a CA-only metric would have returned a
false null on ricin's sidechain-gated motion, and an all-atom-only metric blurs backbone
signal with rotamer noise.

## 5. Hypotheses

- **H1 -- sampling reaches alternative states.** Predicts `directional_gain` > 0 (both
  metrics) and partial apo->holo gap recovery in docking.
- **H2 -- mode collapse toward the dominant (likely apo-like) conformation.** Predicts
  `directional_gain` ~ 0 (both metrics) and near-zero gap recovery.
- A third outcome -- the sampler reaches holo-like states but Vina fails to score them
  correctly -- is distinguished from H2 by a positive `directional_gain` co-occurring with
  near-zero gap recovery; this is a scoring-function finding, not a sampling finding.

Neither H1 nor H2 is the desired outcome; both are informative and will be reported as
observed.

## 6. Ligand set (fixed, see `ligands.csv` and `targets.yaml`)

- **Primary:** abnormal-cannabidiol (`Cmpd-P01`), the native ligand of 6E5L, MW 314.5. Has a
  crystallographic pose, so pose RMSD is measurable against ground truth in `holo_rigid`.
- **Power set:** 10 ligands (`Cmpd-S01`-`Cmpd-S10`), real CRBP1-family co-crystal ligands
  (>=95% sequence identity to the pilot construct) pulled from other structures of this same
  protein: two cannabinoid analogs from the same series as the primary ligand, retinylamine
  (a retinoid, CRBP1's natural ligand class), palmitic acid, and six fragment-series
  oxadiazole/piperidine compounds from a CRBP1 fragment-screening deposition. MW 256-390, all
  RDKit-valid. 6E5L is a **non-cognate** holo receptor for these (each was solved with its own
  ligand bound) -- their `holo_rigid` ceiling is approximate, and this is reported as a
  limitation, not glossed over.

11 ligands total, fixed before any docking.

## 7. Tests to be run

- Directional test: `directional_gain` (CA and all-atom) vs. isotropic-displacement null
  (Monte Carlo), one-sided, for each metric.
- Paired t-test and Wilcoxon signed-rank on `apo_ensemble_score` vs. `apo_seeded_score` across
  the 11-ligand set, with bootstrap 95% CI (10,000 resamples) on mean `noise_corrected_delta`.
- Mean `gap_recovery_fraction` and `noise_recovery_fraction`, bootstrap CIs, degenerate-gap
  ligands (|gap| < 0.5 kcal/mol) excluded and reported separately.
- Same analysis repeated with pose RMSD as the outcome variable.

**Explicitly deferred, not computed here:** the cross-target dose-response correlation between
apo->holo pocket displacement and ensemble benefit. That requires the full ten-target panel.
No correlation coefficient spanning multiple targets is to be reported from this pilot.

## 8. Stopping rule

Analysis runs once, after all four docking conditions complete for all 11 ligands. No interim
peeking at partial results to decide whether to continue, add ligands, or change parameters.

## Guardrails (carried forward from Experiment 1)

- Never fabricate results.
- Never tune parameters (box, exhaustiveness, ConforMix twist-target-stop, etc.) to enlarge an
  effect after seeing partial results.
- Log every Vina seed and exhaustiveness value used.
- Report the outcome whichever way it goes; a null result is reported as plainly as a positive
  one, with no post-hoc search across subgroups or alternative metrics to find significance.

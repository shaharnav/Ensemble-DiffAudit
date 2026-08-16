# Pre-registration — Experiment 2 Pilot: Ricin A-Chain Apo/Holo Cross-Docking

This document is committed before any ConforMix generation or docking run for this
experiment. Any number reported for this pilot must postdate the commit hash of this file.

## 0. Target swap (before any data was generated)

The original pilot target was DHFR (2W9T/2W9S). It was swapped for ricin A chain before any
ConforMix generation or docking ran, because the DHFR pair's confound was not just noise but
an alternative explanation for the pilot's primary finding: 2W9S's pocket rearrangement (the
Met20 loop) closes over NADPH's nicotinamide ring and is cofactor-driven, not
inhibitor-driven, so an apo-conditioned sampler that never sees NADPH could not reach that
conformation regardless of whether the sampling method works. A null result on that pair would
have been ambiguous between H2 (mode collapse — informative) and "the target conformation
requires a cofactor" (uninformative, true by construction) — indistinguishable, which defeats
the pilot's purpose. See `targets.yaml`'s `target_swap_history` for the full reasoning. This is
a stricter bar than "document it or pick a matched pair": a confound that supplies an
alternative explanation for the primary finding must be removed, not merely disclosed.

## 1. Pilot target

Ricin A chain (*Ricinus communis*), apo **1RTC**, holo **1BR6**. Verified against RCSB:

- 1RTC: 268-residue single-chain construct (chain A), full 268/268 residues modeled, no gaps.
  HETATM records contain only crystallographic waters — genuinely ligand-free and
  cofactor-free.
- 1BR6: identical 268-residue sequence (verified, no mutations), full 268/268 residues modeled,
  no gaps, chain A. Contains pteroic acid (ligand code `PT1`, MW 312.3, drug-like) and waters
  only — **no cofactor in either structure**, so the only difference between apo and holo is
  the inhibitor.
- Mechanism: Tyr80 sidechain swings out of the way to open the active site on inhibitor
  binding — documented in the CryptoSite cryptic-pocket benchmark, so this is a published
  benchmark pair rather than one selected post hoc for a favorable outcome.
- Chain A used from each structure.

## 2. Commitment

**This target stays in the final ten-target panel and in the final writeup regardless of
outcome.** A null result does not remove ricin A chain from the panel; a positive result does
not exempt it from the panel's full statistical treatment. This line is what makes this a
pilot rather than a search over targets for a favorable result.

## 3. Co-primary metrics

1. **`directional_gain`** (Phase 2b): `apo_holo_pocket_ca_rmsd − best_conformer_holo_rmsd`,
   tested against an isotropic-displacement null (not against zero), one-sided p-value.
2. **`gap_recovery_fraction`** (Phase 4): `(apo_rigid_score − apo_ensemble_score) / (apo_rigid_score − holo_rigid_score)`,
   per ligand, excluding ligands with |gap| < 0.5 kcal/mol.

Neither metric is designated primary over the other; both are reported regardless of sign or
significance.

## 4. Hypotheses

- **H1 — sampling reaches alternative states.** Predicts `directional_gain` > 0 and partial
  apo→holo gap recovery in docking.
- **H2 — mode collapse toward the dominant (likely apo-like) conformation.** Predicts
  `directional_gain` ≈ 0 and near-zero gap recovery.
- A third outcome — the sampler reaches holo-like states but Vina fails to score them
  correctly — is distinguished from H2 by a positive `directional_gain` co-occurring with
  near-zero gap recovery; this is a scoring-function finding, not a sampling finding.

Neither H1 nor H2 is the desired outcome; both are informative and will be reported as
observed.

## 5. Ligand set (fixed, see `ligands.csv` and `targets.yaml`)

- **Primary:** pteroic acid (`Cmpd-P01`), the native ligand of 1BR6, MW 312.3. Has a
  crystallographic pose, so pose RMSD is measurable against ground truth in `holo_rigid`.
- **Power set:** 10 ligands (`Cmpd-S01`–`Cmpd-S10`), known ricin A-chain inhibitors drawn from
  other ricin A-chain co-crystal structures in the PDB (pteridine-7-carboxamide and
  pyrimidinone-based active-site binders, plus one thiophene-carboxylic-acid scaffold for
  diversity), MW 207–411, all RDKit-valid. 1BR6 is a **non-cognate** holo receptor for these
  (each was solved with its own inhibitor bound, not pteroic acid) — their `holo_rigid`
  ceiling is approximate, and this is reported as a limitation, not glossed over.

11 ligands total, fixed before any docking.

## 6. Tests to be run

- Directional test: `directional_gain` vs. isotropic-displacement null (Monte Carlo), one-sided.
- Paired t-test and Wilcoxon signed-rank on `apo_ensemble_score` vs. `apo_seeded_score` across
  the 11-ligand set, with bootstrap 95% CI on mean `noise_corrected_delta`.
- Mean `gap_recovery_fraction` and `noise_recovery_fraction`, bootstrap CIs, degenerate-gap
  ligands (|gap| < 0.5 kcal/mol) excluded and reported separately.
- Same analysis repeated with pose RMSD as the outcome variable.

**Explicitly deferred, not computed here:** the cross-target dose-response correlation between
apo→holo pocket displacement and ensemble benefit. That requires the full ten-target panel.
No correlation coefficient spanning multiple targets is to be reported from this pilot.

## 7. Stopping rule

Analysis runs once, after all four docking conditions complete for all 11 ligands. No interim
peeking at partial results to decide whether to continue, add ligands, or change parameters.

## 8. Guardrails (carried forward from Experiment 1)

- Never fabricate results.
- Never tune parameters (box, exhaustiveness, ConforMix twist-target-stop, etc.) to enlarge an
  effect after seeing partial results.
- Log every Vina seed and exhaustiveness value used.
- Report the outcome whichever way it goes; a null result is reported as plainly as a positive
  one, with no post-hoc search across subgroups or alternative metrics to find significance.

# Pre-registration — Experiment 2 Pilot: DHFR Apo/Holo Cross-Docking

This document is committed before any ConforMix generation or docking run for this
experiment. Any number reported for this pilot must postdate the commit hash of this file.

## 1. Pilot target

DHFR (*Staphylococcus aureus*), apo **2W9T**, holo **2W9S**. Verified against RCSB:

- 2W9T: 161-residue single construct, chains A/B, HETATM records contain only crystallographic
  waters (136 HOH) — genuinely ligand-free and cofactor-free.
- 2W9S: same 161-residue sequence (verified identical, no mutations), chains A–F, contains
  trimethoprim (ligand code `TOP`, MW 290.3, drug-like) plus NADPH (`NDP`) and glycerol (`GOL`,
  cryoprotectant, ignored).
- Chain A used from each structure.

**Confound, retained rather than resolved by re-selecting a target:** 2W9S carries NADPH: 2W9T
does not. Any measured apo→holo pocket rearrangement may be driven in part or in whole by
cofactor binding rather than solely by trimethoprim binding. This is reported alongside every
geometry and docking result in this pilot, not silently absorbed into the headline number.

## 2. Commitment

**This target stays in the final ten-target panel and in the final writeup regardless of
outcome.** A null result does not remove DHFR from the panel; a positive result does not
exempt it from the panel's full statistical treatment. This line is what makes this a pilot
rather than a search over targets for a favorable result.

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

- **Primary:** trimethoprim (`Cmpd-P01`), the native ligand of 2W9S, MW 290.3. Has a
  crystallographic pose, so pose RMSD is measurable against ground truth in `holo_rigid`.
- **Power set:** 10 ligands (`Cmpd-S01`–`Cmpd-S10`), known DHFR inhibitors drawn from other
  DHFR co-crystal structures in the PDB (methotrexate, trimetrexate, brodimoprim-derivative,
  and several pyridopyrimidine/phthalazinone-series inhibitors), MW 266–497, all RDKit-valid.
  2W9S is a **non-cognate** holo receptor for these — their `holo_rigid` ceiling is
  approximate, and this is reported as a limitation, not glossed over.

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

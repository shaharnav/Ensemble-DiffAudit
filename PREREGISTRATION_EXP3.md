# Pre-registration — Experiment 3: Xylellain De Novo Design Campaign

This document is committed before any ConforMix or DiffSBDD generation, and before any
docking, for this experiment. Any number reported for this experiment must postdate the
commit hash of this file.

## Claims boundary — reproduce verbatim in the writeup, near the top, not in an appendix

This is **not** a validation study. Experiments 1 and 2 asked whether ensemble docking works,
measured against crystallographic ground truth. This experiment has no ground truth: there is
no holo structure of Xylellain, no known potent inhibitor, and no experimental measurement of
anything produced here.

**Supportable without a holo structure:**
- Relative docking scores across receptor states (occluded vs. ensemble vs. propeptide-removed)
- Magnitude of generation home-field bias, in kcal/mol
- Whether ensemble conformers open the pocket, measured as cavity volume
- Whether ensemble best-of-N exceeds a seed-matched rigid baseline
- Physicochemical and validity properties of generated compounds

**Not supportable, and must not be implied:**
- That any generated compound binds Xylellain
- That any docked pose is correct
- That the ensemble conformers correspond to conformations the protein actually populates
- Any statement about potency, selectivity, or efficacy
- That this constitutes a drug candidate, hit, or lead

**Language rule:** write "scored favorably in docking," never "binds," "inhibits," or "is a
hit."

## Receptor set and ligand counts (fixed in advance)

- **Receptors (8 total):** R1 (zymogen, occluded, as deposited), R2 (propeptide-removed,
  unrelaxed deletion, a modeled ceiling proxy, not a measured structure), C0-C5 (6 ConforMix
  conformers generated from R1). See `receptors.yaml`.
- **Generated ligands (45 total):** 15 each conditioned on R1, C_mid (median-cavity-volume
  conformer), and R2. Every ligand tagged with `conditioning_receptor` — this tag is what
  makes the home-field bias measurement (Phase 7) possible.
- **Reference compounds (5):** E-64 (epoxide), leupeptin (aldehyde), a vinyl-sulfone-derived
  cruzain inhibitor, a nitrile-derived cathepsin K inhibitor, an azanitrile-derived SmCB1
  inhibitor — real PDB-deposited ligands, not ground truth for Xylellain, used only as a
  scale check. See `reference_compounds.csv`. Three of the five (vinyl sulfone, nitrile,
  azanitrile) are deposited in their post-covalent-reaction form, not the free pre-reaction
  warhead — recorded per-compound in `reference_compounds.csv`'s `note` column, since this
  pipeline docks them non-covalently regardless.
- **Docking matrix:** 45 ligands x 8 receptors (360) + 45 ligands x R1 x 6 seeds (270) + 5
  reference compounds x 8 receptors (40) = 670 jobs.

## Structural verification (Phase 1a, see `verify_3ois.py` / `receptors.yaml`)

- Ribonucleotide (UDP) is 19.98-24.81 A from the active site — well clear of the 8 A gate.
  Active site confirmed apo, consistent with the paper's own abstract.
- Catalytic triad: Cys78, His237, Asn255, confirmed by 3D proximity (Cys78 SG - His237 NE2:
  4.45 A; His237 ND1 - Asn255 OD1: 3.58 A), not assumed from papain's own numbering.
- Propeptide boundary: residues 23-66, determined from which residues actually approach the
  active-site centroid (occluding loop at 26-36, 3.9-9.8 A) rather than assumed as a fixed
  count from the abstract. Full reasoning in `receptors.yaml`.
- Pocket volume gate (Phase 2, grid-based proxy — fpocket unavailable in this environment):
  R1 2454.2 A^3, R2 2744.5 A^3, both well above the ~150 A^3 threshold. Gate passed.

## Metric definitions (Phase 7)

```
home_advantage        = mean(score in own conditioning receptor) - mean(score in the other 7 receptors)   # per ligand, paired
noise_corrected_delta = R1_seeded_max - ensemble_best                                                       # per ligand, positive = ensemble better
opening_recovery      = (R1_rigid - ensemble_best) / (R1_rigid - R2_rigid)                                  # degenerate-gap guard: exclude |R1_rigid - R2_rigid| < 0.5 kcal/mol
```

Sign convention from Experiment 1: affinities are negative; a positive delta means the
ensemble scored more negative (better) — asserted in a unit test before use.
`opening_recovery` always carries R2's modeling-assumption caveat wherever it's reported.

## Reporting commitments

1. All 45 generated compounds are reported, not just top scorers.
2. The propeptide-deletion modeling assumption (Phase 1c) and its limitation are stated
   everywhere `opening_recovery` or R2 appears, not once and forgotten.
3. Structural alerts (covalent warheads) are reported but interpreted in context — Michael
   acceptors, nitriles, epoxides are legitimate cysteine-protease chemotypes here (E-64 is an
   epoxide), not automatic liabilities.
4. Attrition ledger (generated -> parseable -> embedded -> PoseBusters-valid -> docked ->
   in final results) reported with per-stage failure reasons.
5. Spearman correlation between conformer cavity volume and mean docking score is computed
   and reported explicitly as a score-inflation check, whichever way it comes out.

## Stopping rule

Analysis runs once, after all docking (Phase 5) completes for all receptors and all ligands.
No interim peeking at partial results to decide whether to continue, add ligands, or change
parameters.

## Guardrails (carried forward from Experiments 1 and 2)

- Never fabricate results.
- Never tune parameters (box, exhaustiveness, ConforMix twist-target-stop, generation
  conditioning) to enlarge an effect after seeing partial results.
- Log every Vina seed and exhaustiveness value used.
- Report negative and null outcomes as plainly as positive ones; no post-hoc search across
  subgroups or alternative metrics for significance.

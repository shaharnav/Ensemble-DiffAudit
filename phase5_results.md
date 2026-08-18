# Phase 5 Results: Ensemble Cross-Docking Audit — LasB (1EZM)

## Context

Phase 4 found a headline delta of +1.42 kcal/mol (best DiffSBDD candidate ensemble max vs. crystal),
with 20/20 positive deltas and a noise-corrected mean delta of +0.990 kcal/mol across 140 jobs
(20 candidates × 7 receptors; β4.0 failed silently and was excluded).

Phase 5 tests whether that signal is real via three pre-registered controls.

Pre-registration: `phase5_preregistration.md`, committed `1415229` before any Phase 5 docking.

---

## Step 0 — Reanalysis of Phase 4 Data

All statistics on the 140-job Phase 4 matrix, no new docking.

| Metric | Value |
|--------|-------|
| Headline delta (max-over-5 conformers) | +1.42 kcal/mol |
| Mean noise-corrected delta (mean-over-5) | +0.990 kcal/mol |
| Conformer win distribution | β0.0: 6, β0.8: 4, β1.6: 3, β2.4: 4, β3.2: 3 |
| Per-conformer means (candidates) | −6.07 to −6.29 kcal/mol (range: 0.22 kcal/mol) |
| Crystal mean (candidates) | −5.21 kcal/mol |

**Flag:** Per-conformer means are uniformly ~−6.2 kcal/mol regardless of β value — no gradient,
no directional gain. All five conformers show nearly identical affinity, suggesting a receptor-level
scoring offset rather than geometric recognition of holo-like geometry.

**Table 1 — Conformer means (Phase 4, candidates)**

| Receptor | Mean affinity (kcal/mol) |
|----------|--------------------------|
| 1EZM crystal | −5.214 |
| β0.0 | −6.286 |
| β0.8 | −6.259 |
| β1.6 | −6.118 |
| β2.4 | −6.289 |
| β3.2 | −6.068 |

---

## Control A — Seeded Crystal (Noise Floor)

**Protocol:** 20 DiffSBDD candidates × 6 seeds (42, 123, 456, 789, 1337, 2024) in 1EZM crystal.
120 jobs total. Exhaustiveness 16, same pocket center and box.

**Purpose:** Measure Vina stochasticity to separate real signal from seed noise.

| Metric | Value |
|--------|-------|
| Mean seed noise (best-of-6 vs seed-42) | +0.046 kcal/mol |
| Max single-candidate seed range | +0.109 kcal/mol |
| Noise-corrected advantage (cand mean) | +1.372 kcal/mol |
| Bootstrap 95% CI | [+1.191, +1.542] |

**Verdict: SURVIVES.** The signal (+1.372 kcal/mol) is 30× the noise floor (+0.046 kcal/mol).
Vina stochasticity is not the explanation. CI lower bound (+1.191) far exceeds the pre-registered
falsification threshold of ≤ +0.30 kcal/mol.

---

## Control B — Property-Matched Decoys

**Protocol:** 20 PubChem decoys (MW ±25 Da, cLogP ±0.5, HBD ±1, HBA ±1, RotB ±2,
Tanimoto < 0.30 vs. all 20 candidates) × 6 receptors (5 conformers + 1EZM crystal) = 120 jobs.
Same settings as Phase 4 throughout.

**Purpose:** Distinguish holo geometric recognition from a receptor-level scoring offset.

**Table 2 — Aggregate comparison**

| Metric | Candidates | Decoys |
|--------|-----------|--------|
| Mean advantage over crystal (max-based) | +1.372 kcal/mol | +1.616 kcal/mol |
| Mean delta_mean | +0.990 kcal/mol | +1.217 kcal/mol |
| discrimination_gap (max, cand − decoy) | −0.244 kcal/mol | — |
| Bootstrap 95% CI | [−0.597, +0.062] | — |

**Table 3 — Per-receptor means: decoys vs. candidates**

| Receptor | Candidates | Decoys | Diff (cand − decoy) |
|----------|-----------|--------|---------------------|
| 1EZM crystal | −5.214 | −5.475 | +0.261 |
| β0.0 | −6.286 | −6.896 | +0.610 |
| β0.8 | −6.259 | −6.659 | +0.401 |
| β1.6 | −6.118 | −6.526 | +0.408 |
| β2.4 | −6.289 | −6.624 | +0.335 |
| β3.2 | −6.068 | −6.754 | +0.686 |

Decoys score better than candidates on every receptor.

**Pre-registered threshold:** discrimination_gap < +0.30 → FALSIFIED
**Observed:** −0.244 [CI: −0.597, +0.062]
**Verdict: FALSIFIED.** The conformer advantage applies equally (or more) to property-matched
strangers. The Boltz conformers apply a uniform scoring uplift to any molecule with appropriate
physical properties — they are not recognizing holo-like geometry specific to DiffSBDD candidates.

---

## Control C — True Holo Crystal Ceiling

**Protocol:** Candidates and decoys already docked into 1EZM crystal in Phase 4 and Control B.
Extracted scores for both populations in the holo-rigid crystal (no new docking).

**Purpose:** Establish whether the pocket can discriminate at all in the best possible structure.
If decoys match or beat candidates in the true holo crystal, no conformer ensemble over this
pocket could ever have passed Control B — the discrimination failure is upstream of the conformers.

| Population | Mean | Median | SD |
|------------|------|--------|----|
| DiffSBDD candidates (n=20) | −5.214 | −5.153 | 0.687 |
| Decoys (n=20) | −5.475 | −5.579 | 0.754 |
| Gap (cand − decoy) | +0.261 | — | — |
| Bootstrap 95% CI | [−0.182, +0.686] | — | — |

Decoys beat candidates by 0.261 kcal/mol on average. 13/20 decoys (65%) beat the candidate median.
CI spans zero — no statistically significant discrimination in either direction.

**Verdict:** The LasB pocket cannot distinguish DiffSBDD-generated candidates from property-matched
strangers even in the experimentally determined holo crystal. This is a ceiling problem, not a
conformer problem.

---

## Summary

| Control | Question | Threshold | Observed | Verdict |
|---------|----------|-----------|----------|---------|
| A (noise floor) | Is signal > Vina noise? | noise_corrected_delta ≤ +0.30 → falsified | +1.372 [+1.191, +1.542] | SURVIVES |
| B (decoys) | Is signal specific to candidates? | discrimination_gap < +0.30 → falsified | −0.244 [−0.597, +0.062] | FALSIFIED |
| C (holo ceiling) | Can pocket discriminate at all? | — | gap +0.261 [−0.182, +0.686], CI spans zero | NO DISCRIMINATION |

## Conclusion

The +1.37 kcal/mol noise-corrected advantage from Phase 4 is real (survives noise floor) but
non-specific (decoys gain equally or more). The Boltz conformers produce a uniform scoring uplift
attributable to receptor softening, not holo-geometry recognition. Control C confirms the failure
is baked into the pocket itself: the 1EZM/holo crystal cannot discriminate DiffSBDD candidates
from property-matched strangers. No ensemble over this pocket could have passed Control B.

The discriminability problem sits at the level of either the pocket geometry or what DiffSBDD
generates against it — both are upstream of the conformer step. Vina scores in this pocket
carry no candidate-specific signal.

**Next step (Phase 6):** Calibrate whether Vina scores carry *any* affinity information in this
pocket by correlating with published IC50/Ki for known LasB inhibitors. If no signal at all
(|r| < 0.4), the null result explains the Control B failure at its mechanistic root.

# Phase 2b — ConforMix Conformer Quality & Directional Gain Results

**Committed before DiffSBDD generation starts for LasB.**
Date: 2026-08-17

This document records the Phase 2b gate outcome for both candidate targets (GluR and LasB)
as evaluated by MSA-enabled ConforMix/Boltz runs. It is committed to the repository before
any generation or docking numbers are observed, so the verdict cannot be reverse-engineered
from a downstream docking outcome.

---

## Pre-registered 2×2 framework

Phase 2b evaluates two independent dimensions for each candidate ensemble:

| | **Directional gain < 0** (conformers move *away* from holo) | **Directional gain ≥ 0** (conformers move *toward* holo) |
|---|---|---|
| **Global RMSD ≤ 2.5 Å** (structurally plausible) | H0: conformers are valid but don't sample holo-like states. Generation proceeds on the apo pocket only; ensemble advantage not expected. | **H1 holds**: conformers are structurally plausible AND biased toward the holo pocket. Ensemble cross-docking expected to be informative. |
| **Global RMSD > 2.5 Å** (collapsed/misfolded sample) | Eliminated by hard gate; not counted. | Eliminated by hard gate; not counted. |

The hard gate (global CA RMSD to apo crystal structure ≤ 2.5 Å) was pre-registered during
this session before smoke-test results were observed. It is evaluated first; conformers that
fail it do not contribute to the directional gain verdict.

---

## Methodology note — MSA

All ConforMix/Boltz runs in Phase 2b used `use_msa_server=True` (MMseqs2 server via the
ColabFold API), applied via a sed patch to `boltz/run_twisted.py`. This was required after
the single-sequence baseline produced global RMSD of 8–9 Å for GluR (far outside the hard
gate), demonstrating that single-seq Boltz cannot reliably reproduce the crystal structure
frame for these targets. MSA-enabled runs brought global RMSD to 1.0–2.5 Å for all
passing conformers. The single-seq GluR run is retained as a documented methodology negative
in the notebook but is not used for any Phase 2b verdict.

---

## GluR (S. pyogenes Glutamate Racemase)

**Pair:** 2OHG (apo) / 2OHV (holo, NHL inhibitor)
**Pocket:** 49 residues, subset `9-15,33-40,46,49-57,71-77,79-80,116-121,148-149,152,182-185,239,264`
**Apo→holo baseline pocket CA RMSD:** 1.397 Å (n=49)
**Run settings:** `--twist-target-start -2.0 --twist-target-stop 2.0 --num-twist-targets 6 --samples-per-target 1 --structured-regions-only`

| β (Å) | Global RMSD | Hard gate | Pocket→apo (Å) | Pocket→holo (Å) | Directional gain (Å) | Pass |
|------:|------------|-----------|---------------|----------------|---------------------|------|
| −2.00 | 2.257 | ✓ | 2.654 | 2.007 | −0.610 | PASS |
| −1.20 | **4.973** | **✗** | 5.490 | 5.840 | −4.443 | **FAIL** |
| −0.40 | 2.347 | ✓ | 2.637 | 1.895 | −0.498 | PASS |
| +0.40 | 2.267 | ✓ | 2.574 | 1.860 | −0.463 | PASS |
| +1.20 | 2.451 | ✓ | 2.728 | 2.040 | −0.643 | PASS |
| +2.00 | 2.225 | ✓ | 2.538 | 1.792 | −0.395 | PASS |

**Passing conformers:** 5/6 (β=−1.20 collapsed, eliminated by hard gate).
**Directional gain across all 5 passing conformers:** −0.395 to −0.643 Å (all negative).

**2×2 verdict: Top-left cell — H0.**

The 5 structurally plausible conformers all move the pocket *away* from the holo state
rather than toward it. The twist guidance acts orthogonal or opposite to the apo→holo
direction. Root cause: GluR's target motion (pocket CA RMSD 1.397 Å, qualifies_=True in
`target_screen.csv`) is annotated as motion_type=`loop` — the conformational change is
loop-mediated rather than helix/sheet rearrangement. ConforMix's RMSD guidance acts on
secondary-structure elements (consistent with `structured_regions_only=True`), so it cannot
directly steer loop-localized motions. The 5 passing conformers are structurally valid by
the hard gate but do not sample the holo-like pocket geometry. GluR is retained in the
repository as a documented negative with a known mechanistic cause; no generation is run
against it.

---

## LasB (P. aeruginosa LasB elastase)

**Pair:** 1EZM (apo) / 3DBK (holo, RDF inhibitor)
**Pocket:** 40 residues, subset `110-116,121-122,128-130,132-134,136-142,144,155,160,163-164,167-168,186-187,190-191,197-198,221-224,226`
**Apo→holo baseline pocket CA RMSD:** 1.414 Å (n=40)
**Run settings:** `--twist-target-start 0.0 --twist-target-stop 4.0 --num-twist-targets 6 --samples-per-target 1 --structured-regions-only`

| β (Å) | Global RMSD | Hard gate | Pocket→apo (Å) | Pocket→holo (Å) | Directional gain (Å) | Pass |
|------:|------------|-----------|---------------|----------------|---------------------|------|
| +0.00 | 1.087 | ✓ | 1.270 | 0.372 | +1.042 | PASS |
| +0.80 | 1.341 | ✓ | 1.442 | 0.408 | +1.006 | PASS |
| +1.60 | 0.990 | ✓ | 1.096 | 0.491 | +0.923 | PASS |
| +2.40 | 1.099 | ✓ | 1.206 | 0.529 | +0.885 | PASS |
| +3.20 | 1.248 | ✓ | 1.399 | 0.449 | +0.965 | PASS |
| +4.00 | 1.131 | ✓ | 1.180 | 0.730 | +0.683 | PASS |

**Passing conformers:** 6/6.
**Directional gain across all 6 conformers:** +0.683 to +1.042 Å (all positive).
**Pocket→holo RMSD:** 0.372–0.730 Å across all conformers (all substantially closer to holo than apo is to holo at baseline 1.414 Å).

**2×2 verdict: Top-right cell — H1 holds.**

All 6 conformers are structurally plausible (global RMSD <1.4 Å, well within the hard gate)
and all move the pocket toward the holo reference. The twist guidance is working as intended:
positive and consistent directional gain across the full β range from 0.0 to 4.0 Å. The
pocket residues land 0.37–0.73 Å from the holo reference after alignment, compared to 1.41 Å
for the unperturbed apo — the ensemble spans a range between apo and beyond-holo geometry.
This is the precondition for ensemble cross-docking to be informative: ligands that score
better in the holo-like conformers than in the apo conformer will reveal a pocket-geometry
preference that the rigid-receptor baseline cannot capture.

**LasB is selected for all subsequent steps.** DiffSBDD generation proceeds against the
LasB apo pocket (1EZM). Ensemble cross-docking will be run against all 6 conformers.

---

## Stopping rule as pre-registered

If both targets had landed in the top-left cell (H0), generation would not have proceeded
and this document would have reported a methodology negative: ConforMix does not reliably
steer toward holo geometry for loop-dominated motions, regardless of MSA quality.
GluR demonstrates that case. LasB demonstrates the positive case. Both are reported without
hedging; neither result is selected after the fact.

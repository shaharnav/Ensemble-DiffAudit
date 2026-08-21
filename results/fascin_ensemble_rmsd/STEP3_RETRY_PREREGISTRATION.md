# Fascin Step 3 retry: diagnostic + pre-registered generation-settings correction

Committed before generating anything for the retry, per instruction: a settings
change made after seeing a failure is legitimate only if it corrects a
demonstrated misconfiguration and its acceptance criterion is fixed in advance,
not tuned to whatever comes back.

## Diagnostic: why did all 16 conformers fail?

475 clashes (<2.0A between heavy atoms of residues >=4 apart in sequence,
excluding the crystallographically-unresolved N-terminal residues 1-3) across
16 conformers. Two hypotheses tested directly, not assumed:

**Hypothesis A -- repeat confusion.** Fascin has 4 tandem beta-trefoil domains
(CATH 2.80.10.50, boundaries confirmed via the PDBe CATH mapping for 3LLP:
D1=7-138, D2=139-261, D3=262-386, D4=387-493). If RMSD-guided sampling can't
distinguish which structural repeat it's pulling toward the target, clashes
should concentrate at *structurally equivalent* positions across different
domains (residue N in D1 clashing with residue N in D2, etc).

Tested by mapping every clash to (domain, offset-within-domain) and comparing
the offset-difference distribution against a random-pair null:
- Inter-domain clashes with offset-difference <=10 (i.e. "same" position in
  two domains): **8/69 (12%)**
- Same statistic under a random-pairs null: **222/1499 (15%)**

No enrichment -- indistinguishable from chance. **Hypothesis A rejected.**

**Hypothesis B -- local packing collapse under aggressive guidance.** Of all
475 clashes, **355 (75%) are intra-domain**, not inter-domain. The dominant
failure mode is a single domain's own side chains/backbone breaking their own
packing, not domains overlapping each other. This is consistent with the
guidance force being too large relative to what a domain of this size can
absorb without breaking local geometry -- not a repeat-specific effect, and
not something a global CA RMSD gate would ever catch (13/16 conformers passed
global RMSD <=2.5A while every one of the 16 had severe local clashes; global
RMSD is measuring the right ballpark fold and is blind to this failure mode
entirely -- worth carrying into the final writeup regardless of retry outcome).

## The actual misconfiguration

`--subset-residues` covered **139/493 residues (28.2%)** of fascin, carried
over unchanged in *absolute* radius (pocket residue +/-10) from the LasB run.
LasB's own validated subset covered **40/298 residues (13.4%)** of a
single-domain, smaller protein. Fascin's run used a guidance footprint over
**2x LasB's validated fractional coverage**, on a protein with four repeated
domains where aggressive per-residue guidance has more surface area to break.
That's an untested extrapolation of a setting, not a hypothesis tested and
rejected -- fixing it is a legitimate correction.

## Pre-registered retry

**Fixed, unchanged from the original run and from each other:**
- Clash gate: <2.0A between heavy atoms of residues >=4 apart in sequence
  (excluding unresolved N-term 1-3), same as just applied.
- Global CA RMSD gate: <=2.5A to apo crystal (3LLP), same as just applied.
- Reference structure (6I11), twist-target range (0.0-6.0A).

**Changed, both scaled together, not just one:**
- `--subset-residues`: shrunk to pocket +/-5 (was +/-10) --
  `9-21,43-53,55-65,88-108,129-139,209-229`, **88/493 residues (17.8%)** --
  close to LasB's validated 13.4%, not tuned further than that match.
- `--twist-strength`: ladder of **{15.0 (baseline, already have this data --
  0/16 clash-free, not re-run), 7.5, 3.75}**.

**Ladder protocol** (bounded, not open-ended parameter search): at each new
twist_strength value (7.5, 3.75), generate a **4-target smoke test**
(beta = 0.0, 2.0, 4.0, 6.0 -- sparse coverage of the same range, not the full
16) with the corrected subset window. Record clash count per conformer at
each setting. This is capped at 2 new settings x 4 targets = 8 additional
generations total, not an open search.

**Acceptance criterion, fixed now:**
- If either new twist_strength setting produces **>=3 of 4** smoke-test
  conformers passing both gates (global RMSD gate is already usually passing;
  the real bar is clash-free), regenerate the full 16-target sweep at that
  setting and proceed to the rest of Step 3/4 as originally pre-registered
  (target >=8 of 16 survive; report actual attrition either way).
- If neither setting clears that bar, **this is the final result for fascin
  under ConforMix/Boltz twist-guidance**: report it as a genuine scope
  limit (large, multi-domain, internally-repeated folds are not well served
  by this guidance approach at any tested setting), together with the
  repeat-confusion-rejected / local-packing-collapse finding above, and the
  global-RMSD-blind-to-local-clash observation. No further settings will be
  tried beyond this ladder.

Committed before generating any of the 8 ladder structures.

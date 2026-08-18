# Phase 7 Results: Retrospective Rescoring of LasB Docking Data

No new docking was performed anywhere in this phase. Every score below was computed from poses
already generated in Phases 4-6 and Control A/B (426 pose files with coordinates on disk).

**Scatter plots and ROC curves**: https://claude.ai/code/artifact/7e27aec0-53e5-46a6-8e0d-3cd68b5bc7ab

**Deliverable files**: `results/phase7_tidy_scores.csv` (328-row tidy score table), `results/phase7_vinardo_scores.csv`
(all 426 raw Vinardo scores), `results/phase7_posebusters.csv`, `results/phase7_strain.csv`,
`results/phase7_reproducibility.csv`.

**Environment**: macOS arm64, CPU only, no CUDA/GPU (per instruction).

**Interpretive framing (stated up front per instructions)**: LasB is zinc-confounded (Phase 6).
A positive ΔvinaXGB or Vinardo result would be informative — it would show a correction captures
something Vina misses. A null result here is **not** evidence against those methods; it is
indistinguishable from the zinc-coordination limitation and from N=11/1.49-log-unit statistical
power, both established in Phase 6. Nulls below are reported as uninformative, not disconfirming.

---

## Step 0/3 — DeltaVina Install: INFEASIBLE (documented attempt)

Cloned the official repo (`github.com/jenniening/deltaVinaXGB`, the ΔvinaXGB implementation from
Lu et al. 2019). Its full dependency chain — a modified AutoDock Vina fork (`vina4dv`, Linux
x86_64 binary only, needed to compute the "Vina58" 58-term pairwise energy decomposition that
standard Vina's `--score_only` does not expose), MGLTools (Linux/Mac tarball, Python 2, needed for
`prepare_receptor4.py`/`prepare_ligand4.py`), and MSMS (Linux x86_64 binary, SASA features) — has
no native ARM build path. All three external download servers (mgltools.scripps.edu twice, GitHub
for vina4dv) were reachable, so a Docker build (`--platform linux/amd64`, emulated) was attempted
as the only viable path, using the repo's own Dockerfile (Ubuntu 16.04 base).

Three sequential, distinct build failures were diagnosed and patched, each a genuine
version-compatibility dead end rather than a typo:

1. **GLIBC mismatch**: current Miniconda3 installer requires GLIBC ≥2.28; Ubuntu 16.04 ships 2.23.
   Patched to pin Miniconda3 4.5.4 (contemporary with the 2019 codebase).
2. **`conda run` unavailable**: Miniconda 4.5.4's conda predates the `conda run` subcommand the
   Dockerfile's `SHELL` directive requires. Patched to `source activate` in each `RUN` instead.
3. **Malformed version string from conda-forge**: with the old conda client patched in, the
   `conda install -c conda-forge xgboost=0.80.0 ...` step failed with `CondaValueError: Malformed
   version string '~'` — modern conda-forge/repodata metadata uses syntax the 2018-era conda
   client cannot parse. This is a hard dead end without vendoring a period-correct package index,
   which was judged out of scope for this session.

**Verdict: ΔvinaXGB and ΔvinaRF20 (which shares the same vina4dv/MGLTools/MSMS dependency) could
not be installed in this environment.** No score_only run was ever attempted against a PDBbind
verification complex, so the mandatory verification gate was never reached — this is reported as
an install failure, not a verification failure. Per instructions, this is not counted as a null
result against ΔvinaXGB; it is simply unavailable evidence.

**Fallback executed**: Vinardo, which the local `vina_1.2.7_mac_aarch64` binary supports natively
via `--scoring vinardo --score_only` — no smina install needed. All 426 existing poses (Phase 4
ensemble, Control A, Control B, Phase 6 known actives) were rescored this way with zero failures.
Scores in `results/phase7_vinardo_scores.csv`.

---

## Step 1 — Free Reanalysis (no rescoring needed)

**1a. Heavy-atom-count (HAC), candidates vs decoys** (Phase 4/Control B population, n=20 each):

| | Mean HAC | SD | n |
|---|---|---|---|
| Candidates | 20.30 | 5.11 | 20 |
| Decoys | 19.65 | 5.21 | 20 |

t = 0.398, p = 0.693, Cohen's d = 0.126. **No significant HAC difference** between the two
populations — the decoy-matching procedure worked as intended on this axis, so the Control B
result is not an artifact of the decoys being systematically smaller or larger.

**1b. Control B in ligand-efficiency units (score / HAC)**:

| Metric | Raw kcal/mol | Ligand efficiency |
|---|---|---|
| Candidate mean advantage (crystal → best conformer) | +1.417 | +0.075 |
| Decoy mean advantage | +1.615 | +0.085 |
| discrimination_gap (cand − decoy) | −0.198 | −0.0101 |
| Bootstrap 95% CI | [−0.553, +0.118] | [−0.0296, +0.0090] |

(Note: this segment's raw-kcal/mol recompute of −0.198 differs slightly from Phase 5's reported
−0.244 because Phase 5 used the Control-A seed-corrected candidate mean; both agree well within
CI and support the same conclusion.)

**The decoy advantage persists in ligand-efficiency units — it does not reverse or shrink to
insignificance.** The Control B falsification is not an HAC artifact. No amendment to Phase 5's
conclusion is required.

---

## Step 2 — Pose Inventory

| Pose set | Expected | On disk (with coordinates) | Status |
|---|---|---|---|
| Phase 4 (20 cand × 6 receptors incl. crystal) | 120 | 120 | complete |
| Control A (20 cand × 6 seeds, crystal) | 120 | 120 | complete |
| Control B (20 decoys × 6 receptors) | 120 | 120 | complete |
| Phase 6 (11 actives × 2 receptors × 3 seeds) | 66 | 66 | complete |
| **Holo test — 20 cand + 20 decoys × 3DBK rigid** | **40** | **0** | **MISSING — see below** |
| **Total** | | **426** | |

All 426 pose files carry full coordinates (not just scores), sufficient for rescoring without
re-docking.

### Critical correction to Phase 5 "Control C"

**Phase 5's "Control C — True Holo Crystal Ceiling" was never actually run against 3DBK.**
Inventory confirms zero candidate or decoy poses exist for 3DBK anywhere on disk — the only 3DBK
poses that exist are the 33 Phase 6 known-active poses (11 compounds × 3 seeds). Re-reading
`phase5_results.md`'s own Control C protocol line: *"Candidates and decoys already docked into
**1EZM crystal** in Phase 4 and Control B... Extracted scores for both populations in the
holo-rigid crystal (no new docking)."* — this reused **1EZM (apo)** scores and mislabeled them as
the "true holo crystal ceiling." 1EZM has no bound ligand (its only HETATM records are ZN, CA, and
water); 3DBK is the actual holo structure. Per this phase's explicit rule ("No new docking. If a
needed pose is missing, report it rather than re-docking"), this gap is reported rather than
silently patched by docking now.

**Consequence**: the "ceiling" claim in Phase 5 ("the pocket cannot discriminate candidates from
decoys even in the experimentally determined holo crystal") was demonstrated on the apo structure,
not the holo one. The qualitative conclusion (no discrimination) likely still holds — 1EZM and
3DBK share the same catalytic pocket and Phase 6 showed no candidate/decoy-relevant Zn discipline
in either structure — but the specific "true holo ceiling" framing was not actually tested and
should not be cited as such going forward.

### Zinc/calcium presence — a second, more severe finding than Phase 6 reported

Phase 6 established that Vina models the Zn2+ ion with a generic Lennard-Jones potential, not
explicit coordination chemistry. Checking all receptor files used in this phase surfaces something
Phase 6 did not check: **the Boltz/ConforMix-generated conformer receptors (`lasb_conformer_beta*`)
contain no zinc, no calcium, and no HETATM records of any kind.**

| Receptor | Atom count | ZN present | CA present |
|---|---|---|---|
| 1EZM crystal (used in Phase 4/5/6/Control A) | 3220 | yes (charge +2.000) | yes (charge +2.000) |
| lasb_conformer_beta{0.0,0.8,1.6,2.4,3.2} (all 5) | 2865 | **no** | **no** |
| 3DBK holo (Phase 6) | — | yes | yes |

Traced to source: `results/lasb_payload/ensemble_receptors_aligned/lasb_conformer_beta0.0.pdb`
contains zero `HETATM` records at all — not stripped during pdbqt prep, but never present in the
ConforMix/Boltz-generated conformer structure in the first place.

**This means every "ensemble advantage" measured in Phases 4-5 (the entire premise motivating the
conformer-ensemble approach) was computed against a zinc metalloprotease pocket with the catalytic
zinc entirely absent in 5 of 6 receptors.** This is a stronger and more mechanistic explanation for
Control B's failure than the Phase 6 finding (Vina mistyping Zn) — the conformers do not
mis-model the zinc's contribution, they have no zinc to model at all. The uniform per-conformer
scoring uplift documented in Phase 5 (Table 1: all conformers ~1.4 kcal/mol better than crystal,
regardless of β) is consistent with a generic pocket-volume/softening effect from receptor
regeneration, unrelated to holo-geometry recognition, and is now additionally explained by the
absence of a buried, sterically demanding ion that the crystal structure has and the conformers
do not.

---

## Step 4 — Gate A Retrospective (known actives, N=11, pooled and thiol-stratified)

Vina values are carried over from Phase 6 for direct comparison. All values are `partial
r(pIC50, score | HAC)` as the primary metric per instructions — raw r is also reported.

| Receptor | Function | Pooled r | p | CI | Pooled ρ | r(score,HAC) | **partial r \| HAC** | thiol-only r (n=10) |
|---|---|---|---|---|---|---|---|---|
| 3DBK holo | Vina | −0.388 | 0.238 | [−0.78, 0.67] | −0.246 | −0.915 | **−0.088** | −0.063 |
| 3DBK holo | Vinardo | −0.433 | 0.184 | [−0.77, 0.54] | −0.223 | −0.858 | **−0.209** | −0.161 |
| 1EZM apo | Vina | −0.326 | 0.327 | [−0.77, 0.45] | −0.232 | −0.914 | **+0.078** | +0.015 |
| 1EZM apo | Vinardo | −0.007 | 0.983 | [−0.82, 0.80] | −0.041 | −0.499 | **+0.234** | +0.236 |

**No function or stratum clears the pre-registered |r| ≥ 0.4 threshold.** Vinardo does not recover
signal where Vina found none — if anything, Vinardo's raw pooled r on 3DBK (−0.433) looks nominally
larger than Vina's, but its partial correlation after HAC correction (−0.209) is *worse* (further
from zero in the wrong direction) than Vina's (−0.088), and its raw score is still driven
substantially by size (r(Vinardo,HAC) = −0.858 on 3DBK). **Per instructions: N=11 and 1.49 log-unit
range remain underpowered regardless of function; these results cannot rule out a real weak
correlation, and should not be read as confirming the null either.**

---

## Step 5 — Gate B Retrospective (known actives vs decoys — first screening-power test on this target)

**3DBK is infeasible as specified**: the protocol asks for actives-vs-decoys in 3DBK holo rigid,
but decoys were never docked into 3DBK (Step 2 above). Only 1EZM apo has both populations docked,
so Gate B is run there instead, with this substitution stated explicitly.

**HAC distributions (actives n=11, decoys n=20):** actives mean=21.5 (sd=5.3), decoys mean=19.6
(sd=5.2); t=0.965, p=0.342 — not significantly different, so a substantial fraction of any AUC
signal is not trivially attributable to a size gap between the two populations.

| Function | Metric | AUC | 95% CI | EF10% | EF20% |
|---|---|---|---|---|---|
| Vina | raw score | 0.541 | [0.321, 0.741] | 0.94 | 0.47 |
| Vina | ligand efficiency | 0.423 | [0.226, 0.638] | — | — |
| Vinardo | raw score | **0.727** | **[0.532, 0.894]** | 1.88 | 1.41 |
| Vinardo | ligand efficiency | 0.550 | [0.337, 0.763] | — | — |

Vina shows essentially chance discrimination (AUC ≈ 0.54). **Vinardo's raw AUC (0.727, CI entirely
above 0.5) is the first hint of any discrimination anywhere in this project** — but per the
explicit instruction ("a raw r improvement that vanishes after controlling for HAC is not an
improvement"), the ligand-efficiency-corrected Vinardo AUC collapses to 0.550 (CI spans 0.5), i.e.
**this apparent discrimination is substantially a size effect**, not evidence Vinardo recognizes
LasB-relevant binding character. Reported as-is, not dismissed: a raw AUC of 0.73 in a real
screening context would still enrich a virtual screen, even if the underlying reason is that actives
in this literature set happen to be somewhat larger/more polar than the decoy set on average
(HAC difference did not reach significance at n=11 vs n=20, but the LE-normalized collapse shows
the raw signal is not independent of size).

---

## Step 6 — Effect of Vinardo on Phase 4 / Control B Discrimination

| Metric | Vina (established) | Vinardo (this phase) |
|---|---|---|
| Noise-corrected delta, candidates (max-based) | +1.372 to +1.417 | +1.127 [CI 0.881, 1.345] |
| Decoy advantage (max-based) | +1.616 | +1.667 |
| **discrimination_gap (cand − decoy)** | **−0.244 (Phase5) / −0.198 (this segment)** | **−0.540 [CI −0.952, −0.156]** |
| Threshold (+0.30) | FALSIFIED | FALSIFIED — CI now entirely below zero |

**Does the ensemble show discrimination under Vinardo that it did not under Vina? No — it is
worse.** Vinardo's discrimination_gap CI no longer even touches zero (Vina's CI spanned zero); the
decoy advantage is confidently larger than the candidate advantage under Vinardo. This does not
reframe Control B's failure as scoring-function-limited — both scoring functions agree the
signal is non-specific, and the alternative function makes the case for non-specificity stronger,
not weaker.

Per-conformer means (Vinardo, kcal/mol):

| Receptor | Candidates | Decoys | Diff (cand−decoy) |
|---|---|---|---|
| 1EZM crystal | −3.634 | −3.710 | +0.076 |
| β0.0 | −4.370 | −5.027 | +0.656 |
| β0.8 | −4.324 | −4.758 | +0.434 |
| β1.6 | −4.154 | −4.504 | +0.350 |
| β2.4 | −4.055 | −4.392 | +0.337 |
| β3.2 | −4.380 | −4.995 | +0.615 |

Decoys beat candidates on every conformer under Vinardo too, by larger margins than under Vina.

---

## Step 7 — PDBbind Overlap Audit

The cloned `deltaVinaXGB` repo ships its actual training/validation/test CSVs
(`Dataset/Train_Val/Train_dry_pdbbind.csv`, 3264 PDBbind complexes; plus CSAR and validation
sets), so this could be checked directly against the real training index rather than inferred.

- **1EZM and 3DBK**: absent from every CSV in the repo's `Dataset/` directory (`grep -ril` across
  all Train/Val/Test files, case-insensitive, zero hits). Neither structure used in this project
  is a ΔvinaXGB training or test complex.
- **LasB/pseudolysin complexes in PDBbind generally**: none identified — LasB is not a registered
  ChEMBL target either (Phase 6), consistent with it being a niche/under-studied target for
  ML-scoring training data.
- **Thermolysin-family (M4 peptidase) complexes**: a spot-check of ~20 known thermolysin-associated
  PDB codes against the 3264-entry training list found **3 confirmed present** (3FED, 4TLN, 1ZDP);
  this is a lower bound from a small hand-picked list, not an exhaustive UniProt cross-reference,
  so the true count is almost certainly higher. **Family-level memorization is plausible** if
  ΔvinaXGB were ever run: thermolysin (M4) and LasB (M4) share the same catalytic fold and Zn-motif,
  so a model trained partly on thermolysin complexes could show apparent LasB signal that is really
  learned thermolysin-family geometry, not LasB-specific binding chemistry. This context matters if
  ΔvinaXGB is revisited in a different environment: any positive correlation found there should be
  discounted for this reason before being read as validating the scoring approach for LasB specifically.

---

## Step 8 — Free Geometric Analyses

### 8a. Interaction fingerprint recovery — INFEASIBLE

Requires comparing candidate/decoy pose fingerprints against the native RDF ligand's fingerprint in
3DBK. Blocked by the same gap as Step 2/5: no candidate or decoy poses exist in 3DBK, and 1EZM (the
receptor both populations were actually docked into) is apo — it has no native ligand to serve as
reference. This analysis cannot be run without new docking against 3DBK, which the rules for this
phase exclude.

### 8b. PoseBusters validity

Run on all 426 poses (0 conversion failures). One check, `no_radicals`, failed on **100% of
poses** — traced to the `obabel` pdbqt→SDF conversion step losing bond-order/formal-charge
information needed for correct radical-electron perception; this is a conversion artifact common
to all poses equally, not a real chemistry finding, and is excluded from the pass-rate below.

**Pass rate on all remaining checks** (sanitization, InChI-convertibility, connectivity, bond
lengths/angles, internal clash, ring/double-bond flatness, internal energy, protein-ligand
distance and clash, volume overlap):

| Population | Pass rate | n |
|---|---|---|
| Phase 4 ensemble (candidates) | 88.3% | 120 |
| Control A (candidates, 6 seeds) | 90.0% | 120 |
| Control B (decoys) | 93.3% | 120 |
| Phase 6 (known actives) | 98.5% | 66 |

**Crystal vs conformer receptor, directly answering the protocol's question** ("if conformer poses
fail validity checks more often than crystal poses, the ensemble advantage was partly rewarding
physically implausible geometry"):

| Receptor type | Pass rate | n |
|---|---|---|
| Crystal (1EZM) | 90.0% | 40 |
| Conformer (β0.0-3.2) | 91.0% | 200 |

**Conformer and crystal pass rates are statistically indistinguishable (91.0% vs 90.0%).** The
ensemble's scoring advantage over crystal (Phase 4/5) is not explained by conformer receptors
systematically accepting less physically valid poses — whatever drives the uplift, it isn't cruder
pose validity. All observed failures (besides the universal radical artifact) cluster by ligand
identity, not by receptor: candidates 0008 and 0012 fail identically across all 6 receptors each,
and decoy 0003 fails identically across all 6 receptors — consistent with those specific molecules
having a structural feature (e.g. an unusual functional group) that trips the InChI-convertibility
check regardless of which receptor they're docked into, not a receptor-dependent effect.

### 8c. Torsional strain (MMFF94, docked − freely minimized)

Computed for candidates (Phase 4 ensemble, 120 poses) and decoys (Control B, 120 poses).

| Population | Crystal mean (kcal/mol) | Conformer mean (kcal/mol) |
|---|---|---|
| Candidates | 204.4 | 204.8 |
| Decoys | 147.2 | 146.6 |

Candidates show significantly higher computed strain than decoys overall (t=6.60, p<0.001), and
strain does not differ meaningfully between crystal and conformer receptors for either population
(t=−0.02, p=0.99) — i.e. whatever separates candidates from decoys here is a ligand-intrinsic
property, not a receptor-geometry effect.

**Caveat, stated prominently**: absolute strain values in the 150-400 kcal/mol range are not
chemically plausible for reasonable ligand conformations (typical torsional/pose strain for
drug-like molecules is single- to low-double-digit kcal/mol). This reflects an artifact of the
`obabel`-based pdbqt→SDF conversion plus `RDKit AddHs(addCoords=True)` protocol used here — almost
certainly poor hydrogen placement/valence perception on the converted poses (RDKit logged repeated
atropisomer-geometry warnings during processing, consistent with imperfect 3D perception) — not a
real finding of extreme strain. **The relative candidates > decoys direction is reported as a data
point; the absolute magnitudes should not be cited or compared to literature strain-energy
benchmarks.**

### 8d. Pose reproducibility (Control A, 6 seeds)

Mean pairwise RMSD among top poses per candidate across the 6 Control A seeds: mean=3.18 Å,
sd=3.34 Å, range 0.27–11.11 Å across the 20 candidates — high variance, several candidates show
essentially no reproducibility (>10 Å between seeds), consistent with a search-space or scoring
landscape without a clearly dominant pose for at least some ligands.

Exploratory correlation between reproducibility (mean pairwise RMSD) and Vina crystal score:
r=0.019, p=0.938 — no relationship. (No pIC50 is available for the DiffSBDD candidates, so this
correlation could not be run against affinity as specified; it was run against docking score
instead, as the closest available proxy, and is explicitly exploratory per instructions.)

---

## Summary and Verdict

1. **Step 1**: HAC is balanced between candidates/decoys; the Control B decoy advantage is not an
   HAC artifact and persists in ligand-efficiency units.
2. **Step 2**: 426/466 requested poses exist with coordinates; the 40-pose 3DBK candidate/decoy
   set does not exist — Phase 5's "Control C true holo ceiling" was actually run on 1EZM (apo),
   not 3DBK (holo), and should be re-labeled accordingly. **New finding**: all 5 conformer
   receptors used throughout Phases 4-5 contain no zinc, calcium, or any HETATM at all — a more
   severe and more direct explanation for the non-specific "ensemble advantage" than Vina's
   generic Zn-typing (Phase 6).
3. **Step 0/3**: ΔvinaXGB/ΔvinaRF20 could not be installed after three distinct, diagnosed
   dependency-version dead ends; genuinely infeasible in this environment, not a null result.
   Vinardo substituted successfully via the native Vina binary, no install needed.
4. **Step 4**: Vinardo does not recover affinity signal in the known-actives set (still fails
   |r|≥0.4 pooled and thiol-stratified, both receptors); underpowered dataset caveat applies to
   both functions equally.
5. **Step 5**: First real actives-vs-decoys screening test (1EZM, 3DBK infeasible). Vina: chance
   AUC (0.54). Vinardo: raw AUC 0.73 (CI above chance) but collapses to 0.55 after ligand-efficiency
   correction — apparent discrimination is a size effect, not validated screening power.
6. **Step 6**: Vinardo does not rescue Control B — the discrimination_gap gets more negative and
   more confidently non-zero, not less.
7. **Step 7**: 1EZM/3DBK confirmed absent from ΔvinaXGB's actual training set (checked directly
   against the shipped training CSV); thermolysin-family (M4) complexes confirmed present in that
   set, so any future ΔvinaXGB signal on LasB should be discounted for family-level memorization risk.
8. **Step 8**: Interaction-fingerprint recovery blocked by the same 3DBK pose gap as Step 5.
   PoseBusters validity is statistically identical for conformer (91.0%) and crystal (90.0%)
   receptors — the ensemble advantage is not explained by conformers accepting less physically
   valid poses. Torsional strain shows a candidates>decoys directional difference but absolute
   values are not physically interpretable (conversion artifact, stated explicitly). Pose
   reproducibility is highly variable across candidates and uncorrelated with docking score.

**Overall verdict**: nothing in this retrospective rescoring changes Phase 4-6's conclusions.
Vinardo, tested as the one alternative scoring function that could actually be run in this
environment, does not find candidate-specific signal anywhere Vina failed to — if anything it
makes the non-discrimination finding more confident (Step 6) — and the one place it showed
above-chance raw discrimination (Step 5 AUC) evaporates under a ligand-efficiency size correction.
The strongest new finding of this phase is not about scoring functions at all: the Boltz/ConforMix
conformer receptors used throughout Phases 4-5 never contained the catalytic zinc in the first
place, which is a more direct and more severe explanation for the ensemble's non-specific scoring
uplift than anything about how Vina treats zinc once present.

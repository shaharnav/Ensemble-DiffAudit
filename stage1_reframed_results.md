# Stage 1 (reframed) + Stage 2: Four-Model Comparison and Decorrelation

Supersedes the reproduction-gate framing in `stage1_results.md` for purposes of deciding what to
report going forward — that file's diagnostic investigation stands and is not repeated here. This
document follows through on it: the ≥0.70 gate was calibrated against ΔvinaRF20 trained on ~4,000
PDBbind complexes and evaluated on 285 held-out ones. At n≈230 under leave-one-family-out CV, a
model landing well below 0.70 mostly measures the tenfold-plus difference in training data, not
whether the reimplementation works. Reporting that as a "failed reproduction" would overstate what
was learned. The replacement bar: **beat the Vina baseline on the identical held-out split, with
bootstrap CIs, so "beat" is judged against real uncertainty rather than a point estimate.**

---

## Confirming the two root causes (cheap checks, run before choosing a path)

**Root cause #1 — cancellation, confirmed directly.** `r(predicted_residual, −vina_score) = 0.846`
(p=1.8e-66). The delta model's out-of-fold residual predictions are almost the same thing as
`−vina_score`. Adding them back together doesn't correct Vina's score — it partially erases it.
This is now a clean, quotable finding rather than an inference from indirect symptoms.

**Root cause #2 — is B-F signal independent of Vina, or just "trust Vina, adjust slightly"?**
Ran direct pKd regression on B-F features **with vina_score excluded entirely**:

| Model | Scoring power (Pearson R) |
|---|---|
| B-F features only, no vina_score | **0.598-0.602** |
| Vina baseline | 0.581 |

**B-F alone, with no access to Vina's own score at all, already edges out the Vina baseline.** The
signal in the buried-SASA/contacts/H-bond-satisfaction/ligand-descriptor features is genuinely
independent, not a repackaging of Vina's own prediction. Adding `vina_score` back as one more input
feature gives a further, small bump (0.618) — consistent with Vina carrying some additional
information on top of B-F, not with the model's win being pure pass-through.

---

## Four-model comparison, identical LOGO-CV split (n=238, 71 UniProt families)

| Model | Scoring power (Pearson R) | 95% CI | Ranking power (mean within-family ρ) | r(score, HAC) |
|---|---|---|---|---|
| Vina baseline | 0.581 | [0.496, 0.657] | 0.498 | 0.593 |
| Delta model (residual target, all 44 features) | 0.435 | [0.328, 0.534] | 0.406 | 0.622 |
| **Direct regression (B-F + vina_score)** | **0.618** | **[0.533, 0.694]** | **0.596** | **0.766** |
| Direct regression (B-F only, no vina_score) | 0.598 | [0.512, 0.676] | 0.568 | 0.785 |

**The honest gate — beat Vina baseline — is nominally cleared (0.618 > 0.581) but the 95% CIs
overlap substantially** (Vina's upper bound 0.657 sits inside the direct model's CI; the direct
model's lower bound 0.533 sits inside Vina's CI). **At n=238 under leave-one-family-out CV, the
improvement from 0.581 to 0.618 is not statistically distinguishable from noise.** Ranking power
shows a numerically larger relative gap (0.498 → 0.596) but the same small-N caveat applies; no
bootstrap CI was computed for ranking power in this pass (would require resampling within each of
47 families separately — flagged as a gap, not computed here given time-box).

**A complication worth stating plainly, not burying**: the model that nominally beats Vina also has
a *higher* r(score, HAC) than Vina itself (0.766 vs 0.593). Numerically "improving on Vina" here
comes bundled with leaning more heavily on molecular size, not less — the opposite of what a
size-decorrelated scoring function should look like. This is exactly why Stage 2 (below) is run
before drawing any conclusion about whether this is a real improvement worth using.

---

## Stage 2 — Size decorrelation (run on the direct model, not the failed delta model)

Both variants tested on the identical LOGO-CV split, against the direct-regression baseline
(B-F + vina_score, R=0.618, r(HAC)=0.766):

| Model | Scoring power | 95% CI | r(score, HAC) |
|---|---|---|---|
| Baseline (direct regression, undecorrelated) | 0.618 | [0.533, 0.694] | 0.766 |
| (i) Features residualized against HAC | 0.513 | [0.417, 0.601] | 0.539 |
| (ii) Target = ligand efficiency (pKd/HAC), converted back to pKd scale | 0.497 | [0.360, 0.629] | 0.599 |

**Both decorrelated variants reduce r(score,HAC) substantially (0.766 → 0.54-0.60) but at a real
cost to scoring power — both fall below even the plain Vina baseline (0.581).** Reported plainly,
per instructions: a size-invariant function that is honest about underperforming is a legitimate
result, not a failure to hide.

Neither decorrelation approach achieves full independence from HAC (residual r(HAC) is still
~0.54-0.60, not ~0). This is consistent with per-feature linear residualization only partially
removing a real, likely nonlinear and multivariate size dependency — a known limitation of this
specific decorrelation method, not evidence the underlying features have no non-size information.

### Net picture across all six configurations

| | Scoring power | r(score,HAC) |
|---|---|---|
| Vina baseline | 0.581 | 0.593 |
| Delta model (failed) | 0.435 | 0.622 |
| Direct model (best raw performance) | 0.618 | 0.766 |
| Direct model, features decorrelated | 0.513 | 0.539 |
| Direct model, LE target | 0.497 | 0.599 |

No configuration simultaneously beats Vina's scoring power *and* reduces its size dependence. The
best-performing model (direct regression) achieves its (statistically inconclusive) edge over Vina
by being *more* size-driven, not less; the decorrelated variants trade away scoring power for a
genuine reduction in size dependence but land below the plain Vina baseline. **This is itself the
headline finding of the bias-audit exercise**: on this feature set and this sample size, there is
no free lunch between "beats Vina" and "less size-biased than Vina" — you can pick one axis to
improve, not both, which is precisely the tension the whole delta-learning-reimplementation project
set out to probe.

---

## Data-quality note: the 6 dropped complexes (not silently dropped)

`3u8k`, `3u8n`, `1g2k`, `2w66`, `2vw5`, `2xys` — all 6 produced `NaN` in the buried-SASA feature
group and were excluded from modeling (244 → 238). Root cause: their ligand pdbqt files contain
`CG0`/`G0` atom types — Meeko's "glue pseudo-atoms," emitted when a macrocyclic ring is opened for
flexible-macrocycle docking prep. These aren't real elements and this pipeline's SASA path doesn't
handle them (falls back to a default radius that produces a downstream NaN in FreeSASA's per-atom
area calculation for at least one atom in the complex). This affects only macrocyclic ligands
specifically, so the exclusion is not random with respect to ligand chemotype — a real, stated
limitation of this n=238 result, not a random 2.5% attrition.

---

## Deliverables

- `results/casf2016/stage1_four_model_comparison.csv` — the four-model table above
- `results/casf2016/stage2_decorrelation.csv` — the Stage 2 table above
- `casf_stage1_reframed.py` — reproducible script for both tables

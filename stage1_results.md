# Stage 0.5 + Stage 1: Feature Profiling, Delta Model, Reproduction Gate

## Stage 0.5 — Profiling (time box: 15 min)

Measured on 50 complexes (in-process, no subprocess spawning — features computed directly from
parsed pdbqt coordinate arrays via numpy/scipy):

| Group | ms/complex | % of total |
|---|---|---|
| pdbqt parsing | 78.5 | 84.8% |
| A: Vina terms | 1.3 | 1.4% |
| B: buried SASA (FreeSASA) | 10.3 | 11.1% |
| C: contacts | 1.2 | 1.3% |
| D/E: H-bonds + satisfaction | 0.25 | 0.3% |
| F: ligand descriptors | 1.05 | 1.1% |
| **Total** | **92.5 ms/complex** | |

Projected time for full 244-complex set: **22.6s**. Gate (must be <2hr): **PASS**, by more than
three orders of magnitude — no scaling concerns.

Actual full-set run: 244/244 featurized in **3.5-4.1s** (parallelized, 8 workers). 6 complexes
produced NaN in the buried-SASA group (traced to meeko macrocycle "glue" pseudo-atoms — `CG0`/`G0`
atom types from ring-opened flexible macrocycles — not a real element, not handled by this
pipeline's SASA path) and were dropped. **n = 238** for modeling.

Feature-vs-HAC correlation table: `results/casf2016/feature_hac_correlation.csv`. Strongest
size-correlated features: MW (r=0.98, expected — it's essentially HAC), `contact_C_6_8` (r=0.84),
`vina_gauss2` (r=0.83), `vina_hydrophobic` (r=0.74), `buried_sasa_C` (r=0.74). Weakest: `vina_hbond`
(r=0.01), most `buried_sasa_{N,S}` (r<0.05), most polar-contact bins. As expected, size correlation
concentrates in carbon-contact and hydrophobic-term features; polar/H-bond features are close to
size-independent by construction.

---

## Stage 1 — Delta Model

**Split**: no PDBbind refined set (registration wall) — leave-one-UniProt-family-out CV within the
238-complex CASF-2016-derived set itself, per the Stage -1 fallback. 71 UniProt groups (via
`rcsb_polymer_entity_container_identifiers.uniprot_ids`, one GraphQL query per complex). This is
a **materially different and smaller training regime** than the published papers, which train on
~1300-5000 PDBbind complexes and evaluate on a fully separate CASF-2016 holdout. Flagged as small-N
throughout, as instructed.

**Target**: `pKd_experimental - vina_score` (delta learning, as specified).

**Reproduction gate, stated before running**: published ΔvinaRF20 (Wang & Zhang 2017) CASF-2016
scoring power ≈ 0.803; ΔvinaXGB (Lu et al. 2019) ≈ 0.82-0.86. Gate: Pearson R ≥ 0.70 (within ~0.1
of 0.80). Caveat stated up front: with ~230 training examples instead of 1300-5000, and an
approximate from-scratch feature set (geometric Vina terms, geometric H-bonds, FreeSASA not MSMS),
landing at 0.70 was flagged as optimistic before running — a lower number would not by itself prove
broken features, given these confounds.

### Result: GATE FAILS, and the delta model underperforms the raw Vina baseline

| Model | Scoring power (Pearson R) | Ranking power (mean within-family ρ, 47 families n≥3) |
|---|---|---|
| Vina baseline | **0.581** | 0.498 |
| Vinardo baseline | 0.513 | 0.434 |
| Delta model (all 44 features, default hyperparams) | 0.435 | 0.406 |

The delta model does not beat Vina, let alone approach the 0.70 gate.

### Debugging (per the hard rule: do not report failure without investigating)

**Diagnostic 1 — is the model learning anything?** Checked correlation between predicted residual
and true residual directly (before adding back `vina_score`): **r = 0.77-0.80** across several
hyperparameter settings. The model clearly learns real structure in the residual target
out-of-fold. But scoring power on the reconstructed pKd (`pred_residual + vina_score`) collapsed to
0.25, or even -0.005 for a heavily regularized configuration — reconstruction destroys the signal
the model demonstrably has.

**Diagnostic 2 — collinearity with the base score.** Group A features (the 5 Vina terms) are, by
construction, the raw components `vina_score` is a weighted sum of. Checked directly:
`r(vina_hydrophobic, vina_score) = -0.785`, `r(vina_gauss2, vina_score) = -0.719`,
`r(vina_gauss1, vina_score) = -0.683`. A model with access to these can trivially reconstruct
something close to `-vina_score`; added back to `vina_score`, this cancels the informative part of
the base score and leaves noise. **Root cause candidate #1: Group A features are collinear with the
quantity being corrected.**

Re-ran LOGO-CV excluding Group A (39 features, B/C/D/E/F only):

| Model | Scoring power |
|---|---|
| XGBoost, original hyperparams | 0.421 |
| XGBoost, regularized (fewer trees, shallower, higher reg_lambda) | 0.120 |
| XGBoost, very shallow (max_depth=2, n_estimators=30) | -0.040 |
| Ridge (α=1) | 0.419 |
| Ridge (α=10) | 0.412 |
| Ridge (α=50) | 0.370 |
| Ridge (α=200) | 0.249 |

Removing Group A helps (0.421-0.435 range vs Vina's 0.581) but does **not** close the gap, and more
regularization makes it *worse*, not better — the opposite of what overfitting-driven collinearity
alone would predict, and a sign the residual target is close to unlearnable at this N regardless of
model complexity.

**Diagnostic 3 — is it the delta-target formulation, not the features?** Ran the same B-F features
(plus `vina_score` itself as one plain input feature, not an additive baseline to subtract from) in
a **direct pKd regression** instead of residual regression:

| Model | Scoring power |
|---|---|
| Ridge (α=10), direct pKd, vina_score as a feature | 0.597 |
| Ridge (α=50), direct pKd, vina_score as a feature | **0.618** |
| Ridge (α=200), direct pKd, vina_score as a feature | 0.623 |

**This modestly beats the Vina baseline (0.581).** The same underlying features carry real,
usable signal — the problem is specific to the residual-regression (`pKd - vina_score`) target
construction combined with small N and leave-one-family-out CV, not the feature extraction itself.

### Conclusion

**The reproduction gate as specified (delta/residual model, ≥0.70) FAILS: 0.435.** This is reported
plainly, per instructions, rather than silently substituting the direct-regression formulation that
happens to look better — that would be moving the goalposts after seeing results.

The debugging investigation identifies the mechanism, though: delta learning's residual target is
fragile at ~230 training examples under leave-one-family-out CV, and features collinear with the
base score actively hurt it via reconstruction cancellation. A closely related formulation (direct
regression with the base score as one input among several) does carry positive signal over the
Vina baseline on the identical feature set. This localizes the failure to the training-data-scale /
target-construction interaction that Stage -1's fallback plan explicitly anticipated ("small-N
throughout"), not to a bug in SASA, contacts, H-bond geometry, or ligand descriptors individually.

**Per the hard rule ("do not proceed on unvalidated features"), Stage 2 (size decorrelation) is not
run.** The gate is not met, and proceeding to build decorrelated variants on top of a model that
already fails its own reproduction check would compound an unvalidated result rather than test a
validated one.

### What would plausibly fix this (not attempted — out of current scope)

- Direct regression (score = f(features, vina_score)) instead of pure residual regression, which
  Diagnostic 3 shows already helps on identical features/data.
- More training data — the core, structural limitation. This is exactly what the PDBbind
  registration wall prevents in this environment.
- A published protein-family clustering (real Pfam, not UniProt-only) to avoid any residual
  leakage between near-identical family assignments across complexes.

---

## Deliverables

- `results/casf2016/stage1_features.csv` — 244×47 raw feature matrix (238 used after NaN drop)
- `results/casf2016/feature_hac_correlation.csv` — feature-vs-HAC correlation table
- `results/casf2016/stage1_oof_predictions.csv` — out-of-fold predictions from the (failed-gate)
  delta model, for transparency
- `casf_pipeline/features.py`, `pdbqt_atoms.py` — reusable feature extraction package
- `casf_stage05_profile.py`, `casf_featurize_all.py`, `casf_stage1_train.py` — pipeline scripts

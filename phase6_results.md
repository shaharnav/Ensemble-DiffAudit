# Phase 6 Results: Affinity Signal Calibration in LasB

Pre-registration: `phase6_preregistration.md`, committed `d9b094a` before docking.

---

## Step 0 — Receptor Prep Audit

**Catalytic Zn2+ (residue 302):**
Present in all prepared pdbqt files. Meeko mk_prepare_receptor converts the HETATM record
to an ATOM record, assigns atom type `Zn`, charge +2.000. The zinc is structurally retained.

**Structural Ca2+ (residue 400):**
Present. Atom type `Ca`, charge +2.000.

**What Vina does with Zn:**
AutoDock Vina uses no explicit metal coordination term. The Zn is treated as a generic
heavy atom with a Lennard-Jones 4-8 potential. Metal coordination geometry, ligand-field
effects, and the entropic cost of displacing the catalytic water are entirely unmodeled.
Coordination energy for tight-binding ZBGs (hydroxamate: ΔG ~5 kcal/mol; thiol: ~3 kcal/mol)
is approximated by generic Lennard-Jones overlap — a known inadequacy for metalloprotease targets.

**Implication for all prior phases (Phases 4–5):**
Every docking score in this project was computed against a pocket where the primary
affinity driver (Zn coordination) is unmodeled. Phase 5's Control B failure (decoys matching
candidates) and Phase 6's correlation result (below) both follow mechanistically from this.

---

## Step 1 — Known-Active Dataset

ChEMBL does not register LasB (UniProt P14901) as a target. Dataset assembled from primary
literature. Gate failures were pre-registered before docking; the experiment was labeled
underpowered in the preregistration file.

**Final dataset:** n = 11 (10 thiol mercaptoacetamides + 1 phosphoramidate)

| ID | pIC50 | ZBG | MW | HAC | IC50 (μM) | Source |
|----|-------|-----|----|-----|-----------|--------|
| Engel_3 | 6.319 | thiol | 271 | 19 | 0.48 | PMC9112332 |
| Engel_5 | 5.921 | thiol | 257 | 18 | 1.20 | PMC9112332 |
| Engel_11 | 6.000 | thiol | 302 | 21 | 1.00 | PMC9112332 |
| Engel_12 | 6.155 | thiol | 287 | 20 | 0.70 | PMC9112332 |
| Engel_13 | 6.222 | thiol | 273 | 19 | 0.60 | PMC9112332 |
| Engel_15 | 5.131 | thiol | 289 | 20 | 7.40 | PMC9112332 |
| Engel_16 | 5.602 | thiol | 318 | 22 | 2.50 | PMC9112332 |
| Engel_17 | 5.959 | thiol | 291 | 19 | 1.10 | PMC9112332 |
| Engel_19 | 6.620 | thiol | 321 | 21 | 0.24 | PMC9112332 |
| Engel_23 | 6.222 | thiol | 314 | 21 | 0.60 | PMC9112332 |
| Phosphoramidon | 6.602 | phosphoramidate | 543 | 37 | 0.25 (Ki) | PMC9404851/Morihara1978 |

---

## Step 2 — Gate Check Results

**a) N ≥ 15?** NO (n = 11). Gate FAILS.
LasB is not in ChEMBL; key potency data (hydroxamate Ki = 0.002 μM, carboxylate Ki = 0.16/14.8 μM)
from ACS/Elsevier paywalls were inaccessible. Whole-cell GFP-reporter assays excluded per protocol.

**b) pIC50 range ≥ 2.0 log units?** NO (range = 5.131–6.620 = 1.489 log units). Gate FAILS.
The missing hydroxamate compound would have added 2+ log units. This is the binding data gap
that limits interpretation.

**c) ZBG class distribution:** 1 class with n ≥ 5 (thiol, n = 10). Phosphoramidate n = 1.
No hydroxamate or carboxylate with verified SMILES.

**d) Cross-source consistency:** No overlap between sources; cannot estimate inter-lab noise floor.
Intra-source SD from PMC9112332: ~5–15% of IC50 (CV ~10%), ≈ 0.05 log units per compound.

**Both primary gates fail. Experiment is underpowered. All results below should be interpreted
with this caveat.**

---

## Step 4 — Docking Results

Receptors: 3DBK (holo, Zn/Ca/HOH retained, RDF+SO4 stripped) and 1EZM apo. 3 seeds per compound,
mean best affinity reported. Box: 24×24×24 Å, exhaustiveness 16. Same prep pipeline as Phases 4–5.

**3DBK box center**: (18.721, −5.093, 23.685) — RDF ligand centroid in 3DBK native coordinate frame.
**1EZM box center**: (55.521, 35.882, 20.807) — same as Phases 4–5.

| ID | pIC50 | 3DBK (kcal/mol) | 1EZM (kcal/mol) |
|----|-------|-----------------|-----------------|
| Engel_19 | 6.620 | −7.514 | −5.675 |
| Phosphoramidon | 6.602 | −9.850 | −6.817 |
| Engel_3 | 6.319 | −7.780 | −5.443 |
| Engel_12 | 6.155 | −7.890 | −5.477 |
| Engel_13 | 6.222 | −7.311 | −5.319 |
| Engel_23 | 6.222 | −7.695 | −5.967 |
| Engel_11 | 6.000 | −8.329 | −5.727 |
| Engel_17 | 5.959 | −7.664 | −5.361 |
| Engel_5 | 5.921 | −7.432 | −5.166 |
| Engel_16 | 5.602 | −8.281 | −6.061 |
| Engel_15 | 5.131 | −7.186 | −5.460 |

---

## Step 5 — Correlation Analysis

### 5.1 Pooled Pearson / Spearman (all 11 compounds)

| Receptor | Pearson r | p | 95% CI | Spearman ρ | p |
|----------|----------|---|--------|------------|---|
| 3DBK holo | −0.389 | 0.238 | [−0.778, +0.673] | −0.246 | 0.466 |
| 1EZM apo | −0.326 | 0.327 | [−0.761, +0.445] | −0.232 | 0.492 |

Neither reaches the pre-registered threshold of |r| ≥ 0.40. Both CIs span zero. No significant
correlation in either direction.

### 5.2 Stratified: Thiol class only (n = 10)

| Receptor | Pearson r | p | 95% CI | Spearman ρ | p |
|----------|----------|---|--------|------------|---|
| 3DBK holo | −0.063 | 0.862 | [−0.712, +0.786] | −0.079 | 0.828 |
| 1EZM apo | +0.015 | 0.967 | [−0.568, +0.646] | −0.024 | 0.947 |

Within the thiol class — where Zn-coordination strength is approximately constant and should
cancel — Vina scores carry essentially zero affinity information. r ≈ 0 for both receptors.

This is the most informative result: even after removing the between-ZBG Zn-coordination
variation, Vina cannot rank compounds by tail contacts.

### 5.3 Size Controls

| Metric | 3DBK holo | 1EZM apo |
|--------|-----------|----------|
| r(Vina, HAC) | **−0.915** | **−0.914** |
| r(pIC50, HAC) | +0.389 | +0.389 |
| partial r(pIC50, Vina \| HAC) | −0.088 | +0.078 |

Vina scores are almost entirely a function of heavy atom count (r = −0.91). After controlling
for HAC, the partial correlation between pIC50 and Vina score collapses to ≈ 0 (−0.09 and +0.08
for 3DBK and 1EZM respectively). Any apparent pooled correlation (r = −0.39) is entirely
explained by size: phosphoramidon (HAC = 37, the largest compound) happens to be the highest-
pIC50 compound AND scores best — but this is coincidence, not discrimination.

### 5.4 Mechanistic Interpretation

The Zn coordination energy is the dominant affinity driver for all LasB inhibitors. It is entirely
unmodeled by Vina. What Vina does model (Lennard-Jones contact area, H-bonds, rotatable bond
desolvation) is correlated with molecular size — hence r(Vina, HAC) = 0.91 — but not with the
potency variation within a congeneric series sharing the same ZBG.

---

## Step 6 — Rescoring

**Infeasible as stated in protocol.** gnina and vinardo are not installed in this environment.
Could not rescore poses.

If gnina were available, the expectation based on Step 5 results is that CNNaffinity (which
models coordination geometry via learned features) would show better correlation than Vina,
particularly across ZBG classes. This is the key test that cannot be run without gnina.

---

## Step 7 — Calibration

**Threshold not met for any function or stratum.** No calibration performed.
|r| = 0.389 (pooled 3DBK, closest to threshold) does not clear 0.40, collapses to 0.09 after
HAC correction, and carries CI [−0.78, +0.67] that is consistent with the true r being zero.

No correction was applied to Phase 4–5 candidate scores. The DiffSBDD candidate pIC50
predictions would require a valid calibration to be meaningful, and none exists.

---

## Verdict

**NULL: Vina scores in the LasB pocket carry no affinity information.**

| Metric | 3DBK holo | 1EZM apo | Pre-reg threshold | Status |
|--------|-----------|----------|-------------------|--------|
| Pooled r | −0.389 | −0.326 | ≥ 0.40 | FAILS |
| Thiol r | −0.063 | +0.015 | ≥ 0.40 | FAILS |
| Partial r (HAC-corrected) | −0.088 | +0.078 | — | ≈ 0 |

This null result is mechanistically explained, not merely statistical noise:
1. LasB is a zinc metalloprotease where Zn-coordination is the dominant affinity driver.
2. Vina lacks an explicit Zn-coordination term.
3. Within a congeneric thiol series (constant ZBG, n = 10), r ≈ 0 — confirming the scoring
   function cannot rank tail-group contacts either.
4. Vina scores are HAC proxies (r = 0.91), not affinity proxies.

The null result explains the entire Phase 5 failure chain:
- Phase 4: candidates score better than crystal → real signal, but driven by size/softening, not selectivity
- Control B: decoys match candidates → scoring offset applies to any molecule, confirmed null discrimination
- Control C: holo crystal cannot discriminate → ceiling problem confirmed
- Phase 6: known actives show no pIC50–score correlation → the ceiling problem has a mechanistic root in unmodeled Zn coordination

---

## Caveats

1. **Underpowered**: n = 11, range = 1.49 log units. Cannot rule out weak real correlation.
2. **Missing hydroxamates**: The most potent class (Ki 0.002–17.4 μM, pIC50 4.8–8.7) was absent.
   If included, a pooled r driven by ZBG-class scoring differences might appear and could be
   mistaken for affinity signal.
3. **SMILES reconstruction**: Thiol compound SMILES were reconstructed from structural descriptions
   (PMC9112332 Table 1–2). Structures of heteroaryl derivatives (Engel_19, Engel_23) are based
   on the stated substitution pattern; minor regiochemical uncertainty exists.
4. **Step 6 not run**: gnina/CNNaffinity rescoring could not be performed.

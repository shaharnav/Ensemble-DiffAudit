# Phase 6 Pre-registration: Affinity Signal Calibration in LasB

**Committed before any Phase 6 docking.**

---

## Background

Phase 5 showed that AutoDock Vina scores in the LasB (1EZM/3DBK) pocket cannot distinguish
DiffSBDD-generated candidates from property-matched decoys — not in the ConforMix conformer
ensemble (Control B, discrimination_gap = −0.244 kcal/mol, FALSIFIED) and not in the true
holo crystal (Control C, gap = +0.261 kcal/mol, CI [−0.182, +0.686]).

Phase 6 asks: does Vina carry ANY affinity information in this pocket at all? We test this by
correlating Vina scores against published experimental IC50/Ki values for known LasB inhibitors.

---

## Step 1 — Dataset Assembled (pre-docking)

**Source**: PMC9112332 (Engel et al., ACS Infect Dis 2022) — α-substituted mercaptoacetamides.
Assay: purified LasB enzyme, FRET fluorogenic substrate (BODIPYfl-DQ elastin or Mca-PLGL-Dpa-AR-NH2).
All values are exact IC50 ± SD from Table 1 (aryl series) and Table 2 (heteroaryl, exact values only).

Additionally: Phosphoramidon (PubChem CID 445114), Ki = 0.25 μM from Morihara & Tsuzuki 1978,
cited by PMC9404851 review. Assay: purified LasB, substrate hydrolysis.

ChEMBL status: LasB (pseudolysin, UniProt P14901) is NOT registered as a ChEMBL target.
No ChEMBL bioactivity data available. Dataset assembled entirely from primary literature.

**Final dataset after all filters:**

| ID | pIC50 | ZBG class | Source |
|----|-------|-----------|--------|
| Engel_3 | 6.319 | thiol | PMC9112332 |
| Engel_5 | 5.921 | thiol | PMC9112332 |
| Engel_11 | 6.000 | thiol | PMC9112332 |
| Engel_12 | 6.155 | thiol | PMC9112332 |
| Engel_13 | 6.222 | thiol | PMC9112332 |
| Engel_15 | 5.131 | thiol | PMC9112332 |
| Engel_16 | 5.602 | thiol | PMC9112332 |
| Engel_17 | 5.959 | thiol | PMC9112332 |
| Engel_19 | 6.620 | thiol | PMC9112332 |
| Engel_23 | 6.222 | thiol | PMC9112332 |
| Phosphoramidon | 6.602 | phosphoramidate | PMC9404851/Morihara1978 |

---

## Step 2 — Gate Checks (pre-docking, report here)

**a) N ≥ 15?**
NO. N = 11. Gate FAILS.
- Hydroxamate and carboxylate compounds excluded: SMILES unavailable for structures from behind
  journal paywalls; hydroxamate compound with Ki = 0.002 μM could not be verified.
- Whole-cell assay compounds (LasB-IN-1, IC50 = 8.7 μM via GFP reporter) excluded per protocol.
- Implication: the experiment is underpowered. Proceed per protocol but report this limitation.

**b) pIC50 dynamic range ≥ 2 log units?**
NO. Range = 5.131–6.620 = 1.489 log units. Gate FAILS.
- Missing hydroxamate (pIC50 ≈ 8.7) would add 2.1 log units; carboxylates (pIC50 ≈ 4.8) would add
  0.3 log units lower. Inaccessibility of full paper data is the limiting factor.
- **Critical implication**: with range < 1.5, even a real underlying r of 0.6 would not be
  distinguishable from r = 0 at n = 11 (90% CI width ≈ ±0.50). The experiment is substantially
  underpowered. A null result here cannot rule out a true weak correlation.

**c) ZBG class distribution:**
- Thiol: n = 10 (the only class with n ≥ 5)
- Phosphoramidate: n = 1
- Hydroxamate: n = 0 (structure unavailable)
- Carboxylate: n = 0 (structure unavailable)
Only one class has n ≥ 5, so within-class stratification has limited informativeness.

**d) Cross-source consistency:**
No compound appears in more than one source. Cannot estimate experimental noise floor from
cross-source agreement. From the ACS Infect Dis 2022 paper, typical assay SD ≈ 5-15% of IC50
(coefficient of variation ~10%), corresponding to ≈ 0.05 log units noise floor per compound.

---

## Step 3 — Pre-registered Metrics

**Primary**: Pearson r (Vina score vs pIC50), pooled, for 3DBK holo and 1EZM apo separately.

**Secondary**: Spearman rho (rank correlation, robust to non-linearity).

**Stratified**: Same, within thiol class only (n = 10). Phosphoramidate excluded (n = 1).
Rationale: within a congeneric series sharing a thiol ZBG, the unmodeled Zn-coordination
contribution is approximately constant and partially cancels, leaving tail-group contacts
that Vina can in principle score.

**Size control**: Pearson r(Vina score, heavy atom count) and r(pIC50, heavy atom count).
Partial r(pIC50, Vina score | HAC) tests whether any pooled correlation survives HAC correction.

**Rescoring**: Vinardo and gnina (CNNaffinity, CNNscore) applied to same poses.

---

## Decision Thresholds (pre-registered before docking)

| Condition | Action |
|-----------|--------|
| \|r\| ≥ 0.4 (3DBK or 1EZM, pooled or stratified) | Calibration defensible; fit linear regression, report slope/intercept/RMSE/LOO-RMSE |
| \|r\| < 0.4 for all scoring functions and strata | Null result; do NOT fit calibration; report as null |

---

## Explicit Prediction

**Expected outcome: NULL (|r| < 0.4, pooled and stratified)**

Mechanistic rationale: LasB is a zinc metalloprotease. All known inhibitors (thiol, hydroxamate,
phosphoramidate, carboxylate) bind primarily via direct coordination of the catalytic Zn2+ ion.
This coordination energy is the dominant affinity driver. Vina uses no explicit metal
coordination term — it approximates Zn interactions via a generic Lennard-Jones 4–8 potential
(atom type Zn, +2.000 charge in pdbqt). This approximation is known to underestimate
metal-coordination free energy by 2–4 kcal/mol relative to explicit coordination.

Consequence: across the dataset, the dominant affinity term (Zn-chelation strength) varies
systematically with ZBG class — hydroxamate >> thiol ≈ phosphoramidate >> carboxylate — while
Vina will assign nearly equal weight to all ZBGs (all are treated as generic polar groups).
Within the thiol class, where chelation strength is approximately constant, Vina may in principle
rank compounds by tail contacts; however, the thiol pIC50 range is only 1.49 log units, making
detection of a real signal unlikely at n = 10.

Null prediction is scored CORRECT if |r| < 0.4 for all functions and strata.
Null prediction is scored WRONG if |r| ≥ 0.4 for any function/stratum.

---

## Scoring Details

- Receptor: 3DBK (holo, zinc retained, SO4 and RDF ligand stripped, waters retained)
  and 1EZM (apo, zinc retained)
- Box center (3DBK): 18.721, −5.093, 23.685 (RDF ligand centroid in 3DBK native frame)
- Box center (1EZM): 55.521, 35.882, 20.807 (same center used throughout Phases 4–5)
- Box size: 24 × 24 × 24 Å (identical to Phases 4–5)
- Exhaustiveness: 16, 3 seeds per compound (42, 123, 456), mean best score reported
- Prep pipeline: Meeko mk_prepare_receptor.py, identical to all prior phases

---

## Deliverables

- `phase6_results.md` with all correlations, CIs, scatter plot, and plain verdict
- `results/phase6_compounds.csv` with SMILES, pIC50, source, ZBG, all scores
- Section on zinc: what prep did with it and what this implies for all prior phases

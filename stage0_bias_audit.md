# Stage 0: CASF-2016 Bias Audit

**Delta-learning reimplementation, minimum viable scope.** DeltaVina distribution confirmed
unbuildable (see Phase 7: three diagnosed dependency dead-ends). This reimplements the approach
from scratch, CPU-only, hours not days.

**Scatter plots**: https://claude.ai/code/artifact/99537563-48a6-48a5-b9e4-81027fbaa420

---

## Data access (Stage -1)

PDBbind's own files require registration. Per direction from the user, this used PDB IDs and
pKd/pIC50 labels from a third-party reproduction of the CASF-2016 core set composition
(`test2016.csv` from github.com/guaguabujianle/GIGN — GIGN, IGN, and related repos independently
reproduce this list from the Wang lab's published CASF-2016 papers). Structures were downloaded
directly from RCSB (no registration). This is a **reproduction of CASF-2016, not the official
PDBbind-processed benchmark** — see caveats below.

Verified as a genuine CASF-2016 reproduction: 285 unique PDB IDs, pKd range 2.07–11.82 (matches
published CASF-2016 core-set characteristics: 285 complexes spanning ~10 log units of affinity).

## Extraction pipeline

No PDBbind mol2/pocket files were available, so ligand and protein had to be separated from the
raw deposited structure:

1. **Ligand identification**: among HETATM groups, exclude a blocklist (water, common
   crystallization additives/buffers/ions, glycans, covalently-modified amino acids) and pick the
   largest remaining group by atom count.
2. **Bond order assignment**: fetch the RCSB Chemical Component Dictionary's ideal SMILES for the
   picked residue code, then `RDKit.AllChem.AssignBondOrdersFromTemplate` against the
   crystallographic-coordinate atoms (explicit H's stripped from the raw PDB group first — crystal
   H placement is unreliable and confused template matching).
3. **Protein isolation**: standard ATOM records plus a limited set of metals (Zn, Mg, Mn, Fe, Ca,
   Cu, Na, K — the set meeko's receptor prep has covalent radii for).
4. **PDBQT prep**: protein via `mk_prepare_receptor.py --allow_bad_res --default_altloc A`
   (arbitrary altloc A chosen for every ambiguous residue — a real, stated divergence from
   PDBbind's per-residue curation); ligand via Meeko's `MoleculePreparation.prepare(mol,
   conformer_id=-1)`, which preserves crystal 3D coordinates rather than re-embedding.
5. **Scoring**: native `vina_1.2.7_mac_aarch64 --scoring {vina,vinardo} --score_only`, box = ligand
   bounding box + 8 Å padding (floor 20 Å). No smina needed — the local Vina 1.2.7 binary supports
   both scoring functions natively.

**This is a reproduction, not a replication.** PDBbind applies its own protonation, altloc, and
occupancy curation per complex; this pipeline uses one uniform automated choice throughout.
Absolute scores will not exactly match published ΔvinaRF20/XGB numbers for this reason, independent
of the FreeSASA-for-MSMS substitution that will also apply once SASA features are added in Stage 1.

## Attrition (reported, not hidden)

| Stage | In | Out | Yield |
|---|---|---|---|
| Structure download (RCSB) | 285 | 285 | 100% |
| Ligand extraction | 285 | 266 | 93.3% |
| Receptor + ligand pdbqt prep | 266 | 244 | 91.7% |
| **Final Stage 0 n** | 285 | **244** | **85.6%** |

Ligand-extraction failures (19/285): 3 complexes have no non-blocklisted HETATM group at all
(peptide-chain inhibitors bound as a separate polymer, not a HETATM ligand — out of scope for this
extraction method); 16 failed CCD template bond-order matching (protonation-state mismatches
between the crystal coordinates and the CCD's neutral/ideal representation, or multi-fragment
subcomponent ligands RDKit's template matcher couldn't resolve cleanly).

Receptor-prep failures (22/266): unresolvable altloc/template conflicts in `mk_prepare_receptor`
even with `--allow_bad_res --default_altloc A`, or ligand bounding-box/receptor mismatches.

Full logs: `results/casf2016/extraction_log.csv`, `results/casf2016/stage0_scores.csv`.

## Timing (measured, not estimated)

- Structure download: 285/285 in 7.3s (parallel curl, 12 workers)
- Ligand extraction + bond-order assignment: 266/285 in 0.8s (parallel, 8 workers; CCD lookups cached)
- UniProt mapping (for later family-split use): 285/285 in 3.9s (parallel GraphQL, 12 workers)
- Receptor/ligand pdbqt prep + Vina + Vinardo scoring: 244/266 in 306.4s (1.15s/complex, 8 workers)

**Total Stage 0 wall time: well under the 1-hour box** (dominated by the 5-minute scoring step).

---

## Bias audit (the headline result)

n = 244.

| Metric | Vina | Vinardo |
|---|---|---|
| r(score, HAC) | **−0.592** (p=2.0e-24, CI [−0.663,−0.512]) | **−0.474** (p=4.5e-15, CI [−0.562,−0.380]) |
| Pearson r(score, pKd) | −0.554 (p=5.1e-21) | −0.479 (p=2.0e-15) |
| Spearman ρ(score, pKd) | −0.558 (p=2.2e-21) | −0.469 (p=9.4e-15) |
| **partial r(pKd, score \| HAC)** | **−0.370** | **−0.318** |

(Sign convention: more negative score = tighter predicted binding; higher pKd = tighter measured
binding, so the expected correlation direction is negative — confirmed above.)

### Comparison to this project's LasB finding

| Metric | Vina (LasB, n=11) | Vina (CASF-2016, n=244) | Vinardo (LasB, n=11) | Vinardo (CASF-2016, n=244) |
|---|---|---|---|---|
| r(score, HAC) | −0.915 | −0.592 | −0.858 | −0.474 |
| partial r(pKd, score \| HAC) | −0.088 | **−0.370** | +0.234 | **−0.318** |

**The LasB "Vina is just a size detector" finding does not generalize.** On the diverse CASF-2016
benchmark, roughly 65-70% of Vina's raw pKd correlation survives HAC-correction (partial r = −0.37
vs. raw r = −0.55); on LasB, essentially none of it did (partial r = −0.09 vs. raw r = −0.39). This
localizes the earlier null result to LasB's pocket — most plausibly the zinc-coordination problem
established in Phase 6 (LasB's dominant affinity driver is a Zn2+ interaction Vina cannot model;
most complexes in CASF-2016 do not depend on an unmodeled metal coordination term for their
affinity signal) — rather than a general property of Vina's scoring function. Size dominance
(r(score,HAC) ≈ −0.5 to −0.6) is present on both, so molecular size is always a confound worth
correcting for — but it is not the *whole* story outside LasB the way it appeared to be on LasB.

Vina's raw scoring-power correlation here (Pearson r = −0.554, i.e. |R| ≈ 0.55) is in the range
reported in the literature for Vina-family scoring functions on CASF-2016 (~0.5-0.6), which is a
reasonable sanity check that this reproduction pipeline — despite the registration workaround, the
custom ligand-extraction heuristic, and the FreeSASA-adjacent caveats to come — is producing
scores in the right ballpark, not an artifact of broken structure prep.

---

## Next steps (not yet run, time-boxed per instructions)

- Stage 0.5: profile full Stage 1 feature set (groups A-D, F — no strain) on 50 complexes,
  extrapolate throughput before committing to full-set featurization.
- Stage 1: delta model (XGBoost on pKd − vina_score residual), leave-one-UniProt-family-out CV
  (UniProt mapping already fetched: `results/casf2016/uniprot_map.csv`, 285/285), reproduction gate
  against published ΔvinaRF20 numbers.
- Stage 2: size-decorrelation variants, only after Stage 1's gate is evaluated.

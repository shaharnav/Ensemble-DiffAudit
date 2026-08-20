# Zn²⁺ transplant onto ConforMix/Boltz LasB conformers: the pocket didn't just lose zinc, it collapsed into the vacancy

## Background

Phase 7 established that all 6 `lasb_conformer_beta*` receptors (ConforMix/Boltz-generated,
apo-conditioned) contain zero HETATM records — no Zn²⁺, no Ca²⁺ — while the 1EZM/3DBK crystal
structures retain both. This phase attempted the direct fix: transplant the crystallographic Zn²⁺
and Ca²⁺ coordinates from 1EZM into each conformer (already aligned to the 1EZM frame via a prior
alignment step — `results/lasb_payload/ensemble_receptors_aligned/lasb_conformer_beta*.pdb`, not
re-superposed here) and validate the resulting geometry before accepting any receptor for docking.

## Coordinating residues (identified directly from 1EZM, not assumed)

Geometric search for protein atoms within 3.0 Å of each ion in `1EZM_apo.pdb`:

- **Zn²⁺** (49.202, 37.708, 19.451): His140 NE2 (2.09 Å), His144 NE2 (2.05 Å), Glu164 OE2 (1.84 Å)
  — the HExxH zinc-binding motif of thermolysin-family (M4) metalloproteases, consistent with
  LasB/pseudolysin's known catalytic mechanism.
- **Ca²⁺** (51.525, 45.467, 31.069): Glu172 OE1/OE2, Glu175 OE1/OE2, Asp183 OD1, Asp136 OD2,
  Leu185 backbone O.

## Result: all 6 conformers fail validation, via the same mechanism

| Conformer | His140–Zn (Å) | His144–Zn (Å) | Glu164–Zn (Å) | Ca²⁺ ligands (Å) | Verdict |
|---|---|---|---|---|---|
| β0.0 | 3.20 | **1.26** | 1.87 | 2.51–2.80 | FAIL |
| β0.8 | 3.17 | **1.35** | 1.53 | 2.15–3.42 | FAIL |
| β1.6 | 2.84 | **1.79**† | 2.05 | 2.53–2.84 | FAIL |
| β2.4 | 3.34 | **1.22** | 2.27 | 2.20–3.62 | FAIL |
| β3.2 | 3.10 | **1.43** | 1.61 | 2.24–3.14 | FAIL |
| β4.0 | 2.48 | **1.09** | 1.87 | 2.35–3.33 | FAIL |

†β1.6's closest clash is His144 CE1 (ring carbon, 1.42 Å), not NE2, but the mechanism is the same
side chain.

A real Zn–N(His) coordination bond is ~2.0–2.2 Å. **His144's imidazole nitrogen sits closer to the
crystallographic Zn²⁺ position than a real bond length, in all six conformers independently.** This
is not a borderline geometric call — it's physical atomic overlap. His140 is simultaneously at the
edge of plausible range or beyond it in half the conformers (2.48–3.34 Å). Carboxylate rotamer
symmetry (Glu/Asp O1↔O2 naming ambiguity) was checked and ruled out as an artifact — using
whichever named oxygen is closer changes nothing.

Applying the pre-registered exclusion rule ("exclude any conformer where the coordinating residues
have moved enough that the metal position is chemically implausible") **eliminates all 6
conformers.** Full per-conformer numbers, including every ligand-atom distance and the single
closest non-ligand protein atom to each ion: `results/lasb_ensemble_rmsd/metal_transplant_validation.csv`.

## Why this is a stronger finding than Phase 7's, not a duplicate of it

Phase 7 showed the conformers are missing zinc. This shows *why that's worse than a passive gap*:
with nothing occupying or constraining the site, His144 — one of the three residues whose entire
job is to hold that ion — relaxed into the vacated space. The catalytic site didn't just lose a
ligand, it locally collapsed. Note this also **complicates**, rather than confirms, the earlier
Phase 4 framing ("missing zinc opens volume that ligands can freely fill"): His144 occupying part
of the vacated volume means the net accessible volume change is not simply "more room." The pocket
was reshaped, not just enlarged, and the two effects can't be told apart from this data alone — the
planned zinc-strip control on 1EZM (delete Zn only, hold the rest of the structure fixed) is the
clean way to isolate the pure-volume effect, since it doesn't confound volume loss with an
independent conformational response.

## Generalizable point

Roughly a third of enzymes require a bound metal cofactor for correct active-site geometry.
Structure-prediction models do not emit cofactors by default — they have to be told the metal is
present as part of the input specification. Apo-conditioned generative conformational sampling on
a metalloenzyme, without explicitly specifying the metal, risks generating structures whose active
sites are geometrically invalid in a way that is not visible from sequence, fold, or even backbone
RMSD to the crystal — only from checking the specific coordinating side chains against where the
missing cofactor should be.

## What comes next (not run in this phase)

1. Regenerate the LasB conformer ensemble via ConforMix/Boltz with Zn²⁺ explicitly specified as a
   cofactor in the input (Boltz supports cofactors via CCD code). Whether ConforMix's wrapper
   passes non-protein entities through to the underlying Boltz call, or strips them during guided
   sampling, needs to be checked in the Colab tool itself before assuming this is a drop-in fix.
2. If regeneration with Zn still produces clashing geometry, that is a materially stronger,
   harder-to-dismiss claim than the current one: the method cannot maintain metal-site geometry
   even when explicitly told the cofactor is present.
3. Local side-chain repack of His140/His144 (holding backbone fixed) is available as a labeled
   secondary sensitivity arm — "does the ensemble work if the metal site is manually repaired?" —
   but should not carry the primary RMSD claim, since it repairs the exact geometry the experiment
   is meant to measure.
4. A pure zinc-strip control on 1EZM itself (delete Zn²⁺, change nothing else, redock) isolates the
   volume-only effect that the conformers cannot cleanly isolate, since their pocket collapse
   confounds volume change with an independent conformational response.

## Deliverables

- `lasb_metal_transplant.py` — transplant + validation script (reusable; will re-run cleanly on a
  regenerated conformer set)
- `results/lasb_ensemble_rmsd/metal_transplant_validation.csv` — full per-conformer distance table

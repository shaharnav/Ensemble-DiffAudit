# Pre-registration — Experiment 3: Candida auris DHFR De Novo Design Campaign

This document is committed before any generation or docking for the selected target. It
supersedes the AQP5/DYRK1A screen below only in outcome — the screen itself is retained in
full as the evidence trail for how the target was chosen; nothing in it is deleted or
revised after the fact.

## History — four prior rejections/eliminations, all structural, all pre-run

1. **Xylellain (PDB 3OIS)**: occluding element is a ~50-residue coil propeptide. Two
   controlled ConforMix runs (`structured_regions_only` on and off) placed the propeptide
   24–177 Å from the domain in 5/6 and then 6/6 conformers. Root cause traced to ConforMix's
   own documented behavior (RMSD guidance acts on secondary-structure regions to avoid
   trivial loop-only sampling) and confirmed by direct measurement (biotite SSE: propeptide
   region 16.7% structured). See `xylellain_archive/`.
2. **CLas IMPDH (PDB 6KCF)**: both candidate sites (IMP substrate site, NAD/cofactor site)
   measured on the liganded ortholog 4QM1 at 15.0% and 25.5% helix/sheet respectively — both
   fail the same DSSP floor, stopped before any GPU time per the pre-registered rule. See
   `impdh_archive/`.
3. **AQP5 (human aquaporin-5, hyperhidrosis)** — **eliminated at the holo gate.** All 4
   deposited human AQP5 entries (3D9S, 5C5X, 5DYE, 7STC — confirmed exhaustive via RCSB
   polymer-description search, not title search) contain only a crystallization lipid (PS6)
   or, in one case, a Ni²⁺ ion, as their sole heteroatoms. No small-molecule-bound human AQP5
   crystal structure exists. Separately and independently disqualifying: AQP5's only
   physiologically relevant site is the channel pore itself — a single-file water pathway
   ~3 Å in diameter — which is physically incompatible with DiffSBDD's pocket-conditioned
   generation mode regardless of the holo-gate outcome. Both reasons are recorded because
   either alone would have been sufficient; no substitute (homology model, computational
   holo structure) was attempted, per the pre-registered rule against doing so.
4. **DYRK1A (human DYRK1A kinase, Alzheimer's / Down syndrome)** — **eliminated at the apo
   gate.** All 92 deposited human DYRK1A structures (screened by UniProt accession Q13627,
   not title search, to avoid missing entries with atypical titling) contain at least one
   drug-like small molecule bound at the ATP site — none is buffer/ion/cryoprotectant-only.
   Kinase domains are frequently too conformationally flexible to crystallize unliganded;
   this is a known, general crystallography reality for this fold class, not a search error
   or an artifact of query construction (confirmed by two independent search strategies
   returning the same result). No clean apo structure exists to start ConforMix generation
   from, and manufacturing one by deleting a ligand from a holo structure would not produce
   an independent conformational state — it reproduces the holo backbone exactly (RMSD ≈ 0
   to its own holo by construction), which defeats the purpose of the apo/holo comparison
   entirely rather than approximating it.
5. **C. auris DHFR, first pass (Experiment 2's target screen)**: apo/holo pair originally
   considered there was rejected for a cofactor difference (NADPH bound in holo, absent in
   apo) that would confound any conformational signal with a cofactor-binding signal. See
   "Selected target" below for how this experiment avoids repeating that mistake.

All prior eliminations are structural properties of the target, found by measurement before
generation, not after. The AQP5/DYRK1A screen below applied every lesson from Xylellain and
CLas IMPDH as an explicit, ordered gate set; both candidates were eliminated by that gate set
itself (gates 1 and 4, respectively) rather than downstream in the campaign.

## Candidates and structural hypotheses (as pre-registered, before screening)

- **AQP5** (human aquaporin-5, hyperhidrosis): eccrine sweat gland water channel. Hypothesis
  going in: the only human AQP5 structure (3D9S) is apo, and no drug-like-ligand-bound human
  AQP5 structure is known to exist — this would fail the holo gate before any pocket
  question is even reached. Preferred target if it clears every gate, both for its
  underexplored-condition narrative and because a novel water-channel binding site would be
  a genuinely different pocket topology than anything tested so far. **Confirmed eliminated
  — see History, item 3.**
- **DYRK1A** (human DYRK1A kinase, Alzheimer's / Down syndrome tau hyperphosphorylation):
  classical bilobal kinase fold, ATP pocket at the hinge between the β-sheet N-lobe and
  α-helical C-lobe. Hypothesis going in: kinase ATP pockets are canonically helix/sheet-rich
  and well precedented for structure-based design, so this is expected to clear the DSSP
  gate comfortably; the open question is whether an apo/holo (or DFG-in/DFG-out) pair exists
  with pocket CA RMSD in the 1.0–2.5 Å window — enough conformational signal to be
  informative, not so much that ConforMix's local-perturbation model is being asked to
  reproduce a domain-scale rearrangement. **Confirmed eliminated — see History, item 4.**

## Selected target: Candida auris DHFR (CauDHFR)

Structural series (all four checked before accepting the pair below — HETATM records
verified, not assumed from entry titles):

| PDB | Description |
|---|---|
| 7ZZX | Apoenzyme — no NADPH, no drug |
| 8A0N | Holoenzyme — NADPH bound, no drug |
| 8A0Z | Ternary complex — NADPH + pyrimethamine |
| 8CRH | Ternary complex — NADPH + cycloguanil |

**Apo/holo pair used: 8A0N (apo reference) vs. 8A0Z (holo).** Both contain NADPH — the
conformational difference between them therefore isolates the drug's effect on pocket
geometry, not the cofactor's. This is the specific fix for the failure mode that eliminated
the original human-DHFR attempt in Experiment 2: using 7ZZX (no NADPH at all) as the apo
reference would conflate NADPH-induced and pyrimethamine-induced pocket changes into a single
RMSD, with no way to attribute the observed motion to either cause. 7ZZX and 8CRH are
retained in the repository (verification only) but are not part of the apo/holo comparison
used for the campaign.

## Six screening gates, in order, as applied to AQP5 and DYRK1A during the screen

1. **Holo gate**: a drug-like ligand (MW 150–600; not buffer/ion/cryoprotectant/lipid/
   detergent) co-crystallized at the same site in some deposited structure of the target.
   No holo structure eliminates `gap_recovery_fraction`, the ceiling, and pose validation —
   the metrics that distinguish this experiment from Experiment 1 — so this gate runs first
   and is a hard stop.
2. **Pocket geometry gate**: the site must be a traditional enclosed cavity DiffSBDD can
   condition on — not a channel pore, protein-protein interface, or an allosteric site
   requiring separate validation.
3. **DSSP gate**: pocket-lining residues (any heavy atom within 8 Å of any ligand heavy
   atom, not ligand centroid — the CLas IMPDH screen found centroid-only distance undercounts
   pocket-lining residues by 5x for an elongated ligand and changes the gate's outcome).
   Floor: ≥60% helix/sheet, measured with biotite.
4. **Cofactor check**: apo and holo must be the same construct with no differing bound
   species other than the candidate drug-like ligand (the DHFR failure mode).
5. **Pocket RMSD gate**: apo-to-holo pocket CA RMSD in **1.0–2.5 Å** — above trypsin's noise
   floor (Experiment 1: 0.20–0.27 Å), inside ConforMix's demonstrated local-sampling shell.
6. **Smoke test (Phase 4, only after target selection)**: 2 ConforMix conformers, default
   settings (`structured_regions_only=True`, `--twist-target-stop 2.0` — not reduced, not
   disabled; Xylellain showed disabling it makes things worse). Plausibility gate: pocket
   region CA RMSD to apo ≤ 15 Å, cavity volume not wildly exceeding the holo reference.

## Selection rule (as pre-registered for the screen)

**AQP5 if it passes gates 1–5. DYRK1A otherwise.** If DYRK1A also fails, stop and report
before spending any GPU time or picking a fourth target without discussion.

Both failed. Per that same "stop and report" instinct, CauDHFR was proposed and adopted as a
fourth target only after explicit discussion — not selected unilaterally.

## CauDHFR-specific gate sequence

Same DSSP/RMSD/smoke-test discipline as the screen above, applied directly to the selected
pair rather than as a comparative screen (there is no second candidate at this stage):

1. **HETATM verification** (all four structures) — confirms the table above rather than
   assuming it from entry titles.
2. **Pyrimethamine site check**: confirm the ligand in 8A0Z is drug-like (MW ~249, inside the
   150–600 window) and sits in the folate/substrate site, not the NADPH site — otherwise the
   "drug isolates from cofactor" logic above would not hold.
3. **DSSP gate**: pocket-lining residues = any heavy atom within 8 Å of any pyrimethamine
   heavy atom in 8A0Z (any-atom distance, not centroid — per the CLas IMPDH correction).
   Floor: ≥60% helix/sheet. Hard stop if failed; do not proceed to generation.
4. **Pocket RMSD gate** (only if DSSP passes): 8A0N vs. 8A0Z, pocket CA RMSD and all-atom
   RMSD, plus their ratio. Required range: **1.0–2.5 Å**.
5. **Smoke test** (only if both above pass): 2 ConforMix conformers from 8A0N, default
   settings (`structured_regions_only=True`, `--twist-target-stop 2.0`). Plausibility gate
   applied per the standing rule (same thresholds as Xylellain/CLas IMPDH). Report back
   before the full campaign runs.

All gate results, for AQP5, DYRK1A, and CauDHFR, are recorded in `target_screen_exp3.csv`,
committed before generation starts.

Whichever target is ultimately used **stays in the final writeup regardless of campaign
outcome** — this is a pre-committed decision, not subject to revision after seeing
downstream results.

## Claims boundary — reproduce verbatim in the writeup, near the top, not in an appendix

**Supportable:**
- Relative docking scores across conformational states
- Generation home-field bias, measured
- Whether ensemble conformers move toward the holo state
- Whether ensemble beats a seed-matched noise baseline
- Physicochemical and validity properties of generated compounds
- For DYRK1A: residue-level comparison to other tau kinases (CDK5, GSK-3β)

**Not supportable, and must not be implied:**
- That any generated compound binds the target protein
- That any docked pose is correct
- Any claim of potency, selectivity, or efficacy
- That anything here is a hit, lead, or drug candidate

**Language rule:** write "scored favorably in docking." Never "binds," "inhibits," "hit," or
"lead."

## Stopping rule and guardrails

Gates are checked in order; each is a hard stop, not a suggestion. Do not adjust a threshold
after seeing a result that would fail it. Report every gate result for both candidates in
`target_screen_exp3.csv`, including the one that isn't selected — the screen itself is a
reportable result. Never fabricate results; log every Vina seed and ConforMix setting used
once generation starts; report negative and null outcomes without hedging.

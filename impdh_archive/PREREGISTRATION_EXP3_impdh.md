# Pre-registration — Experiment 3: CLas IMPDH De Novo Design Campaign

This document is committed before any ConforMix or DiffSBDD generation, and before any
docking, for this experiment. Any number reported for this experiment must postdate the
commit hash of this file.

**Target changed from Xylellain to CLas IMPDH.** Xylellain (PDB 3OIS) was abandoned at
Phase 3: two controlled ConforMix runs (`structured_regions_only` on and off) placed the
occluding N-terminal propeptide 24–177 Å from the domain in 5/6 and then 6/6 conformers.
The propeptide is a ~50-residue coil segment; ConforMix's RMSD guidance is documented to act
on secondary-structure regions specifically to avoid trivial loop-only sampling, and the
mask-expansion experiment (disabling `structured_regions_only`) did not fix the problem —
if anything it was worse. That negative result stands on its own (traced to source, not
guessed) and sets the fold-type criterion for this target: the occluding/gating element
must be predominantly helix or sheet, not coil. See `xylellain_archive/` for the full
record.

## Claims boundary — reproduce verbatim in the writeup, near the top, not in an appendix

This is **not** a validation study. There is no holo structure of CLas IMPDH, no known CLas
IMPDH inhibitor, and no experimental measurement of anything produced here.

**Supportable:**
- Relative docking scores across receptor conformational states
- Magnitude of generation home-field bias, in kcal/mol
- Whether ensemble conformers alter pocket cavity volume
- Whether ensemble best-of-N exceeds a seed-matched rigid baseline
- Physicochemical and validity properties of generated compounds
- Sequence/structural differences between the CLas and human IMPDH pockets

**Not supportable, and must not be implied:**
- That any generated compound binds CLas IMPDH
- That any docked pose is correct
- That ensemble conformers correspond to states the protein populates
- Any claim of potency, selectivity, or efficacy
- That anything here is a hit, lead, or drug candidate

**Language rule:** write "scored favorably in docking." Never "binds," "inhibits," "hit," or
"lead."

## Target and structure

- **PDB 6KCF** — CLas IMPDHΔ98-201 (CBS/Bateman-domain deletion construct), apo, 2.55 Å,
  homotetramer (chains A–D, 390 residues each per SEQRES). Verified locally
  (`verify_6kcf.py`, `results/experiment3/verification_6kcf.json`):
  - **Apo confirmed**: only water HETATMs (42), no drug-like ligand in any chain.
  - **Catalytic Cys located by sequence motif, not by the paper's residue number.** The
    paper (Nan et al., *Molecules* 2020, PMID 32423054) reports "Cys309," which is
    UniProt/full-length numbering (verified: `GSIC` — the conserved single-occurrence IMPDH
    nucleophile motif in UniProt C6XG59 — places the Cys at position 309, matching the
    paper). The deposited coordinate file uses its own numbering, confirmed to differ: the
    same motif places the catalytic Cys at **residue 303** in 6KCF's own numbering,
    consistent across all four chains. Every downstream script uses **Cys303 (deposited
    numbering)**, not 309.
  - The catalytic Cys sits immediately after a disordered "flap loop" (residues 297–302,
    unmodeled in all four chains) — consistent with the paper's own note of missing density
    in "the flap loop and a C-terminal loop." This flap loop forms part of the NAD-binding
    subsite (see Phase 2 pocket selection) and its disorder here is itself informative for
    the DSSP gate.
  - CBS-domain deletion boundary confirmed structurally: modeled sequence jumps from residue
    88 to residue 206 (chain A/B/D) or 209 (chain C) — a ~104–117-residue gap encompassing
    both the engineered deletion (98–201 in WT numbering, replaced by a single residue per
    the paper) and additional flanking loop disorder.
  - Oligomeric state: homotetramer. Catalytic Cys303 SG is 8.09 Å from the nearest atom of
    the closest neighboring chain (B) — too far to be pocket-lining at this project's
    standard 8 Å cutoff, and consistent with IMPDH's catalytic site being intra-subunit
    (the CBS/Bateman domains mediate tetramer contacts and allosteric regulation, not
    catalysis, in every characterized IMPDH structure). **R_apo is chain A alone**,
    consistent with the single-chain-receptor convention used throughout this project.

## Construct limitation

**The CBS-domain deletion (Δ98–201) is a construct limitation, not a modeling choice made
in this experiment.** It was necessary for crystallization in the original study. The
full-length enzyme's CBS/Bateman domains are known allosteric regulatory elements in IMPDH
family members; their absence here means nothing about allosteric conformational states can
be inferred, and this limitation must be restated wherever conclusions are drawn.

## Receptor and ligand plan

- Up to **8 receptors**: `R_apo` + up to 6 ConforMix conformers (contingent on the Phase 3
  smoke test and full-ensemble plausibility gate) + human IMPDH2 is used only for the
  selectivity analysis (Phase 6b), not as a docking receptor.
- **45 generated ligands**: 15 each conditioned on `R_apo`, `C_min`, `C_max`.
- **7 reference compounds** (`reference_compounds.csv`, RDKit-valid, verified against
  PubChem/primary literature): IMP (substrate, PubChem CID 8582), XMP (product, CID 73323),
  mycophenolic acid (uncompetitive NAD-site inhibitor, CID 446541), ribavirin
  monophosphate (active ribavirin metabolite), merimepodib/VX-497 (NAD-site inhibitor),
  a Cryptosporidium/bacterial-type-IMPDH-selective urea inhibitor ("7b," EC50 6 nM, 670-fold
  selective; Macpherson et al., PMC3635066), and tiazofurin (NAD-analog prodrug inhibitor,
  CID 457954). These are a **scale check**, not ground truth for CLas IMPDH — used to sanity
  the score range in Phase 5.

## Pocket choice (Phase 2a) — GATE FAILED, experiment stopped here

Both candidate sites were measured (`impdh_site_geometry.py`,
`results/experiment3/phase2a_dssp_gate.json`) directly on the liganded ortholog PDB 4QM1
(*B. anthracis* IMPDH, same CBS-deletion construct type, co-crystallized with both IMP and
NAD-subsite inhibitor 39H/D67) rather than on apo CLas 6KCF, which is missing density at
both sites and cannot support a fold-composition measurement on its own backbone.

Pocket-lining residues were defined as any residue with a heavy atom within 8 Å of *any*
ligand heavy atom (not the ligand centroid — a centroid-only cutoff undercounted the
elongated NAD-site ligand's pocket fivefold on the first pass, n=3 vs. n=55, and changed the
gate's outcome; caught and corrected before trusting it). The NAD/cofactor site was
confirmed to be a genuine subunit interface (39H found 3.22 Å from chain B) and scored across
both chains.

| Site | n pocket-lining residues | Helix/sheet fraction | Gate (≥60%) |
|---|---|---|---|
| IMP substrate site | 60 (58 chain A + 2 chain C) | 15.0% | **FAIL** |
| NAD/cofactor site | 55 (20 chain A + 35 chain B) | 25.5% | **FAIL** |

**Neither site passed. Per the pre-registered stopping rule, the experiment stops here —
no ConforMix run, no DiffSBDD generation, no docking.** CLas IMPDH has the same defect that
ended the Xylellain attempt: both plausible binding sites are loop-dominated, not
predominantly helix/sheet, so ConforMix's RMSD guidance (which acts on secondary-structure
regions) cannot be expected to sample either site plausibly. This is a second, independent
instance of the same applicability boundary, not a new failure mode.

## Metric definitions (Phase 7)

```
home_advantage        = mean(score in own conditioning receptor)
                        − mean(score in the other receptors)     # per ligand, paired

noise_corrected_delta = apo_seeded_max − ensemble_best            # per ligand
                                                                  # positive = ensemble better
```

Affinities are negative (Vina convention); **positive** delta means the ensemble scored more
negative (better). No gap-recovery / opening-state metric exists for this target — there is
no holo structure, and the Xylellain R2-style propeptide-deletion proxy has no analogue here.

## Stopping rule

Analysis runs once, after all docking completes. Gates (DSSP floor, conformer plausibility,
smoke test) are checked in execution order and are hard stops — if any fails, the phase after
it does not run, and the failure is reported as-is.

## Guardrails

Never fabricate results. Never tune parameters (Vina exhaustiveness, ConforMix twist
strength, chemistry filter thresholds) to enlarge an apparent effect after seeing results.
Log every Vina seed and exhaustiveness value used. Report negative and null outcomes without
hedging or searching for a favorable subset. If a gate fails, stop and report — do not adjust
the threshold unilaterally.

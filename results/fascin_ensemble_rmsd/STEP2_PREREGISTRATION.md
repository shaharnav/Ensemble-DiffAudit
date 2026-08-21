# Fascin Step 2: pre-registration, committed before generation or validation

Per protocol: the LasB experiment ended with n=2 validated conformers because
validation happened after generation and attrition wasn't budgeted for. This
commits the target ensemble size and the fallback plan *before* any conformer
is generated or looked at, so neither can be adjusted after seeing results.

## Pre-registered target

**Generate 16 conformers** via ConforMix/Boltz from the fascin apo structure
(3LLP, chain A), spanning a beta/twist sweep (beta = 0.0, 0.4, 0.8, ..., 6.0 --
16 evenly spaced values, matching the density of the LasB sweep scaled to a
16-point budget instead of 6).

**Target ensemble size after attrition: >= 8 validated conformers.**

## Validation gates (applied immediately after generation, not at docking-prep time)

A conformer survives if it passes all of:
1. Meeko receptor prep succeeds (diagnosed by direct atom-distance measurement
   if it fails, not assumed).
2. No non-adjacent-residue backbone clash (direct distance check, the same
   diagnostic that caught beta2.4/beta4.0's defects in the LasB run).
3. PoseBusters / structure validity checks pass.
4. Chain integrity intact (fascin is a monomer in the biological sense --
   this reduces to: single continuous chain, no gaps introduced by the
   generation process in the region used for docking).
5. Pocket-region backbone stays physically reasonable (no requirement to
   preserve the *apo* pocket geometry specifically -- unlike a
   metal-coordination check, there's no cofactor here to validate against).

No metal cofactor exists for fascin, so the LasB-style Zn/Ca transplant step
does not apply here -- confirmed by reading the ConforMix source path (see
below) and checked directly against 3LLP: zero non-polymer entities of any
kind in the apo structure.

## If fewer than 8 survive

Report the actual number and the specific, measured reason for every
exclusion. Proceed with whatever survives, labeled as small-n / preliminary
throughout downstream steps -- per protocol, do NOT generate replacement
conformers to backfill the target, since that would be post-hoc selection.

## Generation environment

Boltz run via the colab-mcp tool (per instruction: no Kaggle for this phase),
using `run_conformixrmsd_boltz.py`'s existing twist-guided sampling CLI,
already used successfully for the LasB regeneration. No cofactor entity is
added to the fasta input (none needed). No structural water is retained
either (none identified in the fascin literature reviewed during Step 1d).

Committed before opening the Colab connection or writing any generation code.

# Pocket-change gate: definition, pre-registered before screening any new target

Committed before running the CryptoSite/PocketMiner or candidate-list (A-F) screen below,
per instruction: settle this once, so the HIV-1 mistake isn't repeated across ~30 targets.

## What went wrong in the HIV-1 screen

The gate uses pocket CA RMSD (apo->holo) and DSSP structured fraction of that
displacement. The pocket residue set was defined as: apo residues with an atom
within 5A of the ligand, **after aligning the ligand into the apo frame**.

Direct measurement shows why this missed the flap: HIV-1 protease's flap-tip
residue (Ile50) sits 2.85-3.87A from the ligand in every one of the 9 holo
structures checked (genuinely in contact once closed) -- but its **apo**
coordinates sit 5.8-6.7A from the aligned ligand position, because in the open/
semi-open apo state the flap has already swung away from where it ends up.

The bug was not the radius. It was the reference frame: checking whether the
*apo* copy of a residue is near the ligand systematically excludes residues
that move the most, because "moved away in apo" and "far from the apo-frame
ligand" are the same thing. A gate built this way is structurally biased
against detecting the induced fit it's supposed to measure.

## Pre-registered definition

1. **Contact-shell residues**: any residue with >=1 heavy atom within 5A of
   >=1 ligand heavy atom, **computed in the native holo (bound) structure**,
   not from the apo-aligned copy. This is standard binding-site-definition
   practice (LigPlot/PocketMiner/CryptoSite all do this) and is frame-correct:
   it asks "what does this residue touch when bound," not "does the apo copy
   happen to still be nearby."

2. **Access-shell residues**: contact-shell residues UNION any residue within
   +/-2 sequence positions (same chain) of a contact-shell residue. Cheap,
   automatable way to capture the rest of a lid/flap/loop element that
   partially, but not fully, contacts the ligand -- without per-target manual
   annotation of "gating" residues, which doesn't scale across ~30 candidates.

3. Both residue sets are then mapped by residue number onto the **apo**
   structure, and apo->holo CA displacement (after global CA superposition) is
   computed and reported for each. Report both; **access-shell pocket CA RMSD
   is the primary gate metric** (contact-shell is reported alongside for
   comparability/diagnosis, since it's what the LasB baseline and the first
   HIV-1 pass used).

4. **DSSP on the full biological assembly**, generated from the deposited
   BIOMT/assembly records (gemmi `make_assembly`), never a lone chain/subunit.
   Verified this matters: monomer-only DSSP on 1HHP gives 54.5% H/E content;
   dimer-aware DSSP gives 61.6%, and the flap-tip residue itself only reads as
   strand (E) when the dimer partner is present.

5. **Gate thresholds** (unchanged in spirit, restated precisely):
   - Primary: access-shell pocket CA RMSD must **substantially exceed** LasB's
     ~1.5A baseline. Interpreting "substantially exceed" as **>= 2.0A** (real
     margin, not a coin-flip over measurement noise).
   - Secondary: **>= 60%** of that access-shell displacement must come from
     DSSP-structured (H/E) residues.
   - Both must pass. Reported per apo/holo pair, not pre-averaged, so
     per-pair variance is visible.

## What this does NOT change

- The occupancy/B-factor/resolution/mutation inclusion criteria from Step 1b.
- The requirement to superpose on CA atoms and report alignment RMSD.
- Symmetry-corrected RMSD, matched compute, and every other downstream rule
  from the original 7-step protocol -- this document only fixes the Step 1c
  gate mechanics.

## Committed before screening

This file is committed before the CryptoSite/PocketMiner screen and the A-F
candidate list are run, so the definition can't be adjusted retroactively
based on which targets it happens to pass or fail.

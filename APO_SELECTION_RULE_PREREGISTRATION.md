# Apo-deposition selection rule: pre-registered before re-testing TEM-1, HSP90N, BACE1

Committed before fetching or looking at any alternate apo depositions for these
three targets, per instruction: iterating apo structures until a target passes
the gate is structure-shopping (same error class as loosening a threshold).
The legitimate version is a rule chosen for reasons unrelated to the outcome,
pre-registered, and applied uniformly -- including to targets that already
passed with their first-tried structure.

## The rule

For a given target (UniProt accession), among all PDB depositions for that
accession:

1. **Apo set** = entries with no non-polymer entity outside the standard
   crystallization-additive blocklist (`BLOCKLIST_LIGANDS` in
   `screen_ligand_availability.py` -- glycerol, PEG fragments, buffer ions,
   sulfate/phosphate, etc.). This is the same blocklist already used for the
   ligand-availability pass, applied here to classify apo vs. holo instead of
   counting drug-like ligands.
2. **Selection within the apo set**: highest crystallographic resolution,
   full stop. Tie-break (if ever needed) by most complete residue coverage
   (fewest missing/unmodeled residues in the region matched during screening).
3. **Isoform/sequence**: if a target has multiple human isoforms deposited
   separately (e.g. HSP90-alpha vs HSP90-beta), fix the isoform choice to
   whichever has more total PDB depositions for that UniProt accession (a
   quality/availability criterion, not a motion criterion), and apply it to
   both the apo and holo pick.
4. **Holo pick**: among entries with a real (non-blocklisted) ligand, again
   highest resolution, restricted to entries already meeting Step 1b's
   existing inclusion bar (resolution <= 2.5A, sane occupancy/B-factor).

Nothing in this rule references pocket RMSD, flap displacement, or DSSP
structured fraction -- the three quantities the gate is built from. It cannot
be satisfied by trying candidates until one passes, because it does not look
at the outcome at all.

## Why resolution instead of "most open" or "closest to holo"

Both of those would be direct selection on the outcome variable and would
invalidate the gate the moment they were used. Resolution is a data-quality
criterion that is, on priors, uncorrelated with how much a flap/lid happens to
have moved in that particular crystal form. It is also the same criterion
already used to break ties in the original Step 1a apo selection (1HHP vs
2PC0 vs 3PHV), so this isn't a new standard introduced just for this pass.

## What this does NOT do

It doesn't guarantee any of the three targets pass. If the highest-resolution
apo deposition for a target still shows little motion, that's the answer --
not a reason to look for a fourth one. The three results below are reported
whichever way they land, with continuous values (not binary pass/fail at a
hard cutoff), per the separate instruction that a result 0.03A from a
threshold is a rounding error, not a finding.

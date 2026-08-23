# Fascin Step 3: final result — ConforMix/Boltz twist-guidance does not produce viable conformers for this target, at any tested setting

## STATUS: FAILED — experiment terminated at Step 3

Fascin is closed out as a **failed target** for this experiment. No usable
conformer ensemble was produced, so Steps 4-7 (docking, RMSD measurement,
secondary metrics, final analysis) were never started for fascin — there is
nothing downstream of this document. This is a clean negative result, not an
inconclusive or paused one: the pre-registered retry ladder was exhausted per
its own fixed acceptance criterion, and per that pre-registration no further
generation settings will be tried.

## Bottom line

**0 of 16** original conformers, and **0 of 8** ladder-retry conformers
(2 settings x 4 targets each), pass a direct all-heavy-atom clash gate
(<2.0A between heavy atoms of residues >=4 apart in sequence, excluding the
crystallographically-unresolved N-terminus). Per the pre-registered retry
plan (`STEP3_RETRY_PREREGISTRATION.md`), the ladder was capped at two new
settings with a fixed acceptance criterion (>=3/4 smoke-test conformers
clash-free) decided before generating anything. Neither setting cleared it.
No further settings will be tried. This is the final result for fascin under
ConforMix/Boltz twist-guidance.

## What was tried, and the numbers

| Setting | subset-residues | twist_strength | n | mean clashes/conformer | pass rate |
|---|---|---|---|---|---|
| Original (Step 3) | 139/493 (28.2%) | 15.0 | 16 | 29.7 | 0/16 |
| Ladder pt.1 | 88/493 (17.8%) | 7.5 | 4 | 13.0 | 0/4 |
| Ladder pt.2 (final) | 88/493 (17.8%) | 3.75 | 4 | 8.2 | 0/4 |

The dose-response is monotonic and real — halving twist_strength roughly
halves the mean clash count each step — but extrapolating the trend, it does
not reach zero in any regime that would still constitute meaningful
RMSD-guided perturbation. At `twist_strength=3.75` the model is barely being
guided at all, and still produces 7-10 clashes per conformer.

## Diagnostic: why does this happen (not just that it happens)

Two hypotheses tested directly before any retry:

**Repeat confusion (rejected).** Fascin has 4 tandem CATH beta-trefoil
domains (2.80.10.50; D1=7-138, D2=139-261, D3=262-386, D4=387-493). If
guidance couldn't distinguish which repeat it was pulling toward the target,
clashes should concentrate at structurally-equivalent positions across
domains. Measured: 8/69 (12%) of inter-domain clashes fall at equivalent
positions, vs 15% under a random-pair null. No enrichment.

**Local packing collapse (supported).** 355/475 (75%) of clashes in the
original 16-conformer run are intra-domain. The guidance is breaking a single
domain's own internal packing, not causing domains to overlap each other.
This is consistent with the observed dose-response (softer guidance = better,
monotonically) and explains why shrinking the subset window alone (28.2% ->
17.8%, matched to LasB's validated fraction) wasn't sufficient on its own —
the remaining guided residues were still being pulled hard enough to break
local geometry within their own domain.

**A specific recurring defect**: the pair Arg197-Arg201 clashes in all 4
conformers at the softest tested setting (twist_strength=3.75), independent
of beta target. This particular turn appears to be a structural bottleneck
that any perturbation strong enough to move the pocket region also disrupts
locally — a candidate for the "hardest to guide without breaking" site in
this fold, though this wasn't investigated further given the ladder is capped.

**Methodologically important, independent of the retry outcome**: 13 of the
original 16 conformers passed the global CA RMSD gate (<=2.5A to the apo
crystal) while every one of the 16 had severe local clashes. A global RMSD
gate — the kind most pipelines would use as their only structural sanity
check — is completely blind to this failure mode. Local, all-heavy-atom
clash checking is necessary; global backbone RMSD is not sufficient.

## What this does and doesn't establish

- It does **not** show that ConforMix/Boltz twist-guidance is broken in
  general — LasB (298 residues, single domain) produced a genuine, validated
  4/6 pass rate on its cofactor-corrected regeneration.
- It **does** suggest a real scope limit: large (493-residue), multi-domain,
  internally-repeated folds may not be well served by this RMSD-guided
  approach at settings carried over from a much smaller single-domain
  protein, and softening those settings within a reasonable range does not
  fix it — only slows the rate of failure.
- The repeat-confusion mechanism specifically proposed as a candidate
  explanation was tested and rejected; the actual mechanism (local packing
  collapse under guidance, worse in some regions than others) is different
  from what was hypothesized going in, and is itself informative for anyone
  considering this method on a large or multi-domain target.

## Deliverables

- `results/fascin_ensemble_rmsd/step3_validation_geometry.csv` — original 16-conformer gate results
- `results/fascin_ensemble_rmsd/step3_clash_domain_analysis.csv` — all 475 clashes mapped to domains
- `results/fascin_ensemble_rmsd/ladder_s75_results.txt`, `ladder_s375_results.txt` — ladder retry results
- `STEP3_RETRY_PREREGISTRATION.md` — the pre-registered retry plan and acceptance criterion
- Raw conformer CIFs: `conformers_raw/` (original 16), `ladder_s75/`, `ladder_s375/` (4 each)

# New target pre-screen — results

Screened the full 49-pair PocketMiner/CryptoSite apo/holo pool
(`results/target_screen/pocketminer_pairs.csv`) against the six pre-registered gates
(`GATE_DEFINITION_PREREGISTRATION.md`, `PREREGISTRATION_EXP3.md`).

## Bug fix, applied before screening

`cryptic_pocket_screen.py`'s ligand-code parser (`lstrip("#0123456789x")`) was stripping
leading digits from *real* ligand codes (`4CS`→`CS`, `181`→``, `03F`→`F`, `89A`→`A`), not
just the intended `NxCODE` multiplicity annotations (`2xCHD`→`CHD`). This caused 10 of 49
pairs to fail with `ligand_not_found`. Replaced with a parser that only strips a genuine
`^\d+x` prefix and handles `:NNN`/`-ION` suffixes. One source-data row (4V38/4V3B) had a
literally wrong ligand code (`C`); corrected to `C5P` (confirmed from the deposited HETATM
records). Re-running recovered all 10 pairs to real results — none of them pass the
structural gates (RMSD/DSSP), so this doesn't change the shortlist, but the pool is now
correctly characterized (43 of 49 usable; the other 6 have no deposited holo structure and
are not recoverable).

## Gate 3+5 (DSSP / pocket RMSD, already correctly implemented): 9 of 43 pass

## Gates 1/2/4/6 applied to those 9

| Pair | Protein | Ligand | Verdict |
|---|---|---|---|
| 1S2O/1U2S | Sucrose-phosphatase | GLC (glucose) | REJECTED — enzyme's own substrate |
| 1Y1A/1Y1A | CIB1 | GSH (glutathione) | REJECTED — endogenous cofactor |
| 2CEY/6H76 | SiaP | SLB (sialic acid) | REJECTED — transport protein's own cargo |
| 3FVJ/2B03 | Phospholipase A2 | TUD (bile acid) | REJECTED — endogenous lipid |
| 3P53/6I11 | **Fascin** | H0H | Would pass — **excluded, already attempted and terminated (commit 2a0d38e)** |
| 2HQ8/2HPS | Obelin | CTZ (coelenterazine) | BORDERLINE — natural luciferin substrate; not literally blocklisted but not a drug candidate |
| 4R72/4R74 | AfuA | F6P (sugar-phosphate) | REJECTED — transport protein's own cargo |
| 1OK8/1OKE | (detergent soak) | BOG (octyl glucoside) | REJECTED — crystallization detergent; its 8.58Å RMSD is likely a soak artifact |
| **2WGB/2V57** | **LfrR** (TetR-family regulator, *M. smegmatis*) | **PRL (proflavine)** | **VIABLE — passes all six gates** |

## Bottom line

**Only one confirmed-viable new target: LfrR + proflavine (2WGB apo / 2V57 holo).**
- Holo gate: proflavine, MW 209, an actual approved antiseptic — genuinely drug-like.
- Cofactor gate: apo has only HOH/SO4 (blocklisted), holo adds only PRL + blocklisted
  additives — clean.
- Apo-selection gate: trivially satisfied — only 2 PDB depositions exist for this UniProt
  (Q58L87) and 2WGB is the sole apo one.
- Pocket geometry: enclosed hydrophobic regulatory-domain pocket, precedented in the
  TetR-family SBDD literature (EthR, QacR) — flagged, not silently assumed, because it is
  allosteric (ligand pocket ≠ DNA-binding site) and `GATE_DEFINITION_PREREGISTRATION.md`'s
  language explicitly excludes "an allosteric site requiring separate validation."

**One borderline case**: Obelin + coelenterazine passes the *literal* text of the holo gate
(MW in range, not lipid/buffer/cryoprotectant/detergent) but is a natural bioluminescent
substrate, not a drug-candidate-style ligand — same category problem as the four rejected
natural-substrate/cofactor pairs above, which the pre-registered gate text doesn't actually
have a clause for. This is a real gap in the gate definition, not a judgment call I should
make unilaterally.

This is nowhere near 10. The 43-pair PocketMiner/CryptoSite pool is close to exhausted:
only 9 pairs clear the structural gates at all, and 7 of those are natural
substrate/cofactor/detergent complexes that were never going to be viable de novo design
targets regardless of RMSD.

## Widened search (round 2)

Since LfrR (a TetR-family regulator) was the one hit, searched four more well-precedented
TetR/multidrug-regulator-family targets with literature apo/holo pairs and drug-like
effectors:

| Target | Apo/holo | Ligand | Result |
|---|---|---|---|
| TtgR (*P. putida*) | 2XDN/2UXH | Quercetin | REJECTED — pocket RMSD 0.40Å, far under the 2.0Å floor |
| TtgR | 2XDN/2UXP | Chloramphenicol (real antibiotic) | REJECTED — pocket RMSD 0.35Å |
| EthR (*M. tuberculosis*, major TB SBDD target) | 1U9N/3G1M | Synthetic inhibitor BDM31381 | REJECTED — pocket RMSD 0.30Å |
| MtrR (*N. gonorrhoeae*) | 6OF0/8FW3 | Testosterone | REJECTED — pocket RMSD 0.36Å |
| QacR (*S. aureus*) | — | — | REJECTED — no unliganded deposition exists among 14 PDB entries (same failure mode as DYRK1A) |
| SimR (*Streptomyces*) | 2Y2Z/2Y30 | Simocyclinone D8 | REJECTED — ligand MW 932, over the 600 ceiling |

All four RMSD failures are genuine, not gate-tuning artifacts: these regulators' ligand
pockets are essentially rigid (apo already sits in a ligand-competent conformation), so
there's no meaningful conformational change for ConforMix's local-perturbation model to
target — the same floor problem that eliminated LpqN (0.56Å) in Experiment 2. LfrR's 4.58Å
shift is evidently not representative of the TetR family as a whole.

## Widened search (round 3 — large push)

Went further into classic literature cryptic-pocket/induced-fit systems, spanning kinases,
phosphatases, a cytokine PPI hot spot, and a viral polymerase's allosteric sites:

| Target | Apo/holo | Ligand | Result |
|---|---|---|---|
| Interleukin-2 | 1M47/1PY2 | Hot-spot inhibitor FRH | REJECTED — ligand MW 663 (over 600 ceiling) *and* RMSD 0.60Å; also a PPI-interface pocket, not an enclosed cavity |
| PTP1B | 1SUG/1T49 | Allosteric BB-site inhibitor 892 | REJECTED — ligand MW 658 (over ceiling) *and* RMSD 1.67Å |
| HIV-1 reverse transcriptase (NNRTI pocket) | 1DLO/1FK9 | Efavirenz (real approved drug, MW 316) | REJECTED — RMSD 3.53Å **passes primary**, but only 47% DSSP-structured (needs 60%) |
| p38-alpha MAPK (DFG-out pocket) | 1WFC/1KV2 | BIRB-796, MW 528 | REJECTED — RMSD 1.42Å, under the 2.0Å floor |
| HCV NS5B (thumb site 2) | 1C2P/2GIQ | NN2, MW 427 | REJECTED — pocket essentially rigid, RMSD 0.20Å |
| HCV NS5B (palm site) | 1QUV/1YVF | PH7, MW 438 | REJECTED — RMSD 0.52Å |

Efavirenz/HIV-RT came closest — it's the only round-3 candidate to clear the primary RMSD
gate, and misses the secondary DSSP-structured-fraction gate by a real margin (47% vs. 60%
required), not a rounding error.

## Bottom line after three rounds

**Still 1 confirmed-viable new target: LfrR/proflavine.** 21 individual candidate pairs
across roughly 15 distinct protein families — bacterial multidrug regulators, a cytokine PPI
hot spot, a phosphatase allosteric site, a kinase DFG-out pocket, a viral polymerase's two
allosteric sites, and a viral reverse-transcriptase NNRTI pocket — produced exactly one that
clears every gate. This is not a search-execution problem; the gate combination
(≥2.0Å pocket CA RMSD *and* ≥60% DSSP-structured, on top of the holo/cofactor/apo-selection
gates) is a genuinely stringent filter, and it matches this codebase's own established
pattern: most well-chosen candidates fail pre-screen (Xylellain, CLas IMPDH, AQP5, DYRK1A,
LpqN, QacR, SimR, and now 20 more, all failed here). Reaching 10 viable targets this way
would plausibly require screening on the order of 100+ more candidate pairs by hand at this
observed ~5% hit rate — a large, slow, manual-curation-bound effort. Recommend discussing
with the user whether to keep pushing at this rate, systematically pull PocketMiner's full
published dataset (which is larger than the 49-pair subset already vendored here) instead of
hand-picking one family at a time, or accept LfrR as the sole confirmed target and move to
staging it for Phase 4.

# Fascin: Step 1 summary (target and structure selection)

## 1a. Target and apo structure

**Fascin-1** (UniProt Q16658), an actin-bundling protein with a druggable
small-molecule groove between beta-trefoil domains 1/2 (the site targeted by
the NP-G2-044-class oncology inhibitor series). Selected after LasB (null,
gate not built to detect the failure mode) and HIV-1 protease (failed its own
gate at n=9 holo structures) via a systematic screen: 33 curated PocketMiner
apo/holo pairs plus 5 named candidates (BACE1, HSP90N, TEM-1, PPARgamma,
adenylate kinase), scored jointly on pocket-opening magnitude and drug
availability (full screen: `cryptic_pocket_screen_results.csv`,
`cryptic_pocket_screen_with_availability.csv`).

**Apo structure: 3LLP, chain A** (1.8A resolution, wild-type sequence
confirmed by direct diff against canonical UniProt Q16658, zero non-polymer
entities of any kind -- not even crystallization additives). Selected over
the original screening pair's apo (3P53, 2.0A) per the pre-registered,
outcome-blind apo-selection rule (`APO_SELECTION_RULE_PREREGISTRATION.md`):
highest-resolution genuine apo deposition for the UniProt accession. Verified
directly that switching apo structures did not weaken the gate result -- it
strengthened it (see 1c).

Biological assembly: 3LLP's asymmetric unit contains two copies (chains A/B);
fascin functions as a monomer (each chain independently binds two actin
filaments), so this is crystallographic packing, not an obligate oligomer.
Chain A used throughout.

## 1b. Holo ligand set

All PDB entries linked to UniProt Q16658 were enumerated (28 total). Applying
resolution <=2.5A and requiring a genuine (non-crystallization-additive)
ligand: **11 structures survive**, all part of one structure-based drug-design
series (compounds 1-24 plus BDP-13176, from the same published campaign) --
below the 15-25 target range because that's what exists for this less
extensively studied target, reported plainly rather than padded.

| PDB | Ligand | Resolution | MW | Occupancy | Mean B-factor |
|---|---|---|---|---|---|
| 6I0Z | GZQ | 1.77 | 218.1 | 1.00 | 36.5 |
| 6I10 | GZK | 2.10 | 413.3 | 1.00 | 46.0 |
| 6I11 | H0H | 1.67 | 358.4 | 1.00 | 29.6 |
| 6I12 | H08 | 1.65 | 392.8 | 1.00 | 26.8 |
| 6I13 | H0Q | 1.79 | 392.8 | 1.00 | 29.1 |
| 6I14 | GZN | 1.73 | 427.3 | 1.00 | 27.9 |
| 6I15 | GZT | 1.91 | 322.3 | 1.00 | 36.0 |
| 6I16 | H0B | 2.00 | 428.3 | 1.00 | 33.2 |
| 6I17 | GZW | 1.56 | 511.4 | 1.00 | 20.2 |
| 6I18 | H0N | 1.49 | 497.4 | 1.00 | 22.7 |
| 9GS6 | A1IPD | 1.97 | 380.2 | 1.00 | 38.8 |

Full occupancy across all 11 -- no partial/alt-conf ambiguity requiring
exclusion (unlike LasB's V85/IEV). MW range 218-511 Da (GZQ, a small fragment
hit, is the low end; still inside the 150-700 protocol bound).

**Site-consistency check** (fascin's analogue of "distance to catalytic Zn,"
since fascin has no metal cofactor): each ligand's centroid, transformed into
the apo frame, falls 1.7-4.7A from the pocket residue set independently
established during target screening (3LLP->6I11 access-shell residues:
14, 16, 48, 60, 93-95, 101, 103, 134, 214-217, 224). All 11 confirmed binding
at the same site -- none excluded.

**Chemical diversity** (Morgan r=2 Tanimoto, pairwise): mean 0.40, median
0.37, max 0.863 (H0H-H08, same scaffold, expected SAR-series neighbors), no
pair exceeds the 0.9 near-duplicate threshold -- nothing pruned. Diversity is
real but moderate: this is one campaign's structure-activity series, not
independently-discovered chemotypes the way HIV-1 protease's 9 approved drugs
were. Reported as a limitation, not concealed.

## 1c. Conformational change (gate re-confirmation on final apo choice)

Measured directly on the final 3LLP(A)->6I11 pair, using the pre-registered
gate definition (holo-frame contact shell, +/-2 residue access-shell
expansion, DSSP on the deposited biological assembly):

| Metric | Value | Gate | LasB (for comparison) |
|---|---|---|---|
| Pocket CA RMSD (access-shell) | **5.95 A** | >= 2.0A -> PASS | ~1.5A (apo/2PC0 best case) |
| Structured fraction | **71.7%** | >= 60% -> PASS | 27-50% |
| Global CA RMSD | 3.87 A | -- | ~1.0-1.3A |

Both gates pass with real margin, using the higher-resolution apo (3LLP) --
stronger than the original screening pair (3P53->6I11: 5.25A / 68.5%). This
was verified, not assumed, before locking in 3LLP as the apo structure.

## 1d. Provenance

| Structure | Role | Resolution | Mutations | HETATM (non-water) |
|---|---|---|---|---|
| 3LLP (chain A) | apo | 1.8A | none (wild-type, confirmed by sequence diff) | none |
| 6I0Z, 6I10-6I18, 9GS6 | holo (Condition C sources, Step 3/4) | 1.49-2.10A | none listed | one drug-like ligand each (see table above) |

No structural water retained in any receptor at this stage (no conserved
bridging water analogous to LasB's zinc or the HIV-1 flap water has been
identified for fascin in the literature reviewed) -- flagged here for
awareness, to be revisited if Step 3 docking shows systematic pocket-shape
discrepancies.

## Deliverables

- `results/target_screen/cryptic_pocket_screen_results.csv` -- 33-pair PocketMiner gate screen
- `results/target_screen/cryptic_pocket_screen_with_availability.csv` -- + drug-availability columns
- `results/fascin_ensemble_rmsd/holo_ligand_set_final.csv` -- final 11-ligand table
- `results/fascin_ensemble_rmsd/holo_ligands_aligned/*.pdb` -- reference poses in the 3LLP frame
- `GATE_DEFINITION_PREREGISTRATION.md`, `APO_SELECTION_RULE_PREREGISTRATION.md`

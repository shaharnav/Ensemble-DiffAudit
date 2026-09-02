# LfrR: does the ConforMix ensemble beat the rigid crystal on docking affinity?

## Verdict

**No.** Across all 6 candidates that docked successfully, the rigid crystal baseline
(2V57, chain A, native-ligand/waters stripped) scores better than every one of the 3
ConforMix conformers, for every candidate, with no exceptions. Best ensemble-vs-baseline
gap ranges from −0.14 to −1.58 kcal/mol (ensemble always worse). This is a clean, boring
result reported plainly: the induced-fit ensemble generated for LfrR does not surface a
more favorable binding pocket than the native structure for these 7 DiffSBDD-generated
candidates.

A second, more consequential finding surfaced along the way: **even the 3 ConforMix
conformers that passed ConforMix's own built-in clash filter were not receptor-preparable
as generated.** Fixing that was necessary before any of the 3 conformers could be docked
at all.

## Pipeline recap

1. **Conformer generation** (Colab, DiffSBDD + ConforMix): 10 drug candidates generated
   against the LfrR pocket (2V57, PRL/proflavine centroid), 7 passed QED≥0.3/SA≤6.0
   filtering. ConforMix ran Boltz twist-guided sampling (targets 0.0/1.0/2.0 Å RMSD,
   `--subset-residues 67-152`, structured-regions-only) → 6 raw conformers → **3 passed
   ConforMix's own CA-CA/C-N/heavy-atom clash filter** (50% rejected at that stage alone).
2. **Alignment**: ConforMix output is in Boltz's own prediction frame, not the crystal's
   coordinate frame (raw, unfit CA distance ~77 Å) — confirmed empirically before trusting
   any docking box coordinates. Residue *number* couldn't be used to pair atoms (ConforMix
   renumbers 1..N from the FASTA; the crystal keeps native numbering with an internal gap
   from unmodeled residues), so alignment was redone by sequence-identity matching
   (pairwise global alignment, matched-identity positions only), not residue number or
   naive index — the same numbering trap already documented for the LasB payload code.
3. **Receptor-prep failure, then fix**: all 21 docking jobs against the 3 aligned
   ConforMix conformers failed at `mk_prepare_receptor` with "Expected N paddings... but
   got 0" — Meeko inferring spurious covalent bonds between residues that are far apart
   in sequence (e.g. 60↔106, 65↔117) because their side-chain atoms sit close enough in
   space to look bonded. Backbone CA-CA geometry was clean (0 anomalous consecutive
   distances), so this was a side-chain-only distortion. Diagnostic minimization energies
   confirmed it was real, not a fluke: starting potential energy 4.08×10⁶ kJ/mol
   (var_0), 1.17×10⁸ (var_1), 9.47×10⁸ (var_2) — all converging to a consistent ≈−18,000
   to −19,000 kJ/mol after a short heavy-atom-restrained OpenMM minimization
   (`minimize_conformix.py`, PDBFixer + amber14-all.xml, k=100 kJ/mol/nm² restraint so the
   minimization relieves clashes without unfolding the ConforMix-predicted conformation).
   Post-minimization, all 3 conformers passed receptor prep cleanly.

   **This is the one methodologically important takeaway**: ConforMix's own clash filter
   (CA-CA / C-N / heavy-atom-distance based) is not sufficient to guarantee a
   receptor-preparable structure. It passed structures with local side-chain geometry bad
   enough to register as ~10⁶–10⁹ kJ/mol above a relaxed minimum. Any future ConforMix
   ensemble intended for rigid-receptor docking should budget for this repair step rather
   than assume "passed ConforMix's filter" implies "dockable."

4. **Docking**: AutoDock Vina 1.2.7, box center (7.28, −5.27, 6.68) — the chain-A PRL
   (proflavine) ligand centroid in 2V57 — 24 Å cube, exhaustiveness 16, num_modes 9,
   seed 42. 7 candidates × 4 receptors (3 minimized ConforMix conformers + rigid crystal
   baseline) = 28 jobs, 27/28 succeeded.

## Results

| Compound | Best affinity (kcal/mol) | Winning receptor | Baseline (crystal) | Best ensemble conformer | Δ (ensemble − baseline) |
|---|---|---|---|---|---|
| Cmpd-0003 | **−8.92** | crystal baseline | −8.92 | −7.70 | −1.22 |
| Cmpd-0009 | −8.64 | crystal baseline | −8.64 | −7.06 | −1.58 |
| Cmpd-0007 | −8.35 | crystal baseline | −8.35 | −6.80 | −1.55 |
| Cmpd-0008 | −7.56 | crystal baseline | −7.56 | −6.49 | −1.07 |
| Cmpd-0006 | −7.13 | crystal baseline | −7.13 | −6.99 | −0.14 |
| Cmpd-0001 | −5.79 | crystal baseline | −5.79 | −4.73 | −1.06 |
| Cmpd-0004 | FAILED (all 4 receptors) | — | — | — | — |

Cmpd-0004's SMILES (`COC1=CC=C[C@]23[C@H]1N2[C@@H]3[C@@H](CC(=O)Cc1ccc2c(F)cccc2c1)C(=O)O`)
has a strained fused tricyclic/aziridine-like ring system that RDKit's ETKDG conformer
embedding cannot resolve — `prepare_ligand`'s documented failure path, independent of and
prior to any receptor issue. It failed identically against all 4 receptors, consistent
with a ligand-side cause, not a receptor-side one.

Full per-receptor affinity matrix: `results_lfrr.csv` / `results_lfrr.json`.

## What this does and doesn't say about LfrR

- The ~4.5–4.6 Å global/pocket CA RMSD motion documented for LfrR's real apo→holo
  transition (`targets.yaml`) is large and structured, but ConforMix/Boltz's twist-guided
  *prediction* of that motion, even after clash-filtering and post-hoc minimization, did
  not produce a pocket geometry more favorable to any of the 6 dockable candidates than
  the native crystal pocket itself.
- This is not evidence the real induced-fit pocket doesn't exist (targets.yaml's own
  audit already showed the real crystallographic motion is large and DSSP-structured,
  not disordered) — it's evidence that *this generation pipeline's prediction* of that
  motion, at these twist-guidance settings, doesn't yield a better docking target than
  the crystal for this ligand set.
- Consistent with the open question already flagged when LfrR was staged: whether
  ConforMix's local-perturbation sampling can reach a domain-scale (if structured) motion
  was explicitly called out as unresolved. This docking result is one data point against
  it being useful here, not a definitive resolution.

## Files

- `dock_lfrr.py` — driver: unpack payload, sequence-align ConforMix conformers to the
  chain-A-only crystal reference, minimize, build the rigid baseline, run the N×M matrix.
- `minimize_conformix.py` — PDBFixer + restrained OpenMM minimization used to repair
  ConforMix output before receptor prep.
- `results_lfrr.json` / `results_lfrr.csv` — full ranked results with per-receptor
  affinity breakdown, QED, SA.
- `pdbs/2V57.pdb` — chain-A-only reference (rebuilt from the full biological assembly;
  the original 4-chain file broke alignment against ConforMix's single-chain output).

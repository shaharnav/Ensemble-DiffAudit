# Ensemble-DiffAudit

Ensemble-DiffAudit generates *de novo* ligands with DiffSBDD, predicts a receptor
conformational ensemble with ConforMix/Boltz, and cross-docks every ligand against every
conformer with AutoDock Vina. The question it was built to answer: does docking against a
predicted conformational ensemble find better binding poses than docking against a single
rigid crystal structure — or does it just look that way because taking the best of several
noisy samples inflates the result on its own?

The pipeline itself (below) was validated on a single pilot target, trypsin (3PTB). Since
then the project has grown into a multi-target program that pre-registers gates before
touching a new target, runs the same generation → audit → noise-correction pipeline, and
reports every outcome — including the negative and terminated ones — without hedging. See
[Project status](#project-status) for what's been run and what each result means, and
[Repository layout](#repository-layout) for how the code is organized.

**[Live demo](https://ensemble-diff-audit.vercel.app/)** — trypsin (PDB 3PTB), DiffSBDD-generated
candidates docked against a 6-conformer breathing receptor ensemble, rendered as an interactive
4D denoising-trajectory viewer. Static build, no backend required (see
[Frontend & Live Demo](#frontend--live-demo) below).

## Project status

| Target | Stage | Outcome | Docs |
|---|---|---|---|
| Trypsin (3PTB) — pilot | Full pipeline run, N=112 | **Negative, statistically significant**: ensemble underperforms noise-matched rigid baseline | [Key findings](#key-findings) below |
| LasB (*P. aeruginosa* elastase) | Phases 4-6 generation + docking | **Null result** (n=15 ligands, CI spans zero). Earlier phase 4-6 numbers were computed on zinc-free receptors and are superseded — see `lasb_ensemble_rmsd_final_report.md` | `lasb_*.py`, `lasb_ensemble_rmsd_final_report.md` |
| Fascin | Step 1-3 generation | **Terminated/failed**: ConforMix twist-guidance produced 0/4 passing conformers at every tested setting | `fascin_archive/`, `results/fascin_ensemble_rmsd/STEP3_FINAL_FINDING.md` |
| LfrR (TetR-family regulator) + proflavine | Pre-screen | **Viable** — only target to clear all 6 pre-registered gates across ~30 candidates in 3 screening rounds. Staged, generation not yet run | `NEW_TARGET_SCREEN_RESULTS.md`, `targets.yaml` |
| Xylellain, CLas IMPDH, AQP5, DYRK1A | Pre-screen | **Rejected** before any GPU time — DSSP/pocket-RMSD gates failed | `xylellain_archive/`, `impdh_archive/`, `PREREGISTRATION_EXP3.md` |
| CauDHFR (*Candida auris* DHFR) | Experiment 3, gate sequence | In progress | `PREREGISTRATION_EXP3.md`, `caudhfr_*.py` |
| GluR | Pocket/box determination | In progress | `glur_*.py` |

Pre-registration is the load-bearing discipline here: gates (DSSP structuredness, pocket RMSD,
holo ligand druglikeness, cofactor differences, apo-selection) are committed to a file
*before* a target is screened, so a rejection can't be quietly redefined into a pass after
the fact. A target that fails pre-screen gets its working files moved to `<target>_archive/`
and stays documented, not deleted.

## Key findings

Full target: trypsin (PDB 3PTB). N=112 DiffSBDD-generated candidates, cross-docked against a
6-conformer ConforMix/Boltz ensemble plus the rigid crystal structure (784 docking jobs), then
re-docked against the rigid crystal alone across 6 seeds to establish a noise floor (672 more
jobs). 108/112 candidates had complete data across every comparison; 4 failed RDKit conformer
embedding and are excluded from all statistics below (see [Limitations](#limitations)).

**The noise-corrected comparison is the one that matters, and it comes out negative.**
Taking the best of 6 conformer scores inflates the result even with zero real conformational
change, purely from maximizing over 6 noisy Vina samples. To measure that inflation, we docked
each ligand against the *unmodified* crystal structure 6 more times with different seeds and
took the best of those — that's the score pure sampling noise produces on its own. Any real
induced-fit benefit has to beat it.

| Metric | Value |
|---|---|
| Mean `noise_corrected_delta` | **-0.060 kcal/mol** (median -0.060), N=108 |
| 95% bootstrap CI (10,000 resamples) | **[-0.114, -0.005]** — entirely negative, excludes zero |
| Paired t-test | t=-2.146, p=0.034 |
| Wilcoxon signed-rank | p=0.014 |
| Cohen's d (paired) | -0.207 |
| Candidates beating the noise-matched rigid baseline | 42/108 (39%) |
| Candidates where the crystal beat every conformer outright | 44/108 (41%) |

At N=108, this isn't "no detectable benefit" — it's statistically significant evidence that the
predicted ensemble performs *worse* than a rigid structure re-sampled with different seeds.
Whatever apparent gain shows up when comparing ensemble-best-of-6 against a single crystal run
(mean `delta_ensemble_vs_crystal` = +0.074 kcal/mol, N=108, p=0.023) is consistent with sampling
inflation, not induced fit.

Two things this doesn't rule out: a real effect too small to detect at N=108 (the CI's upper
bound is only -0.005, close to zero), or a real effect specific to a target more flexible than
trypsin (see [Limitations](#limitations)). It does rule out treating the raw ensemble-vs-crystal
comparison as evidence of induced fit on its own.

No raw affinity number here should be read without ligand efficiency alongside it — Vina score
scales with molecule size, so a larger molecule can score better purely by having more atoms to
contact the pocket, independent of binding quality. Ligand efficiency (|affinity| / heavy atoms)
across the set: mean 0.282, median 0.271, sd 0.083 (N=108). Best raw score: Cmpd-0015 at
-8.629 kcal/mol, QED 0.494, SA 4.977, winning conformer `conformix_var_3.pdb`.

Pocket-lining conformational displacement across the 6 conformers (22 residues within 8 Å of
the native ligand centroid): CA RMSD 0.20-0.27 Å. Spearman correlation between this
displacement and `delta_ensemble_vs_crystal`: rho=-0.122, p=0.210, N=108 — no significant
relationship, consistent with the pocket geometry not actually moving enough to produce a real
docking effect.

Every number above is reproducible from `results.json` and `rigid_control.csv` via
`python summarize_results.py`, which regenerates `results_summary.md`.

## Validation

Before trusting any candidate affinity, the docking engine itself is validated by redocking
the native co-crystallized ligand (benzamidine, into 3PTB's own stripped receptor) and checking
that the top-scoring pose recovers the crystallographic pose — not just that the score looks
plausible, since a wrong pose can still score well.

| Metric | Value |
|---|---|
| Median top-pose RMSD vs. crystal pose (6 seeds) | 0.393 Å |
| Pass rate at 2.0 Å threshold | 100% (6/6) |
| RMSD method | Symmetry-corrected heavy-atom RMSD (RDKit `CalcRMS`), no post-hoc superposition |

Run it yourself: `./venv/bin/python calibrate.py` (writes `validation/redocking_report.json`).

## Limitations

- **Vina scoring error.** AutoDock Vina's scoring function has an inherent error against
  experimental binding affinity on the order of ±2 kcal/mol. Differences smaller than that
  between candidates, or between ensemble and rigid runs, are within the noise floor of the
  method itself — this is part of why the noise-corrected control above matters.
- **Max-over-N inflation.** Selecting the best score across N samples (conformers, or seeded
  reruns) is guaranteed to look at least as good as, and usually better than, any single sample
  — even with zero real signal. The rigid-seed control exists specifically to measure and
  subtract this inflation; the raw ensemble-vs-crystal comparison in isolation is not valid
  evidence of induced fit.
- **Target rigidity.** 3PTB trypsin is a small, well-ordered serine protease with a shallow,
  rigid active site — close to a worst case for demonstrating induced fit. Measured pocket CA
  displacement across the predicted ensemble is only 0.20-0.27 Å, near the resolution limit of
  what crystallographic coordinates can distinguish from noise. A target with documented
  conformational plasticity (e.g. a kinase with DFG-in/DFG-out states, or a protein with a known
  cryptic pocket) would be a more informative test of the pipeline's actual claim, and the
  negative result here may not generalize to one.
- **Score-size dependence.** Raw Vina affinity scales with molecule size; a larger ligand scores
  better partly just from having more atoms to contact the pocket. Ligand efficiency
  (|affinity| / heavy atoms) is reported alongside every affinity comparison in this README for
  that reason — raw score alone is not a fair ranking across candidates of different sizes.
- **Attrition is only partially explained.** Of candidates generated in this run, some are
  dropped before reaching a valid SMILES (RDKit sanitization/QED/SA filtering) or before
  completing docking (RDKit conformer embedding failure in `rigid_control.py`, affecting 4 of
  112 final candidates). `attrition.py` diagnoses geometric non-convergence where it can, but
  explicitly does not claim to reproduce every original generation-time failure reason —
  a naive local bond-order reconstruction was tested and rejected because it produced false
  positives even on candidates confirmed valid.

## How it works

### 1. Generation (Colab, GPU)
[`diffusion_model/diffsbdd_generation.ipynb`](diffusion_model/diffsbdd_generation.ipynb) runs
pocket-conditioned generation with [DiffSBDD](https://github.com/arneschneuing/DiffSBDD),
filters candidates by RDKit sanitization / QED / synthetic accessibility, and runs
[ConforMix](https://github.com/drorlab/conformix) (built on [Boltz](https://github.com/jwohlwend/boltz))
to predict a small conformational ensemble of the receptor via guided-RMSD sampling. Output is
packaged into `ensemble_payload.zip`: valid candidates (SDF), their 4D denoising trajectories,
the receptor conformers, and a `metadata.json` with the pocket center/radius used.

`POCKET_CENTER` matters more than anything else in this step — a mis-centered pocket conditions
generation on the wrong (often flat, non-binding) surface patch and caps every downstream
affinity regardless of docking-engine tuning. Get it from the true binding site (e.g. the
centroid of a co-crystallized ligand's HETATM coordinates), not a guess.

### 2. Local audit (AutoDock Vina)
`ensemble_auditor.py` unpacks the payload, structurally aligns every conformer to the full
reference crystal structure (sequence-index CA pairing, not residue-number pairing — ConforMix
renumbers residues 1..N), builds a real rigid-crystal baseline (waters and the native ligand
stripped, not the small pocket crop used only for generation conditioning), and cross-docks
every candidate against every conformer plus the baseline. Docking jobs run concurrently (see
[Performance](#performance)) and are cached per-job on (box, exhaustiveness, seed, receptor
mtime), so an interrupted run resumes without redoing completed work.

Per-candidate output columns (`results.csv`): `crystal_affinity`, `ensemble_best_affinity`,
`ensemble_best_conformer`, `ensemble_mean_affinity`/`sd`/`range`, `overall_best_affinity`,
`delta_ensemble_vs_crystal` (positive = ensemble beat the crystal), plus H-bond count, QED, and
SA score. `True_Affinity`/`Baseline_Affinity`/`Winning_Conformation` are kept as deprecated
aliases for `overall_best_affinity`/`crystal_affinity`/`overall_best_structure`.

### 3. Noise correction (`rigid_control.py`)
Docks each candidate against the unmodified crystal structure across the same number of seeds
as there are conformers, using the exact box/exhaustiveness the ensemble run used (read back
from its own logged params and asserted to match — not re-specified, so the two can't silently
drift). `noise_corrected_delta = rigid_max_over_seeds - ensemble_best_affinity`, positive means
the ensemble beat pure sampling noise.

### 4. Chemistry and pose-quality (`chem_metrics.py`)
Adds ligand efficiency, Lipinski descriptors, structural-alert screening (epoxides, Michael
acceptors, alkyl halides, etc.), and [PoseBusters](https://github.com/maabuu/posebusters)
geometry/clash validation against each candidate's actual winning docked pose.

### 5. Pocket geometry (`pocket_rmsd.py`) and attrition (`attrition.py`)
CA and all-heavy-atom RMSD of pocket-lining residues per conformer, and a funnel accounting for
every generated candidate through validity filtering and docking.

### 6. Statistics (`summarize_results.py`)
Recomputes every number in [Key Findings](#key-findings) fresh from `results.json` and
`rigid_control.csv` and writes `results_summary.md`.

## Performance

Sequential docking (original full-matrix run, N=112, 784 jobs): 183 runs/hour, `--cpu 4` per
job, one job at a time.

Concurrent docking (current default): jobs are independent, so `docking_engine.compute_parallel_plan()`
runs several at once via a thread pool, derived from `os.cpu_count()` at runtime (not hardcoded
— on a 4-core machine it falls back to the old 1-job-at-a-time behavior). On the 10-core machine
used for this run: 4 concurrent jobs at `--cpu 2` each, ~430 runs/hour measured on the 672-job
rigid-control run — about 2.3x the sequential rate.

## Local environment setup

1. **Initialize virtual environment & install requirements**:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```
2. **Download the AutoDock Vina binary**:
   ```bash
   mkdir -p bin
   curl -L https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_mac_aarch64 -o bin/vina_1.2.7_mac_aarch64
   chmod +x bin/vina_1.2.7_mac_aarch64
   ```
   (Swap the asset name for your platform — see the
   [v1.2.7 release assets](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7).)
3. **Calibrate the engine**: verify pose recovery (see [Validation](#validation)) before
   trusting any candidate affinity:
   ```bash
   curl -L https://files.rcsb.org/download/3PTB.pdb -o pdbs/3PTB.pdb
   ./venv/bin/python calibrate.py
   ```

## Usage: M×N ensemble pipeline

### 1. Generation (Colab)
Open `diffusion_model/diffsbdd_generation.ipynb` in a GPU Colab runtime.

1. **Cell 1 — Configuration**: set `PDB_ID`, `POCKET_CENTER` (see [How it works](#how-it-works)),
   `POCKET_RADIUS` (≥12 Å to cover the full catalytic pocket), and `N_SAMPLES`.
2. **ConforMix ensemble interlude**: tune `--num-twist-targets` / `--samples-per-target` (their
   product is the conformer count M) and `--twist-target-stop` (max RMSD opening, kept at 2.0 Å
   to capture loop breathing without unfolding the pocket).
3. Run all cells top-to-bottom. Cell 6 (generation) prints per-batch progress (batch count,
   molecule count, elapsed time, ETA) so a large `N_SAMPLES` run isn't silent for hours.
4. Download `ensemble_payload.zip` (Cell 9 stages it and copies it to Google Drive).

### 2. Local audit
```bash
./venv/bin/python ensemble_auditor.py --payload ensemble_input/ensemble_payload.zip
```
Add `--dry-run` first on a large payload to see the job count and a rough time estimate before
committing to a long run.

### 3. Noise correction, chemistry, geometry, statistics
```bash
./venv/bin/python rigid_control.py
./venv/bin/python chem_metrics.py
./venv/bin/python pocket_rmsd.py
./venv/bin/python summarize_results.py
```

### 4. Interpreting `results.csv`
Look up `Cmpd-XXXX` in `results/payload_unpacked/valid_trajectories/mol_XXXX.xyz` for the
generated 3D coordinates, or `ensemble_best_conformer` in
`results/payload_unpacked/aligned_receptors/` for the receptor structure it docked against.

## Frontend & Live Demo

`frontend/` is a Vite + React app with two views:

- **Trajectory Viewer** — replays the DiffSBDD denoising trajectory for each candidate against
  the 3Dmol-rendered breathing receptor ensemble. Reads a static `viz_bundle.json` (built by
  `build_viz_bundle.py`), so it needs no backend. This is the view deployed to Vercel.
- **Docking** — live 1 PDB : 1 SMILES docking (H-bonds, metal coordination, salt bridges,
  halogen bonds) against `app.py`'s Flask/Vina backend. Requires a local backend and AutoDock
  Vina, so it's disabled in the static deploy.

### Run locally
```bash
cd frontend
npm install
npm run dev       # trajectory viewer only
```
To also use the Docking tab locally, run `python app.py` (with the venv + Vina binary from
setup above) and set `VITE_ENABLE_DOCKING=true` in `frontend/.env.local` — Vite proxies `/api`,
`/results`, and `/pdbs` to the Flask server per `frontend/vite.config.js`.

### Deploy the static demo to Vercel
The trypsin trajectory data is committed under `frontend/public/results/payload_unpacked/` so
the build is fully static.

1. In the Vercel dashboard, import this GitHub repo.
2. Set **Root Directory** to `frontend`.
3. Framework preset: Vite (build command `npm run build`, output `dist` — already set in
   `frontend/vercel.json`).
4. Leave `VITE_ENABLE_DOCKING` unset so the Docking tab stays hidden (no backend on Vercel).
5. Deploy, then swap the live demo link at the top of this README for your `*.vercel.app` URL.

To refresh the demo data after regenerating a new ensemble, copy the new `viz_bundle.json` and
`receptor_breathing.pdb` from `results/payload_unpacked/` into
`frontend/public/results/payload_unpacked/` and redeploy.

## Work in progress

- Extending the Docking tab's active-site analysis to cover pi-pi stacking on top of the
  existing H-bond, metal coordination, salt bridge, and halogen bond detection.
- CauDHFR (Experiment 3) and GluR gate sequences, still running — see
  [Project status](#project-status).
- Staging LfrR/proflavine for generation, the one target so far to clear every
  pre-registered pre-screen gate.

## Repository layout

- **Core pipeline** (target-agnostic): `analyzer.py`, `app.py`, `fetcher.py`,
  `docking_engine.py`, `ensemble_auditor.py`, `chem_metrics.py`, `attrition.py`,
  `calibrate.py`, `pocket_rmsd.py`, `rigid_control.py`, `summarize_results.py`,
  `build_viz_bundle.py`. See [How it works](#how-it-works).
- **Target screening** (finds and gates the next target before any generation run):
  `target_screen.py`, `cryptic_pocket_screen.py`, `exp2_screen_enrich.py`, plus the
  generic `exp3_*.py` gate scripts (DSSP, pocket-RMSD, holo checks) applied to whichever
  candidate `targets.yaml` / `pocket_config_*.yaml` points at. Screening criteria are
  pre-registered in `PREREGISTRATION.md`, `PREREGISTRATION_EXP3.md`,
  `GATE_DEFINITION_PREREGISTRATION.md`, `APO_SELECTION_RULE_PREREGISTRATION.md`.
- **Per-target scripts**: prefixed by target (`lasb_*.py`, `glur_*.py`, `caudhfr_*.py`,
  `hiv_step1bc_screen.py`) — receptor prep, ligand transplant, and RMSD/geometry
  validation specific to that target's structures. `casf_pipeline/` and `casf_*.py` are a
  separate CASF-2016 benchmark reproduction, not a target campaign.
  `phase5_*.py`/`phase6_*.py`/`phase7_*.py` are later-stage decoy/docking/reproducibility
  scripts, currently attached to the LasB campaign.
- **`<target>_archive/`**: working files for a target whose pre-screen failed or whose
  campaign was terminated (`xylellain_archive/`, `impdh_archive/`, `fascin_archive/`).
  Results and writeups for these stay under `results/<target>/`; only the driver scripts
  move. `archive/` (no target prefix) holds old payload zips, not target work.
- **`results/`**: one subdirectory per target/experiment (e.g. `results/fascin_ensemble_rmsd/`,
  `results/lasb_ensemble_rmsd/`, `results/experiment3_glur/`), plus `target_screen/` for
  cross-target screening output.
- **`frontend/`**: the Vite/React trajectory viewer and docking UI — see
  [Frontend & Live Demo](#frontend--live-demo).

## License and upstream projects

This project is licensed under the [MIT License](LICENSE).

It builds on, and does not modify the license terms of:
- [DiffSBDD](https://github.com/arneschneuing/DiffSBDD) (MIT)
- [ConforMix](https://github.com/drorlab/conformix) (MIT) and [Boltz](https://github.com/jwohlwend/boltz) (MIT)
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) (Apache-2.0)
- [PoseBusters](https://github.com/maabuu/posebusters), [RDKit](https://github.com/rdkit/rdkit), [Meeko](https://github.com/forlilab/Meeko)

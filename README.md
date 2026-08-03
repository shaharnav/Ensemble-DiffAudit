# Ensemble-DiffSplat

Ensemble-DiffSplat is a high-throughput computational biology pipeline integrating Generative Structure-Based Drug Design (DiffSBDD) with dynamic protein conformational breathing (ConforMix & Boltz). It systematically generates *de novo* ligands and evaluates them against an induced-fit structural ensemble via AutoDock Vina—saving robust cross-docking matrices for dynamic 4D rendering via Gaussian Splatting interfaces.

**[Live demo](https://your-vercel-project.vercel.app)** — trypsin (PDB 3PTB), 3 DiffSBDD-generated candidates docked against a 6-conformer breathing receptor ensemble. Static build, no backend required (see [`frontend/README` setup below](#frontend--live-demo)).

## Local Environment Setup
Before executing the auditor locally, you must hydrate the Python environment and download the AutoDock Vina binary solver.

1. **Initialize Virtual Environment & Install Requirements**:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```
2. **Download AutoDock Vina Binary**:
   Ensure you have the correct executable for your architecture mapped locally so the engine can run the solver matrix:
   ```bash
   mkdir -p bin
   # Example: Downloading macOS ARM binary
   curl -L https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_mac_catalina_arm64 -o bin/vina_1.2.7_mac_aarch64
   chmod +x bin/vina_1.2.7_mac_aarch64
   ```

## Usage: $M \times N$ Ensemble Pipeline

### 1. Cloud Generation (Colab)
To generate your candidate libraries and protein structural states, open the `diffsbdd_generation.ipynb` notebook in a Google Colab GPU-runtime instance.

1. **Configure Parameters (Cell 1):** Set your target `PDB_ID`, select your critical `subset_residues`, and define `N_SAMPLES` to control how many *de novo* molecular geometries DiffSBDD synthesizes.
2. **Tune ConforMix Flexibility:** Scroll down to the `run_conformixrmsd_boltz` cell. 
   - `--num-twist-targets`: Controls the quantity of unique structural variations (e.g., set to 5 for a robust ensemble).
   - `--twist-target-stop`: Controls the maximal root-mean-square deviation (RMSD) opening limit (highly recommended to stay at `2.0` Å to capture natural "loop breathing" without breaking the pocket).
3. **Execute:** Run all cells top-to-bottom. The pipeline will pull crystal structures, generate the matrices, structurally align all variants to the native active-site origin, and package the payload.
4. **Download Payload:** Once execution completes, download `ensemble_payload.zip` locally from the Colab filesystem directory.

### 2. Local Auditor Docking
Your local machine manages the thermodynamic evaluation natively using AutoDock Vina, fully backing up results to disk.

1. **Stage Payload:** Place `ensemble_payload.zip` directly into the `ensemble_input/` directory of this local repository.
2. **Execute Auditor:** Run the following command exactly:
   ```bash
   ./venv/bin/python ensemble_auditor.py --payload ensemble_input/ensemble_payload.zip
   ```
   *The auditor securely checks Vina cache points, auto-constructs ID trackers for SMILES linkage, and parses all states.*

### 3. Interpreting Results
When computation halts, the original payload safely relocates into the timestamp-sealed `archive/` folder. Analyzing the outputs is trivial:

- **`results.csv`**: Contains your full combinatorial matrix. Drag this into any DataFrame to track binding variance or calculate aggregate stability across the protein's conformational sweep.
- **The CLI Rankings:**
  ```text
  --- Top candidates by True Affinity ---
  Rank   ID              Affinity   Baseline   H-Bonds    Winning Conformation     
  ---------------------------------------------------------------------------------
  1      Cmpd-0083       -4.94      -2.85      4          conformix_var_0.pdb     
  ```
  - **Affinity**: The "induced-fit" energy peak achieved against the flexible states.
  - **Baseline**: What that identical drug scored against the rigid immobile crystal default. Use the differential here to prove binding plasticity!
  - **ID**: Look up `Cmpd-0083` strictly within your `results/payload_unpacked/valid_trajectories/` directory. Fetching `mol_0083.xyz` lets you pipe the perfect generated coordinates into your downstream visual render engines linearly with no ambiguity.

## Frontend & Live Demo

`frontend/` is a Vite + React app with two views:

- **Trajectory Viewer** — replays the DiffSBDD denoising trajectory for each candidate against the 3Dmol-rendered breathing receptor ensemble. Reads a static `viz_bundle.json` (built by `build_viz_bundle.py`), so it needs no backend. This is the view deployed to Vercel.
- **Docking** — live 1 PDB : 1 SMILES docking (H-bonds, metal coordination, salt bridges, halogen bonds) against `app.py`'s Flask/Vina backend. Requires a local backend and AutoDock Vina, so it's disabled in the static deploy.

### Run locally
```bash
cd frontend
npm install
npm run dev       # trajectory viewer only
```
To also use the Docking tab locally, run `python app.py` (with the venv + Vina binary from setup above) and set `VITE_ENABLE_DOCKING=true` in `frontend/.env.local` — Vite proxies `/api`, `/results`, and `/pdbs` to the Flask server per `frontend/vite.config.js`.

### Deploy the static demo to Vercel
The trypsin trajectory data is committed under `frontend/public/results/payload_unpacked/` so the build is fully static.

1. In the Vercel dashboard, import this GitHub repo.
2. Set **Root Directory** to `frontend`.
3. Framework preset: Vite (build command `npm run build`, output `dist` — already set in `frontend/vercel.json`).
4. Leave `VITE_ENABLE_DOCKING` unset so the Docking tab stays hidden (no backend on Vercel).
5. Deploy, then swap the live demo link at the top of this README for your `*.vercel.app` URL.

To refresh the demo data after regenerating a new ensemble, copy the new `viz_bundle.json` and `receptor_breathing.pdb` from `results/payload_unpacked/` into `frontend/public/results/payload_unpacked/` and redeploy.

## Work in progress features
- 4D Gaussian Splatting Visualization: the current Trajectory Viewer uses 3Dmol ball-and-stick/cartoon rendering; a true Gaussian Splatting (smooth electron-cloud) renderer is still planned.
- Extending the Docking tab's active-site analysis to cover pi-pi stacking on top of the existing H-bond, metal coordination, salt bridge, and halogen bond detection.
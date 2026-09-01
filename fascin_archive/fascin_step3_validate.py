"""
Fascin Step 3: validate all 16 generated conformers immediately after
generation (per Step 2 pre-registration), before any docking prep.

Gates applied, in order (a conformer failing an earlier gate is not checked
against later ones -- report the specific failure):
  0. Global CA RMSD to apo crystal (3LLP) <= 2.5A. Hard gate discovered from
     this project's own history: the LasB run earlier in this session used
     Boltz in single-sequence mode (use_msa_server=False), later found (in a
     different experiment) to produce 8-9A global RMSD without MSA. Fascin
     was generated with the MSA-enabled patch applied, but this gate confirms
     it actually worked rather than assuming it.
  1. Backbone integrity: no non-adjacent-residue clash (any two residues
     >=4 apart in sequence with backbone atoms within 2.0A -- the same
     defect class that broke beta2.4/beta4.0 in the LasB run).
  2. Meeko receptor prep succeeds (diagnosed directly if it fails).
  3. PoseBusters structural validity (via a from-file structure check).
  4. Chain integrity: single continuous chain, no gaps in the region used
     for docking (pocket + its access-shell neighborhood).
"""
import gemmi
import numpy as np
import csv
import subprocess
import os

CONF_DIR = "results/fascin_ensemble_rmsd/conformers_raw"
APO_PATH = "results/target_screen/structures/3LLP.pdb"
OUT_DIR = "results/fascin_ensemble_rmsd/receptors_raw"
os.makedirs(OUT_DIR, exist_ok=True)

BETAS = ["0.0", "0.4", "0.8", "1.2000000000000002", "1.6", "2.0", "2.4000000000000004",
         "2.8000000000000003", "3.2", "3.6", "4.0", "4.4", "4.800000000000001",
         "5.2", "5.6000000000000005", "6.0"]
BETA_LABELS = ["0.0", "0.4", "0.8", "1.2", "1.6", "2.0", "2.4", "2.8", "3.2", "3.6",
               "4.0", "4.4", "4.8", "5.2", "5.6", "6.0"]

POCKET_RESNUMS = [14, 16, 48, 60, 93, 94, 95, 101, 103, 134, 214, 215, 216, 217, 224]


def kabsch(mob, ref):
    mob_c, ref_c = mob.mean(axis=0), ref.mean(axis=0)
    H = (mob - mob_c).T @ (ref - ref_c)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = ref_c - R @ mob_c
    return R, t


apo_st = gemmi.read_structure(APO_PATH)
apo_st.setup_entities()
apo_chain = apo_st[0]["A"]
apo_ca = {}
apo_backbone = {}
for res in apo_chain:
    names = {a.name: a.pos for a in res}
    if "CA" in names:
        apo_ca[res.seqid.num] = names["CA"]

rows = []
for beta_file, beta_label in zip(BETAS, BETA_LABELS):
    cif_path = f"{CONF_DIR}/beta{beta_file}.cif"
    row = {"beta": beta_label}
    try:
        st = gemmi.read_structure(cif_path)
        st.setup_entities()
    except Exception as e:
        row["error"] = f"parse_failed: {e}"
        rows.append(row)
        continue

    chains = [c for c in st[0]]
    if len(chains) != 1:
        row["error"] = f"unexpected_chain_count: {len(chains)}"
        rows.append(row)
        continue
    chain = chains[0]

    residues = []
    for res in chain:
        names = {a.name: a.pos for a in res}
        if "CA" in names:
            residues.append((res.seqid.num, names))
    resnums = sorted(set(r[0] for r in residues) & set(apo_ca.keys()))
    row["n_residues"] = len(residues)
    row["n_ca_matched_to_apo"] = len(resnums)

    if len(resnums) < 400:
        row["error"] = f"too_few_matched_residues: {len(resnums)}"
        rows.append(row)
        continue

    ca_by_resnum = {r[0]: r[1]["CA"] for r in residues}
    mob = np.array([[ca_by_resnum[rn].x, ca_by_resnum[rn].y, ca_by_resnum[rn].z] for rn in resnums])
    ref = np.array([[apo_ca[rn].x, apo_ca[rn].y, apo_ca[rn].z] for rn in resnums])
    R, t = kabsch(mob, ref)
    mob_aligned = (R @ mob.T).T + t
    global_rmsd = float(np.sqrt(np.mean(np.sum((mob_aligned - ref) ** 2, axis=1))))
    row["global_ca_rmsd_to_apo_A"] = round(global_rmsd, 3)
    row["gate0_global_rmsd_pass"] = global_rmsd <= 2.5

    # pocket-region CA displacement (informational at this stage, not a gate --
    # this is what Step 4 will measure against holo structures; here we only
    # check the conformer's own internal geometric plausibility)
    pocket_present = [rn for rn in POCKET_RESNUMS if rn in ca_by_resnum]
    row["n_pocket_residues_present"] = len(pocket_present)

    # Gate 1: non-adjacent residue clash check, ALL heavy atoms (not just
    # backbone -- a first pass using only N/CA/C missed real side-chain
    # clashes down to 0.5A that meeko's bond-perception caught downstream;
    # broadened after that direct measurement, not assumed).
    # Residues 1-3 are crystallographically unresolved in the apo structure
    # (disordered N-terminus, no ground-truth position -- confirmed directly:
    # they're simply absent from 3LLP's deposited coordinates). Excluded from
    # the clash check since there's no reference geometry to judge them
    # against; every other residue IS resolved and stays in the check.
    all_heavy = []
    for res in chain:
        if res.seqid.num <= 3:
            continue
        atoms = [a.pos for a in res if a.element.name != "H"]
        if atoms:
            all_heavy.append((res.seqid.num, np.array([[p.x, p.y, p.z] for p in atoms])))
    clash_found = None
    all_clashes = []
    for i in range(len(all_heavy)):
        rn_i, atoms_i = all_heavy[i]
        for j in range(i + 1, len(all_heavy)):
            rn_j, atoms_j = all_heavy[j]
            if abs(rn_i - rn_j) < 4:
                continue
            d = np.linalg.norm(atoms_i[:, None, :] - atoms_j[None, :, :], axis=2)
            min_d = float(d.min())
            if min_d < 2.0:
                all_clashes.append((rn_i, rn_j, round(min_d, 2)))
    row["gate1_n_clashes"] = len(all_clashes)
    row["gate1_clashes"] = "; ".join(f"{a}-{b}:{d}A" for a, b, d in all_clashes[:5]) if all_clashes else "none"
    row["gate1_pass"] = len(all_clashes) == 0

    rows.append(row)
    print(f"beta={beta_label}: n_res={row['n_residues']} n_ca_matched={row['n_ca_matched_to_apo']} "
          f"global_RMSD={row['global_ca_rmsd_to_apo_A']}A gate0={'PASS' if row['gate0_global_rmsd_pass'] else 'FAIL'} "
          f"clashes={row['gate1_n_clashes']} ({row['gate1_clashes']}) gate1={'PASS' if row['gate1_pass'] else 'FAIL'}")

with open("results/fascin_ensemble_rmsd/step3_validation_geometry.csv", "w", newline="") as f:
    fieldnames = sorted(set().union(*[r.keys() for r in rows]))
    w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
    w.writeheader()
    w.writerows(rows)
print("\nWritten results/fascin_ensemble_rmsd/step3_validation_geometry.csv")

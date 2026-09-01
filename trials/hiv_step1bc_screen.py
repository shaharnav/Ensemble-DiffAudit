"""
Step 1b/1c: check inclusion criteria for the 9 approved-inhibitor HIV-1 protease
holo structures, superpose each onto the 1HHP apo frame (CA atoms), and measure
the conformational change (global CA RMSD, pocket CA RMSD, flap tip separation,
DSSP-classified structured-vs-loop fraction of pocket displacement).

Mutations were already confirmed by direct sequence diff against the wild-type
HXB2 protease sequence (see conversation record); recorded here as constants.
"""
import gemmi
import numpy as np
import csv
import pydssp

STRUCTS = {
    "1HXB": {"ligand": "ROC", "drug": "saquinavir", "mutations": "I3V"},
    "1HXW": {"ligand": "RIT", "drug": "ritonavir", "mutations": "S37N"},
    "1HSG": {"ligand": "MK1", "drug": "indinavir", "mutations": "none (wild-type)"},
    "1OHR": {"ligand": "1UN", "drug": "nelfinavir", "mutations": "none (wild-type)"},
    "1HPV": {"ligand": "478", "drug": "amprenavir", "mutations": "none (wild-type)"},
    "2O4S": {"ligand": "AB1", "drug": "lopinavir", "mutations": "Q7K"},
    "2O4K": {"ligand": "DR7", "drug": "atazanavir", "mutations": "Q7K"},
    "2IEN": {"ligand": "017", "drug": "darunavir", "mutations": "Q7K,L33I,L63I,C67A,C95A"},
    "2O4P": {"ligand": "TPV", "drug": "tipranavir", "mutations": "Q7K"},
}
HOLO_DIR = "results/hiv_protease_ensemble_rmsd/holo_structures"
APO_PATH = f"{HOLO_DIR}/1HHP.pdb"

# HIV-1 protease flap tip: Ile50/Ile50' (canonical numbering, both chains of the dimer)
FLAP_TIP_RESNUM = 50


def kabsch(mob, ref):
    mob_c, ref_c = mob.mean(axis=0), ref.mean(axis=0)
    H = (mob - mob_c).T @ (ref - ref_c)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = ref_c - R @ mob_c
    return R, t


def norm_resnum(n):
    """Some entries offset residue numbering per chain/copy: 2IEN's chain B
    continues from chain A (101-199), 2PC0 offsets its single ASU chain by
    1000 (1001-1099). Normalize back to the 1-99 monomer numbering used by
    apo chain A so residues match across chains/entries."""
    if n > 999:
        return n - 1000
    if n > 99:
        return n - 100
    return n


def get_ca_dict(chain):
    """resnum -> CA gemmi.Position, only for residues that have a CA atom."""
    out = {}
    for res in chain:
        atom = next((a for a in res if a.name == "CA"), None)
        if atom is not None:
            out[norm_resnum(res.seqid.num)] = atom.pos
    return out


def get_backbone_coords(chain, resnums):
    """(N,4,3) array of N/CA/C/O per residue in resnums, in order. Skips residues
    missing any of the four backbone atoms (returns the subset actually present)."""
    coords, kept = [], []
    for rn in resnums:
        res = next((r for r in chain if norm_resnum(r.seqid.num) == rn), None)
        if res is None:
            continue
        names = {a.name: a.pos for a in res}
        if all(n in names for n in ("N", "CA", "C", "O")):
            coords.append([[names[n].x, names[n].y, names[n].z] for n in ("N", "CA", "C", "O")])
            kept.append(rn)
    return np.array(coords), kept


# ---- load apo structure (chain A used as the reference monomer) ----
apo_st = gemmi.read_structure(APO_PATH)
apo_st.setup_entities()
apo_chain = apo_st[0]["A"]
apo_ca = get_ca_dict(apo_chain)

# HIV-1 protease is an obligate homodimer with an interdigitated N/C-terminal
# beta-sheet at the dimer interface -- DSSP on an isolated monomer misses those
# H-bonds and undercounts real beta-strand content (verified directly: 1HHP
# monomer-only DSSP gives 54.5% H/E, dimer-aware gives 61.6%, and the flap tip
# residue 50 itself only registers as strand when the partner chain is present).
# Build the true biological dimer via the file's own BIOMT operators and run
# DSSP on both chains together, then keep only chain A1's per-residue calls.
apo_asmb = apo_st.assemblies[0]
apo_dimer = gemmi.make_assembly(apo_asmb, apo_st[0], gemmi.HowToNameCopiedChain.AddNumber, None)
_coordsA, _keptA = get_backbone_coords(apo_dimer["A1"], list(range(1, 100)))
_coordsB, _keptB = get_backbone_coords(apo_dimer["A2"], list(range(1, 100)))
_ss_dimer = pydssp.assign(np.concatenate([_coordsA, _coordsB], axis=0), out_type="c3")
apo_ss_by_resnum = dict(zip(_keptA, _ss_dimer[:len(_coordsA)]))
apo_resnums_sorted = sorted(apo_ca.keys())

rows = []
for pdbid, info in STRUCTS.items():
    ligcode = info["ligand"]
    st = gemmi.read_structure(f"{HOLO_DIR}/{pdbid}.pdb")
    st.setup_entities()
    model = st[0]

    # find the chain carrying the ligand (protein chain, not a separate HETATM chain).
    # Keep only the majority-occupancy altloc (blank or 'A') when the ligand is
    # modeled in multiple alternate conformations, so the reference pose/centroid
    # isn't a nonsensical blend of two distinct positions.
    lig_atoms, lig_chain_name = [], None
    for chain in model:
        for res in chain:
            if res.name == ligcode:
                atoms = [a for a in res if a.altloc in ("", "A")]
                lig_atoms = atoms if atoms else [a for a in res]
                lig_chain_name = chain.name
                break
        if lig_atoms:
            break

    if not lig_atoms:
        rows.append({"pdbid": pdbid, "ligand": ligcode, "error": "ligand_not_found"})
        continue

    occupancies = [a.occ for a in lig_atoms]
    bfactors = [a.b_iso for a in lig_atoms]

    # Superpose CA atoms of the ligand-bearing chain onto apo chain A (by residue number)
    prot_chain = model[lig_chain_name]
    holo_ca = get_ca_dict(prot_chain)
    common_resnums = sorted(set(holo_ca.keys()) & set(apo_ca.keys()))
    mob = np.array([[holo_ca[rn].x, holo_ca[rn].y, holo_ca[rn].z] for rn in common_resnums])
    ref = np.array([[apo_ca[rn].x, apo_ca[rn].y, apo_ca[rn].z] for rn in common_resnums])
    R, t = kabsch(mob, ref)
    mob_aligned = (R @ mob.T).T + t
    global_ca_rmsd = float(np.sqrt(np.mean(np.sum((mob_aligned - ref) ** 2, axis=1))))

    # ligand atoms transformed into the apo frame -> pocket = apo residues with
    # ANY atom within 5A of ANY ligand atom (not centroid -- these ligands are
    # large/elongated, so a per-atom contact shell is needed, not a distance to
    # the geometric mean position).
    lig_coords = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lig_atoms])
    lig_coords_apo_frame = (R @ lig_coords.T).T + t
    apo_all_atoms = [(norm_resnum(res.seqid.num), a) for res in apo_chain for a in res if a.name != "H"]
    pocket_resnums = set()
    for rn, a in apo_all_atoms:
        p = np.array([a.pos.x, a.pos.y, a.pos.z])
        if np.min(np.linalg.norm(lig_coords_apo_frame - p, axis=1)) < 5.0:
            pocket_resnums.add(rn)
    pocket_resnums = sorted(pocket_resnums & set(common_resnums))

    pocket_idx = [common_resnums.index(rn) for rn in pocket_resnums]
    pocket_rmsd = float(np.sqrt(np.mean(np.sum(
        (mob_aligned[pocket_idx] - ref[pocket_idx]) ** 2, axis=1)))) if pocket_idx else None

    # flap tip (Ile50) displacement -- both chains if dimer present in this entry
    flap_disp = {}
    for cn in model:
        res50 = next((r for r in cn if norm_resnum(r.seqid.num) == FLAP_TIP_RESNUM), None)
        if res50 is None:
            continue
        ca50 = next((a for a in res50 if a.name == "CA"), None)
        if ca50 is None:
            continue
        p = np.array([ca50.pos.x, ca50.pos.y, ca50.pos.z])
        p_aligned = R @ p + t if cn.name == lig_chain_name else p  # only transform the aligned chain properly
        flap_disp[cn.name] = p_aligned

    apo_res50 = next((r for r in apo_chain if norm_resnum(r.seqid.num) == FLAP_TIP_RESNUM), None)
    apo_ca50 = next((a for a in apo_res50 if a.name == "CA"), None) if apo_res50 else None
    flap_tip_disp_A = None
    if lig_chain_name in flap_disp and apo_ca50 is not None:
        apo_p = np.array([apo_ca50.pos.x, apo_ca50.pos.y, apo_ca50.pos.z])
        flap_tip_disp_A = float(np.linalg.norm(flap_disp[lig_chain_name] - apo_p))

    # classify pocket residues as structured (H/E) vs loop using dimer-aware DSSP
    structured_disp, total_disp = 0.0, 0.0
    for rn in pocket_resnums:
        idx = common_resnums.index(rn)
        d = np.linalg.norm(mob_aligned[idx] - ref[idx])
        total_disp += d
        if apo_ss_by_resnum.get(rn, "-") in ("H", "E"):
            structured_disp += d
    frac_structured = (structured_disp / total_disp) if total_disp > 0 else None

    rows.append({
        "pdbid": pdbid, "ligand": ligcode, "drug": info["drug"], "mutations": info["mutations"],
        "n_ca_common": len(common_resnums),
        "global_ca_rmsd_A": round(global_ca_rmsd, 3),
        "n_pocket_residues": len(pocket_resnums),
        "pocket_ca_rmsd_A": round(pocket_rmsd, 3) if pocket_rmsd is not None else "",
        "flap_tip_disp_A": round(flap_tip_disp_A, 3) if flap_tip_disp_A is not None else "",
        "frac_pocket_disp_structured": round(frac_structured, 3) if frac_structured is not None else "",
        "n_lig_atoms": len(lig_atoms),
        "mean_occupancy": round(float(np.mean(occupancies)), 2),
        "min_occupancy": round(float(np.min(occupancies)), 2),
        "mean_bfactor": round(float(np.mean(bfactors)), 1),
        "max_bfactor": round(float(np.max(bfactors)), 1),
    })
    print(f"{pdbid} ({info['drug']}, {ligcode}): global CA RMSD={global_ca_rmsd:.2f}A, "
          f"pocket CA RMSD={rows[-1]['pocket_ca_rmsd_A']}A, flap tip disp={rows[-1]['flap_tip_disp_A']}A, "
          f"frac structured={rows[-1]['frac_pocket_disp_structured']}, "
          f"occ[min/mean]={rows[-1]['min_occupancy']}/{rows[-1]['mean_occupancy']}, "
          f"B[mean/max]={rows[-1]['mean_bfactor']}/{rows[-1]['max_bfactor']}")

with open("results/hiv_protease_ensemble_rmsd/step1bc_screen_1HHP.csv", "w", newline="") as f:
    fieldnames = sorted(set().union(*[r.keys() for r in rows]))
    w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
    w.writeheader()
    w.writerows(rows)
print("\nWritten results/hiv_protease_ensemble_rmsd/step1bc_screen.csv")

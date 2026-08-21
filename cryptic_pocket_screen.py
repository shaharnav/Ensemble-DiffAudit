"""
Generalized apo/holo pocket-change screen, applied to the PocketMiner
validation+test set (49 curated cryptic-pocket pairs), using the gate
definition pre-registered in GATE_DEFINITION_PREREGISTRATION.md:

- contact-shell pocket residues defined in the HOLO (bound) frame, not by
  checking whether the apo copy of a residue happens to be near the aligned
  ligand -- that structurally excludes exactly the residues that move most.
- access-shell = contact-shell +/- 2 sequence positions (same chain).
- apo/holo residue correspondence via pairwise sequence alignment (not
  assumed-identical residue numbering -- unlike LasB/HIV, these pairs are not
  guaranteed to share numbering).
- DSSP on the deposited biological assembly (via BIOMT), not an isolated chain.
- primary gate: access-shell pocket CA RMSD >= 2.0 A
- secondary gate: >= 60% of that displacement is DSSP-structured (H/E)
"""
import gemmi
import numpy as np
import csv
import os
import warnings
from Bio import Align
from Bio.PDB.Polypeptide import protein_letters_3to1_extended
import pydssp

STRUCT_DIR = "results/target_screen/structures"
OUT_CSV = "results/target_screen/cryptic_pocket_screen_results.csv"
AA3TO1 = {k.upper(): v for k, v in protein_letters_3to1_extended.items()}


def kabsch(mob, ref):
    mob_c, ref_c = mob.mean(axis=0), ref.mean(axis=0)
    H = (mob - mob_c).T @ (ref - ref_c)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = ref_c - R @ mob_c
    return R, t


def get_chain_residues(chain):
    """List of (resnum, one_letter_aa, res_object) for standard amino acids
    with a CA atom, in chain order."""
    out = []
    for res in chain:
        aa = AA3TO1.get(res.name)
        if aa is None or aa == "X":
            continue
        ca = next((a for a in res if a.name == "CA"), None)
        if ca is None:
            continue
        out.append((res.seqid.num, aa, res))
    return out


ALIGNER = Align.PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.open_gap_score = -10
ALIGNER.extend_gap_score = -0.5
ALIGNER.match_score = 2
ALIGNER.mismatch_score = -1


def align_residues(apo_res, holo_res):
    """Sequence-align apo/holo chain residue lists; return list of
    (apo_resnum, apo_res_obj, holo_resnum, holo_res_obj) for columns where
    both sides are present AND identical amino acid (safe correspondence,
    robust to differing numbering/point mutations elsewhere)."""
    apo_seq = "".join(r[1] for r in apo_res)
    holo_seq = "".join(r[1] for r in holo_res)
    aln = ALIGNER.align(apo_seq, holo_seq)[0]
    matched = []
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for i in range(a1 - a0):
            ai, bi = a0 + i, b0 + i
            if apo_seq[ai] == holo_seq[bi]:
                matched.append((apo_res[ai][0], apo_res[ai][2], holo_res[bi][0], holo_res[bi][2]))
    return matched


def get_backbone_coords(residues_by_resnum, resnums):
    coords, kept = [], []
    for rn in resnums:
        res = residues_by_resnum.get(rn)
        if res is None:
            continue
        names = {a.name: a.pos for a in res}
        if all(n in names for n in ("N", "CA", "C", "O")):
            coords.append([[names[n].x, names[n].y, names[n].z] for n in ("N", "CA", "C", "O")])
            kept.append(rn)
    return np.array(coords), kept


def dssp_ss_for_chain(structure, chain_name):
    """SS per residue number for chain `chain_name`, computed on the full
    deposited biological assembly (BIOMT-expanded) if available and small
    enough, else on the asymmetric unit as-is (already may contain the
    biological unit for monomeric proteins)."""
    model = structure[0]
    try:
        asmb = structure.assemblies[0]
        full = gemmi.make_assembly(asmb, model, gemmi.HowToNameCopiedChain.AddNumber, None)
    except Exception:
        full = model

    all_coords, chain_target_kept = [], None
    offset = 0
    for chain in full:
        residues = get_chain_residues(chain)
        resnums = [r[0] for r in residues]
        res_by_num = {r[0]: r[2] for r in residues}
        coords, kept = get_backbone_coords(res_by_num, resnums)
        if len(coords) == 0:
            continue
        # identify which expanded copy corresponds to the target original chain
        # (make_assembly names copies "<orig><n>", e.g. "A1"); take the first
        # copy whose name starts with the target chain id as the one we report.
        if chain_target_kept is None and chain.name.rstrip("0123456789") == chain_name:
            chain_target_kept = (offset, offset + len(coords), kept)
        all_coords.append(coords)
        offset += len(coords)

    if not all_coords:
        return {}
    if chain_target_kept is None:
        # fallback: target chain by exact name match (no assembly expansion happened)
        residues = get_chain_residues(model[chain_name])
        resnums = [r[0] for r in residues]
        res_by_num = {r[0]: r[2] for r in residues}
        coords, kept = get_backbone_coords(res_by_num, resnums)
        if len(coords) == 0:
            return {}
        ss = pydssp.assign(coords, out_type="c3")
        return dict(zip(kept, ss))

    full_coords = np.concatenate(all_coords, axis=0)
    if len(full_coords) > 3000:
        # too large to be worth expanding (e.g. huge multimeric assembly) --
        # fall back to just the single target chain, still better than nothing
        start, end, kept = chain_target_kept
        ss = pydssp.assign(full_coords[start:end], out_type="c3")
        return dict(zip(kept, ss))

    ss_all = pydssp.assign(full_coords, out_type="c3")
    start, end, kept = chain_target_kept
    return dict(zip(kept, ss_all[start:end]))


def screen_pair(apo_pdb, apo_chain_id, holo_pdb, holo_chain_id, ligcode):
    apo_path = f"{STRUCT_DIR}/{apo_pdb}.pdb"
    holo_path = f"{STRUCT_DIR}/{holo_pdb}.pdb"
    if not (os.path.exists(apo_path) and os.path.exists(holo_path)):
        return {"error": "missing_pdb_file"}

    apo_st = gemmi.read_structure(apo_path)
    apo_st.setup_entities()
    holo_st = gemmi.read_structure(holo_path)
    holo_st.setup_entities()

    if apo_chain_id not in [c.name for c in apo_st[0]] or holo_chain_id not in [c.name for c in holo_st[0]]:
        return {"error": "chain_not_found"}

    apo_chain = apo_st[0][apo_chain_id]
    holo_chain = holo_st[0][holo_chain_id]

    apo_res = get_chain_residues(apo_chain)
    holo_res = get_chain_residues(holo_chain)
    if len(apo_res) < 20 or len(holo_res) < 20:
        return {"error": "chain_too_short"}

    matched = align_residues(apo_res, holo_res)
    if len(matched) < 30:
        return {"error": "too_few_matched_residues", "n_matched": len(matched)}

    apo_ca = np.array([[next(a for a in r[1] if a.name == "CA").pos.x,
                         next(a for a in r[1] if a.name == "CA").pos.y,
                         next(a for a in r[1] if a.name == "CA").pos.z] for r in matched])
    holo_ca = np.array([[next(a for a in r[3] if a.name == "CA").pos.x,
                          next(a for a in r[3] if a.name == "CA").pos.y,
                          next(a for a in r[3] if a.name == "CA").pos.z] for r in matched])
    apo_resnums = [r[0] for r in matched]

    R, t = kabsch(holo_ca, apo_ca)  # mobile=holo -> ref=apo
    holo_ca_aligned = (R @ holo_ca.T).T + t
    global_rmsd = float(np.sqrt(np.mean(np.sum((holo_ca_aligned - apo_ca) ** 2, axis=1))))

    # find ligand atoms in holo (majority altloc only)
    lig_atoms = []
    lig_search_chains = [holo_chain] + [c for c in holo_st[0] if c.name != holo_chain_id]
    for chain in lig_search_chains:
        for res in chain:
            if res.name == ligcode:
                atoms = [a for a in res if a.altloc in ("", "A")]
                lig_atoms = atoms if atoms else list(res)
                break
        if lig_atoms:
            break
    if not lig_atoms:
        return {"error": "ligand_not_found", "n_matched": len(matched), "global_rmsd_A": round(global_rmsd, 3)}

    lig_coords_holo = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lig_atoms])
    occupancies = [a.occ for a in lig_atoms]
    bfactors = [a.b_iso for a in lig_atoms]

    # contact-shell: holo residues (among matched pairs) with an atom within
    # 5A of a ligand atom, IN THE NATIVE HOLO FRAME (frame-correct, per gate doc)
    contact_idx = []
    for i, r in enumerate(matched):
        holo_res_obj = r[3]
        rc = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in holo_res_obj if a.name != "H"])
        if rc.size and np.min(np.linalg.norm(lig_coords_holo[:, None, :] - rc[None, :, :], axis=2)) < 5.0:
            contact_idx.append(i)
    if not contact_idx:
        return {"error": "no_contact_residues", "n_matched": len(matched), "global_rmsd_A": round(global_rmsd, 3)}

    access_idx = sorted(set(contact_idx) | set(
        j for i in contact_idx for j in range(max(0, i - 2), min(len(matched), i + 3))
    ))

    pocket_rmsd_contact = float(np.sqrt(np.mean(np.sum(
        (holo_ca_aligned[contact_idx] - apo_ca[contact_idx]) ** 2, axis=1))))
    pocket_rmsd_access = float(np.sqrt(np.mean(np.sum(
        (holo_ca_aligned[access_idx] - apo_ca[access_idx]) ** 2, axis=1))))

    ss_by_resnum = dssp_ss_for_chain(apo_st, apo_chain_id)
    structured_disp, total_disp = 0.0, 0.0
    for i in access_idx:
        d = np.linalg.norm(holo_ca_aligned[i] - apo_ca[i])
        total_disp += d
        if ss_by_resnum.get(apo_resnums[i], "-") in ("H", "E"):
            structured_disp += d
    frac_structured = (structured_disp / total_disp) if total_disp > 0 else None

    return {
        "n_matched_ca": len(matched),
        "global_ca_rmsd_A": round(global_rmsd, 3),
        "n_contact_residues": len(contact_idx),
        "n_access_residues": len(access_idx),
        "pocket_ca_rmsd_contact_A": round(pocket_rmsd_contact, 3),
        "pocket_ca_rmsd_access_A": round(pocket_rmsd_access, 3),
        "frac_access_disp_structured": round(frac_structured, 3) if frac_structured is not None else "",
        "n_lig_atoms": len(lig_atoms),
        "mean_lig_occupancy": round(float(np.mean(occupancies)), 2),
        "min_lig_occupancy": round(float(np.min(occupancies)), 2),
        "mean_lig_bfactor": round(float(np.mean(bfactors)), 1),
        "gate_primary_pass": pocket_rmsd_access >= 2.0,
        "gate_secondary_pass": (frac_structured is not None and frac_structured >= 0.6),
    }


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    with open("results/target_screen/pocketminer_pairs.csv") as f:
        pairs = list(csv.DictReader(f))

    rows = []
    for p in pairs:
        ligcode = p["ligand"].split(",")[0].lstrip("#0123456789x")  # strip prefix annotations like '#2' or '2x'
        try:
            result = screen_pair(p["apo_pdb"], p["apo_chain"], p["holo_pdb"], p["holo_chain"], ligcode)
        except Exception as e:
            result = {"error": f"exception: {e}"}
        row = {"apo_pdb": p["apo_pdb"], "holo_pdb": p["holo_pdb"], "ligand": ligcode,
               "motion_type": p["motion_type"], "set": p["set"]}
        row.update(result)
        rows.append(row)
        status = row.get("error", f"pocketRMSD(access)={row.get('pocket_ca_rmsd_access_A')} "
                                    f"frac_struct={row.get('frac_access_disp_structured')}")
        print(f"{p['apo_pdb']}->{p['holo_pdb']} ({ligcode}): {status}")

    fieldnames = sorted(set().union(*[r.keys() for r in rows]))
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWritten {OUT_CSV}")

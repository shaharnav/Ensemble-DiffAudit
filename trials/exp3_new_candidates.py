"""
Pair selection + full gate screen for LdG6PD, LdNMT, MtAlaDH.
For each: HETATM audit, cofactor check, DSSP, pocket CA RMSD.
"""
import os, requests, time
import numpy as np
from Bio.PDB import PDBParser, Superimposer

PDB_DIR = "pdbs/exp3_new"
os.makedirs(PDB_DIR, exist_ok=True)

BUFFER_CODES = {
    "SO4","PO4","GOL","EDO","PEG","MPD","BME","DTT","CIT","ACT","ACE","TRS","MES",
    "HEP","TAR","TLA","FMT","IMD","EPE","PGE","PG4","1PE","15P","EOH","ETF","IOD",
    "CL","NA","K","MG","CA","ZN","MN","FE","CU","NI","CO","CD","HG","PT","AU",
    "BR","I","F","IPA","DMS","GLC","FRU","GAL","MAN","XYL","RIB","FUC","NAG",
    "NDG","BMA","MMA","SEP","PTR","TPO","CSO","HOH","WAT","DOD","ACY","NH4"
}

def download_pdb(pdb_id):
    path = os.path.join(PDB_DIR, f"{pdb_id.upper()}.pdb")
    if not os.path.exists(path):
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", timeout=30)
        r.raise_for_status()
        with open(path, "w") as fh:
            fh.write(r.text)
        time.sleep(0.2)
    return path

def get_hetatms(pdb_path, chain_id=None):
    """Return dict of {res_name: count} for non-water HETATM residues."""
    from collections import Counter
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                c = line[21]
                if chain_id and c != chain_id:
                    continue
                res_name = line[17:20].strip()
                res_seq  = line[22:26].strip()
                if res_name not in ("HOH","WAT","DOD"):
                    seen.add((c, res_seq, res_name))
    return Counter(r for _, _, r in seen)

def parse_sse(pdb_path):
    helix_res, sheet_res = set(), set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HELIX "):
                c = line[19]
                try:
                    s, e = int(line[21:25]), int(line[33:37])
                    for r in range(s, e + 1): helix_res.add((c, r))
                except ValueError: pass
            elif line.startswith("SHEET "):
                c = line[21]
                try:
                    s, e = int(line[22:26]), int(line[33:37])
                    for r in range(s, e + 1): sheet_res.add((c, r))
                except ValueError: pass
    return helix_res, sheet_res

def run_gates(label, apo_id, apo_chain, holo_id, holo_chain, lig_codes,
              cofactor_note=""):
    print(f"\n{'='*65}")
    print(f"{label}")
    print(f"  Pair: {apo_id} (apo, chain {apo_chain})  vs  {holo_id} (holo, chain {holo_chain})")
    if cofactor_note:
        print(f"  Cofactor note: {cofactor_note}")
    print(f"{'='*65}")

    apo_path  = download_pdb(apo_id)
    holo_path = download_pdb(holo_id)

    # ── HETATM audit ──
    apo_hets  = get_hetatms(apo_path,  apo_chain)
    holo_hets = get_hetatms(holo_path, holo_chain)
    apo_nonbuf  = {k: v for k, v in apo_hets.items()  if k not in BUFFER_CODES}
    holo_nonbuf = {k: v for k, v in holo_hets.items() if k not in BUFFER_CODES}
    print(f"\n  APO  non-buffer: {dict(apo_nonbuf)}")
    print(f"  HOLO non-buffer: {dict(holo_nonbuf)}")

    # Cofactor difference check: apo and holo non-buffer minus the drug ligands
    drug_set = set(lig_codes)
    apo_cofac  = {k: v for k, v in apo_nonbuf.items()  if k not in drug_set}
    holo_cofac = {k: v for k, v in holo_nonbuf.items() if k not in drug_set}
    cofactor_diff = set(apo_cofac.keys()) != set(holo_cofac.keys())
    print(f"  Cofactor diff (non-drug non-buffer): apo={set(apo_cofac.keys())}  holo={set(holo_cofac.keys())}")
    print(f"  Cofactor gate: {'FAIL (different cofactors)' if cofactor_diff else 'PASS'}")

    # ── DSSP gate ──
    helix_res, sheet_res = parse_sse(holo_path)
    parser = PDBParser(QUIET=True)
    holo_struct = parser.get_structure(holo_id, holo_path)[0]

    try:
        chain_obj = holo_struct[holo_chain]
    except KeyError:
        chain_obj = list(holo_struct.get_chains())[0]
        holo_chain = chain_obj.id
        print(f"  WARNING: using chain {holo_chain}")

    # Ligand heavy atoms (search all chains if needed)
    lig_atoms = []
    for ch in holo_struct.get_chains():
        for res in ch.get_residues():
            if res.get_resname().strip() in lig_codes and res.get_id()[0] != " ":
                for atom in res.get_atoms():
                    if atom.element and atom.element != "H":
                        lig_atoms.append(atom.get_vector().get_array())
    if not lig_atoms:
        print(f"  ERROR: ligand {lig_codes} not found in {holo_id}")
        return None

    lig_coords = np.array(lig_atoms)
    print(f"\n  Ligand heavy atoms: {len(lig_coords)}")

    # Pocket-lining residues
    pocket_keys = set()
    for ch in holo_struct.get_chains():
        for res in ch.get_residues():
            if res.get_id()[0] != " ": continue
            for atom in res.get_atoms():
                if not atom.element or atom.element == "H": continue
                coord = np.array(atom.get_vector().get_array())
                if np.sqrt(((lig_coords - coord[None, :])**2).sum(axis=1)).min() <= 8.0:
                    pocket_keys.add((ch.get_id(), res.get_id()[1]))
                    break

    total = len(pocket_keys)
    h_n = sum(1 for k in pocket_keys if k in helix_res)
    s_n = sum(1 for k in pocket_keys if k in sheet_res)
    c_n = total - h_n - s_n
    hs  = (h_n + s_n) / total if total else 0.0
    dssp_pass = hs >= 0.60
    print(f"  Pocket residues: {total}  |  helix={h_n} sheet={s_n} coil={c_n}  →  {hs:.1%}")
    print(f"  DSSP gate: {'PASS' if dssp_pass else 'FAIL'}")

    if not dssp_pass:
        return {"label": label, "apo": apo_id, "holo": holo_id,
                "cofactor_diff": cofactor_diff, "hs_frac": hs, "dssp_pass": False,
                "pocket_ca_rmsd": None, "rmsd_pass": None}

    # ── Pocket RMSD gate ──
    apo_struct = parser.get_structure(apo_id, apo_path)[0]

    def get_ca(model, chain_id):
        try:
            ch = model[chain_id]
        except KeyError:
            ch = list(model.get_chains())[0]
        atoms = {}
        for res in ch.get_residues():
            if res.get_id()[0] == " " and "CA" in res:
                atoms[res.get_id()[1]] = res["CA"]
        return atoms

    apo_ca  = get_ca(apo_struct,  apo_chain)
    holo_ca = get_ca(holo_struct, holo_chain)
    shared_all = sorted(set(apo_ca.keys()) & set(holo_ca.keys()))

    # Global alignment
    sup = Superimposer()
    sup.set_atoms([holo_ca[r] for r in shared_all], [apo_ca[r] for r in shared_all])
    sup.apply(apo_struct.get_atoms())
    print(f"\n  Global CA RMSD (all shared): {sup.rms:.3f} Å  ({len(shared_all)} residues)")

    # Pocket CA RMSD (post-alignment)
    pocket_rids = {rid for _, rid in pocket_keys}
    common_pocket = sorted(set(apo_ca.keys()) & set(holo_ca.keys()) & pocket_rids)
    if not common_pocket:
        print("  No pocket CA pairs found after alignment")
        return None

    apo_arr  = np.array([apo_ca[r].get_vector().get_array()  for r in common_pocket])
    holo_arr = np.array([holo_ca[r].get_vector().get_array() for r in common_pocket])
    ca_rmsd  = np.sqrt(((apo_arr - holo_arr)**2).sum(axis=1).mean())
    rmsd_pass = 1.0 <= ca_rmsd <= 2.5
    print(f"  Pocket CA RMSD (n={len(common_pocket)}): {ca_rmsd:.3f} Å")
    print(f"  Pocket RMSD gate (1.0–2.5 Å): {'PASS' if rmsd_pass else 'FAIL'}")

    return {"label": label, "apo": apo_id, "holo": holo_id,
            "cofactor_diff": cofactor_diff, "hs_frac": hs, "dssp_pass": dssp_pass,
            "pocket_ca_rmsd": ca_rmsd, "rmsd_pass": rmsd_pass,
            "n_pocket": total, "global_rmsd": sup.rms}


# ── LdG6PD ──────────────────────────────────────────────────────────────────
# 4APN (JV0, 423 Da, pyrazole non-substrate inhibitor) vs.
# 4ADW (apo: has NDP/NADP only — confirms cofactor present in both)
r1 = run_gates(
    "LdG6PD — Leishmania donovani G6PD",
    apo_id="4ADW", apo_chain="A",
    holo_id="4APN", holo_chain="A",
    lig_codes=["JV0"],
    cofactor_note="NADP (NDP) should be present in both; JV0 is non-substrate pyrazole inhibitor"
)

# ── LdNMT ───────────────────────────────────────────────────────────────────
# 2WUU (no drug, has NHW = CoA fragment) vs 6QD9 (HWT 332 Da, thienopyrimidine)
r2 = run_gates(
    "LdNMT — Leishmania donovani NMT",
    apo_id="2WUU", apo_chain="A",
    holo_id="6QD9", holo_chain="A",
    lig_codes=["HWT"],
    cofactor_note="NHW (CoA fragment) may differ; need to check cofactor parity"
)

# ── MtAlaDH ─────────────────────────────────────────────────────────────────
# AlaDH has two binding sites: NAD(H) site + pyruvate/amino acid site.
# Drug-like inhibitors (benzofuran/indole series) occupy the pyruvate site.
# Apo: 2VHX (NAD + PYR) — PYR is substrate and leaves the site open to drugs
# Actually better apo: structure with NAD only, no substrate, no drug.
# Use 1MOP (no cofactor, no drug) vs 3ISJ (A8D, 268 Da):
# But 1MOP has SO4/GOL — need to check if NAD absent matters (cofactor confound).
# Safer pair: 2VHX (NAD only, no drug) vs 3ISJ (NAD + A8D).
# But 2VHX has PYR (substrate) — check if PYR occupies drug binding site or different site.
# Checking 2VHX HETATM: NAD, 2OP (2-oxo-propanoate = pyruvate), MG
# 3ISJ has: A8D (268 Da) — what else?
r3 = run_gates(
    "MtAlaDH — Mtb alanine dehydrogenase",
    apo_id="2VHX", apo_chain="A",
    holo_id="3ISJ", holo_chain="A",
    lig_codes=["A8D"],
    cofactor_note="2VHX has NAD+pyruvate (2OP); 3ISJ has A8D — need to verify NAD parity"
)

print("\n\n" + "="*70)
print("SUMMARY — New candidates")
print("="*70)
for r in [r1, r2, r3]:
    if r:
        pair = f"{r['apo']}/{r['holo']}"
        rmsd_str = f"{r['pocket_ca_rmsd']:.3f}" if r["pocket_ca_rmsd"] else "N/A"
        cofac = "FAIL" if r["cofactor_diff"] else "PASS"
        print(f"{pair:<12}: cofactor={cofac}, DSSP={r['hs_frac']:.1%} {'✓' if r['dssp_pass'] else '✗'}, "
              f"CA RMSD={rmsd_str} {'✓' if r['rmsd_pass'] else ('✗' if r['rmsd_pass'] is not None else 'N/A')}")

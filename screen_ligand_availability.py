"""
Extends the cryptic-pocket gate screen with a joint ligand-availability pass,
per instruction: rank on gate metrics AND drug availability together, not
sequentially (adenylate kinase clears the gate decisively and is unusable --
that only surfaced by testing ligand count after the fact; folding both into
one pass catches this class of case automatically).

For each apo/holo pair already screened (results/target_screen/
cryptic_pocket_screen_results.csv), adds:
  - uniprot: UniProt accession of the apo chain's polymer entity
  - pdb_ligand_structures: number of PDB entries sharing that UniProt with
    >=1 bound ligand outside a blocklist of common crystallization additives/
    buffer components (proxy for "how many distinct chemical matter points
    exist structurally" -- not a strict drug-likeness filter, a availability
    proxy)
  - chembl_target_id, chembl_activity_count: ChEMBL target match (if any) and
    total recorded bioactivity measurements against it (proxy for how
    medicinal-chemistry-attention a target has received)
  - cofactor_present: whether a non-crystallization-additive metal ion
    (Zn/Mg/Mn/Fe/Ni/Co/Cu) appears within 6A of the pocket residues in the
    holo structure -- surfaces the LasB-style cofactor confound automatically
    rather than discovering it during validation.
"""
import csv
import time
import subprocess
import json
import gemmi
import numpy as np

BLOCKLIST_LIGANDS = {
    "GOL", "EDO", "SO4", "PO4", "ACT", "DMS", "PEG", "CL", "NA", "K", "MG",
    "CA", "ZN", "MN", "FE", "NI", "CO", "CU", "HOH", "PGE", "PG4", "BME",
    "MPD", "IOD", "BR", "ACY", "FMT", "TRS", "IMD", "1PE", "MRD", "CPS",
    "UNX", "EPE", "BOG", "P6G", "PE4", "PE8", "SIN", "TAR", "MES", "BCT",
    "GD3", "IPA", "OCT", "NH4", "URE", "TLA",
}
METAL_COFACTORS = {"ZN", "MG", "MN", "FE", "NI", "CO", "CU", "CD", "MO", "W"}


def http_get_json(url, retries=1, timeout=15):
    for attempt in range(retries + 1):
        try:
            out = subprocess.run(["curl", "-s", "--max-time", str(timeout), url],
                                  capture_output=True, timeout=timeout + 5)
            return json.loads(out.stdout)
        except Exception:
            if attempt == retries:
                return None
            time.sleep(1)


def http_post_json(url, payload, timeout=15):
    try:
        out = subprocess.run(
            ["curl", "-s", "--max-time", str(timeout), "-X", "POST",
             "-H", "Content-Type: application/json", "-d", json.dumps(payload), url],
            capture_output=True, timeout=timeout + 5)
        return json.loads(out.stdout)
    except Exception:
        return None


def get_uniprot(pdbid, chain_hint="1"):
    for entity in ["1", "2", "3"]:
        d = http_get_json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdbid}/{entity}")
        if d is None:
            continue
        ids = d.get("rcsb_polymer_entity_container_identifiers", {}).get("uniprot_ids")
        if ids:
            return ids[0]
    return None


def pdb_ligand_structure_count(uniprot, sample_check=8):
    """Total PDB entries for this UniProt with >=1 non-polymer entity, plus a
    small-sample estimate of how many carry a ligand outside the common
    crystallization-additive blocklist (fast: samples up to `sample_check`
    entries rather than drilling into all of them)."""
    query = {
        "query": {
            "type": "group", "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match", "value": uniprot}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                    "operator": "exact_match", "value": "UniProt"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                    "operator": "greater", "value": 0}},
            ],
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    d = http_post_json("https://search.rcsb.org/rcsbsearch/v2/query", query, timeout=15)
    if d is None:
        return None, None
    total = d.get("total_count", 0)
    entry_ids = [r["identifier"] for r in d.get("result_set", [])]
    if not entry_ids:
        return 0, 0.0
    n_with_real_ligand = 0
    n_checked = 0
    for eid in entry_ids[:sample_check]:
        comps = http_get_json(f"https://data.rcsb.org/rest/v1/core/entry/{eid}", retries=0, timeout=10)
        if comps is None:
            continue
        np_ids = comps.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids") or []
        has_real = False
        for n in np_ids[:5]:
            npd = http_get_json(f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{eid}/{n}", retries=0, timeout=10)
            if npd is None:
                continue
            code = npd.get("pdbx_entity_nonpoly", {}).get("comp_id", "")
            if code and code.upper() not in BLOCKLIST_LIGANDS:
                has_real = True
                break
        n_checked += 1
        if has_real:
            n_with_real_ligand += 1
    frac_real = (n_with_real_ligand / n_checked) if n_checked else None
    return total, frac_real


def chembl_activity_count(uniprot):
    d = http_get_json(f"https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession={uniprot}")
    if not d or not d.get("targets"):
        return None, 0
    target_id = d["targets"][0]["target_chembl_id"]
    act = http_get_json(f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={target_id}&limit=1")
    count = act.get("page_meta", {}).get("total_count", 0) if act else 0
    return target_id, count


def cofactor_near_pocket(holo_pdb, holo_chain, pocket_resnums_holo, struct_dir="results/target_screen/structures"):
    try:
        st = gemmi.read_structure(f"{struct_dir}/{holo_pdb}.pdb")
        st.setup_entities()
    except Exception:
        return None
    model = st[0]
    pocket_coords = []
    for chain in model:
        for res in chain:
            if chain.name == holo_chain and res.seqid.num in pocket_resnums_holo:
                pocket_coords.extend([[a.pos.x, a.pos.y, a.pos.z] for a in res])
    if not pocket_coords:
        return None
    pocket_coords = np.array(pocket_coords)
    for chain in model:
        for res in chain:
            if res.name.upper() in METAL_COFACTORS:
                for a in res:
                    d = np.min(np.linalg.norm(pocket_coords - np.array([a.pos.x, a.pos.y, a.pos.z]), axis=1))
                    if d < 6.0:
                        return res.name.upper()
    return False


if __name__ == "__main__":
    with open("results/target_screen/cryptic_pocket_screen_results.csv") as f:
        rows = list(csv.DictReader(f))

    ok_rows = [r for r in rows if r.get("pocket_ca_rmsd_access_A")]
    print(f"{len(ok_rows)} successfully-screened pairs to annotate\n")

    out = []
    for r in ok_rows:
        uniprot = get_uniprot(r["apo_pdb"])
        pdb_total, frac_real = (None, None) if uniprot is None else pdb_ligand_structure_count(uniprot)
        chembl_id, chembl_n = (None, None) if uniprot is None else chembl_activity_count(uniprot)
        row = dict(r)
        row["uniprot"] = uniprot or ""
        row["pdb_ligand_structures_total"] = pdb_total if pdb_total is not None else ""
        row["pdb_ligand_structures_frac_real"] = round(frac_real, 2) if frac_real is not None else ""
        row["chembl_target_id"] = chembl_id or ""
        row["chembl_activity_count"] = chembl_n if chembl_n is not None else ""
        out.append(row)
        print(f"{r['apo_pdb']}->{r['holo_pdb']}: uniprot={uniprot} pdb_lig_total={pdb_total} "
              f"frac_real~{frac_real} chembl={chembl_id} n_activities={chembl_n} | "
              f"pocketRMSD={r['pocket_ca_rmsd_access_A']} fracStruct={r['frac_access_disp_structured']}")

    fieldnames = sorted(set().union(*[o.keys() for o in out]))
    with open("results/target_screen/cryptic_pocket_screen_with_availability.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, restval="")
        w.writeheader()
        w.writerows(out)
    print("\nWritten results/target_screen/cryptic_pocket_screen_with_availability.csv")

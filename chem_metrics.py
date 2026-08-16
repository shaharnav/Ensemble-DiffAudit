"""
Chemistry and pose-quality annotations for ensemble_auditor.py results (Phase 3).

Raw Vina score scales with molecule size, so "better affinity" isn't interpretable
without ligand efficiency, drug-likeness, and structural-alert screening alongside it.
This module adds those columns, plus PoseBusters geometry/clash validation on each
candidate's actual winning docked pose (not just its SMILES).

Usage:
    from chem_metrics import annotate
    df = annotate(pd.read_csv("results.csv"))
"""
import json
import os
import sys

import pandas as pd
from rdkit import Chem, RDConfig
from rdkit.Chem import Descriptors, QED
from meeko import PDBQTMolecule, RDKitMolCreate
from posebusters import PoseBusters

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

ENSEMBLE_WORK_DIR = os.path.join("results", "ensemble_audit")
ALIGNED_RECEPTORS_DIR = os.path.join("results", "payload_unpacked", "aligned_receptors")

STRUCTURAL_ALERT_SMARTS = {
    "epoxide": "C1OC1",
    "aziridine": "C1NC1",
    "thiirane": "C1SC1",
    "allene": "C=C=C",
    "peroxide": "OO",
    "acyl_halide": "C(=O)[F,Cl,Br,I]",
    "michael_acceptor": "C=CC(=O)",
    "alkyl_halide": "[CX4][Cl,Br,I]",
    "free_thiol": "[#6][SX2H]",
    "phosphonic_acid": "P(=O)(O)O",
    "aldehyde": "[CX3H1](=O)[#6]",
}
_ALERT_PATTERNS = {
    name: Chem.MolFromSmarts(smarts) for name, smarts in STRUCTURAL_ALERT_SMARTS.items()
}


def structural_alerts(mol: Chem.Mol) -> list[str]:
    """Names of every structural-alert SMARTS pattern that matches *mol*."""
    return [name for name, patt in _ALERT_PATTERNS.items() if mol.HasSubstructMatch(patt)]


def _find_docked_pose_pdbqt(smiles: str, winning_structure: str, n_candidates: int) -> str | None:
    """
    Locate the out.pdbqt for the docking job that actually produced *winning_structure*
    for this exact SMILES, by reading back docking_engine's per-job params sidecars
    (job_name encodes an index + receptor stem, not the SMILES). Only indices within
    the current candidate count are trusted -- results/ensemble_audit/ accumulates
    sidecars from every ensemble_auditor.py run ever performed in this repo.
    """
    stem = os.path.splitext(winning_structure)[0]
    for i in range(n_candidates):
        sidecar = os.path.join(ENSEMBLE_WORK_DIR, f"ens_s{i:04d}_{stem}_params.json")
        if not os.path.exists(sidecar):
            continue
        with open(sidecar) as f:
            params = json.load(f)
        if params.get("smiles") == smiles:
            out_pdbqt = os.path.join(ENSEMBLE_WORK_DIR, f"ens_s{i:04d}_{stem}_out.pdbqt")
            if os.path.exists(out_pdbqt):
                return out_pdbqt
    return None


def _docked_pose_mol(pdbqt_path: str) -> Chem.Mol:
    """Top-scoring (first) pose, correctly bond-ordered via Meeko."""
    pdbqt_mol = PDBQTMolecule.from_file(pdbqt_path, skip_typing=True)
    return RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]


def _run_posebusters(pb: PoseBusters, docked_mol: Chem.Mol, receptor_path: str) -> tuple[bool, str]:
    result = pb.bust(mol_pred=docked_mol, mol_cond=receptor_path, full_report=False)
    row = result.iloc[0]
    failed = [col for col, passed in row.items() if not bool(passed)]
    return len(failed) == 0, ";".join(failed)


def annotate(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of *df* with chemistry, drug-likeness, and pose-quality columns added."""
    df = df.copy()
    n = len(df)
    pb = PoseBusters(config="dock")
    smiles_col = "smiles" if "smiles" in df.columns else "SMILES"

    heavy_atoms, mw, logp, hbd, hba, rot_bonds, tpsa = [], [], [], [], [], [], []
    lipinski_violations, qed_vals, sa_vals = [], [], []
    alerts_joined, n_alerts = [], []
    pb_pass, pb_failures = [], []

    for _, row in df.iterrows():
        smiles = row[smiles_col]
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            heavy_atoms.append(None); mw.append(None); logp.append(None)
            hbd.append(None); hba.append(None); rot_bonds.append(None); tpsa.append(None)
            lipinski_violations.append(None); qed_vals.append(None); sa_vals.append(None)
            alerts_joined.append(""); n_alerts.append(None)
            pb_pass.append(None); pb_failures.append("invalid SMILES")
            continue

        heavy_atoms.append(mol.GetNumHeavyAtoms())
        mol_wt = Descriptors.MolWt(mol)
        mol_logp = Descriptors.MolLogP(mol)
        n_hbd = Descriptors.NumHDonors(mol)
        n_hba = Descriptors.NumHAcceptors(mol)
        mw.append(mol_wt)
        logp.append(mol_logp)
        hbd.append(n_hbd)
        hba.append(n_hba)
        rot_bonds.append(Descriptors.NumRotatableBonds(mol))
        tpsa.append(Descriptors.TPSA(mol))
        lipinski_violations.append(
            sum([mol_wt > 500, mol_logp > 5, n_hbd > 5, n_hba > 10])
        )
        qed_vals.append(QED.qed(mol))
        sa_vals.append(sascorer.calculateScore(mol))

        alerts = structural_alerts(mol)
        alerts_joined.append(";".join(alerts))
        n_alerts.append(len(alerts))

        winning_structure = row.get("overall_best_structure") or row.get("winning_conformation")
        pdbqt_path = (
            _find_docked_pose_pdbqt(smiles, winning_structure, n) if winning_structure else None
        )
        receptor_path = (
            os.path.join(ALIGNED_RECEPTORS_DIR, winning_structure) if winning_structure else None
        )
        if pdbqt_path and receptor_path and os.path.exists(receptor_path):
            try:
                docked_mol = _docked_pose_mol(pdbqt_path)
                passed, failures = _run_posebusters(pb, docked_mol, receptor_path)
                pb_pass.append(passed)
                pb_failures.append(failures)
            except Exception as exc:
                pb_pass.append(None)
                pb_failures.append(f"error: {exc}")
        else:
            pb_pass.append(None)
            pb_failures.append("no docked pose found")

    df["heavy_atoms"] = heavy_atoms
    affinity_col = "overall_best_affinity" if "overall_best_affinity" in df.columns else "true_affinity"
    df["ligand_efficiency"] = [
        (-affinity / ha) if (affinity is not None and pd.notna(affinity) and ha) else None
        for affinity, ha in zip(df[affinity_col], heavy_atoms)
    ]
    df["mw"] = mw
    df["logp"] = logp
    df["hbd"] = hbd
    df["hba"] = hba
    df["rot_bonds"] = rot_bonds
    df["tpsa"] = tpsa
    df["lipinski_violations"] = lipinski_violations
    df["qed"] = qed_vals
    df["sa_score"] = sa_vals
    df["structural_alerts"] = alerts_joined
    df["n_structural_alerts"] = n_alerts
    df["posebusters_pass"] = pb_pass
    df["posebusters_failures"] = pb_failures

    return df


if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    if not os.path.exists("results.csv"):
        logger.error("results.csv not found. Run ensemble_auditor.py first.")
        sys.exit(1)

    df = pd.read_csv("results.csv")
    annotated = annotate(df)
    annotated.to_csv("results_annotated.csv", index=False)
    logger.info(f"Annotated {len(annotated)} candidates -> results_annotated.csv")

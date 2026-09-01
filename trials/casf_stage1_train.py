"""
Stage 1 -- delta model. No PDBbind refined set (registration wall), so per the
Stage -1 fallback: leave-one-protein-family-out CV WITHIN the 238-complex
CASF-2016-derived set itself (UniProt-based grouping, 71 groups). This is a
materially smaller and differently-structured training regime than the
published papers (trained on ~1300-5000 PDBbind complexes, evaluated on a
fully separate CASF-2016 holdout) -- flagged as small-N throughout, per
instructions.

REPRODUCTION GATE (stated before running, per instructions):
Published ΔvinaRF20 (Wang & Zhang 2017) CASF-2016 scoring power: Pearson R ~ 0.803.
Published ΔvinaXGB (Lu et al. 2019) CASF-2016 scoring power: Pearson R ~ 0.82-0.86.
Vina baseline (this reproduction, Stage 0): Pearson R = 0.554 (n=244).
Gate: Pearson R within ~0.1 of 0.80 (i.e. >= 0.70) is a defensible reproduction.
Caveat stated up front: with an order of magnitude less training data (~230 vs
1300-5000 complexes) and a from-scratch approximate feature set (geometric
Vina terms, geometric H-bonds, FreeSASA not MSMS), landing at 0.70 is optimistic
-- a materially lower number would not by itself prove broken features, given
these confounds. Reported honestly either way.
"""
import pandas as pd, numpy as np
from scipy import stats
from sklearn.model_selection import LeaveOneGroupOut
from xgboost import XGBRegressor

df = pd.read_csv("/tmp/casf_modeling_ready.csv")
feat_cols = [c for c in df.columns if c not in ("pdbid","pKd","vina","vinardo","uniprot_ids","contact_C_0_2","contact_S_0_2")]
print(f"n = {len(df)}, features = {len(feat_cols)}, groups = {df['uniprot_ids'].nunique()}")

X = df[feat_cols].values
y_pkd = df["pKd"].values
vina_score = df["vina"].values
vinardo_score = df["vinardo"].values
y_residual = y_pkd - vina_score
groups = df["uniprot_ids"].values

logo = LeaveOneGroupOut()
oof_pred_residual = np.zeros(len(df))
n_folds = 0
for train_idx, test_idx in logo.split(X, y_residual, groups):
    model = XGBRegressor(n_estimators=200, max_depth=4, learning_rate=0.05,
                          subsample=0.8, colsample_bytree=0.8, random_state=42, n_jobs=-1)
    model.fit(X[train_idx], y_residual[train_idx])
    oof_pred_residual[test_idx] = model.predict(X[test_idx])
    n_folds += 1

oof_pred_pkd = oof_pred_residual + vina_score

print(f"\nCompleted {n_folds}-fold leave-one-family-out CV")

def scoring_power(pred, true):
    r, p = stats.pearsonr(pred, true)
    return r, p

def ranking_power(pred, true, groups):
    rhos = []
    for g in set(groups):
        idx = groups == g
        if idx.sum() >= 3:
            rho, _ = stats.spearmanr(pred[idx], true[idx])
            if not np.isnan(rho):
                rhos.append(rho)
    return np.mean(rhos), len(rhos)

print("\n=== STAGE 1 RESULTS ===\n")
r_delta, p_delta = scoring_power(oof_pred_pkd, y_pkd)
rho_delta, n_groups_delta = ranking_power(oof_pred_pkd, y_pkd, groups)
print(f"Delta model (LOGO-CV, n={len(df)}):")
print(f"  Scoring power (Pearson R) = {r_delta:.3f} (p={p_delta:.2e})")
print(f"  Ranking power (mean within-family Spearman, {n_groups_delta} families n>=3) = {rho_delta:.3f}")

r_vina, p_vina = scoring_power(vina_score, y_pkd)  # note: score is negative-good, so sign flips
r_vina = -r_vina
rho_vina, _ = ranking_power(-vina_score, y_pkd, groups)
print(f"\nVina baseline (same n={len(df)} subset):")
print(f"  Scoring power (Pearson R) = {r_vina:.3f}")
print(f"  Ranking power = {rho_vina:.3f}")

r_vinardo = -stats.pearsonr(vinardo_score, y_pkd)[0]
rho_vinardo, _ = ranking_power(-vinardo_score, y_pkd, groups)
print(f"\nVinardo baseline:")
print(f"  Scoring power (Pearson R) = {r_vinardo:.3f}")
print(f"  Ranking power = {rho_vinardo:.3f}")

print(f"\n=== REPRODUCTION GATE ===")
print(f"Published ΔvinaRF20 CASF-2016 scoring power: ~0.803")
print(f"Published ΔvinaXGB CASF-2016 scoring power: ~0.82-0.86")
print(f"This reproduction (delta model, LOGO-CV): {r_delta:.3f}")
margin = 0.80 - r_delta
print(f"Gap from 0.80 target: {margin:.3f}")
gate_pass = r_delta >= 0.70
print(f"GATE (>=0.70, i.e. within ~0.1 of 0.80): {'PASS' if gate_pass else 'FAIL'}")

results = pd.DataFrame({"pdbid": df["pdbid"], "pKd": y_pkd, "vina": vina_score,
                         "vinardo": vinardo_score, "delta_pred_pkd": oof_pred_pkd,
                         "uniprot_group": groups})
results.to_csv("results/casf2016/stage1_oof_predictions.csv", index=False)
print("\nWritten results/casf2016/stage1_oof_predictions.csv")

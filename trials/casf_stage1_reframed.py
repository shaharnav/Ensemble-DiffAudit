"""
Reframed Stage 1 + Stage 2, per direction after the delta-model gate failure.

The >=0.70 gate (calibrated against ΔvinaRF20 trained on ~4000 PDBbind complexes)
is uninformative at n~230 -- dropped. Replaced with an honest bar: beat the Vina
baseline on the same held-out split, with bootstrap CIs so "beat" is judged
against actual uncertainty, not a point estimate.

Four models compared on the identical LOGO-CV split:
  1. Vina baseline (no model)
  2. Delta/residual model (all 44 features, target = pKd - vina_score) -- the
     formulation that failed
  3. Direct regression, B-F features + vina_score as a plain input
  4. Direct regression, B-F features only (no vina_score) -- isolates whether
     the win is independent signal or "trust Vina, adjust slightly"

Stage 2 decorrelation then runs on model 3 (the strongest validated model),
not on the failed delta formulation.
"""
import pandas as pd, numpy as np
from scipy import stats
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor

df = pd.read_csv("/tmp/casf_modeling_ready.csv")
vina_cols = ['vina_gauss1','vina_gauss2','vina_repulsion','vina_hydrophobic','vina_hbond']
feat_cols_all = [c for c in df.columns if c not in ("pdbid","pKd","vina","vinardo","uniprot_ids","contact_C_0_2","contact_S_0_2")]
feat_cols_bf = [c for c in feat_cols_all if c not in vina_cols]

y_pkd = df["pKd"].values
hac = df["hac"].values
vina_score = df["vina"].values
groups = df["uniprot_ids"].values
logo = LeaveOneGroupOut()

def ridge_oof(X, y, alpha=50.0):
    oof = np.zeros(len(y))
    for tr, te in logo.split(X, y, groups):
        scaler = StandardScaler().fit(X[tr])
        Xtr, Xte = scaler.transform(X[tr]), scaler.transform(X[te])
        m = Ridge(alpha=alpha)
        m.fit(Xtr, y[tr])
        oof[te] = m.predict(Xte)
    return oof

def bootstrap_ci_pearson(x, y, n=10000, seed=0):
    rng = np.random.default_rng(seed)
    x, y = np.array(x), np.array(y)
    N = len(x)
    vals = []
    for _ in range(n):
        idx = rng.integers(0, N, N)
        if np.std(x[idx]) == 0 or np.std(y[idx]) == 0:
            continue
        vals.append(stats.pearsonr(x[idx], y[idx])[0])
    return np.percentile(vals, 2.5), np.percentile(vals, 97.5)

results = {}

# 1. Vina baseline
results["Vina baseline"] = -vina_score  # sign-flip: higher predicted = tighter binding

# 2. Delta/residual model (all 44 features) -- reuse from prior run
X_all = df[feat_cols_all].values
y_residual = y_pkd - vina_score
oof_resid = np.zeros(len(df))
for tr, te in logo.split(X_all, y_residual, groups):
    m = XGBRegressor(n_estimators=200, max_depth=4, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, random_state=42, n_jobs=-1)
    m.fit(X_all[tr], y_residual[tr])
    oof_resid[te] = m.predict(X_all[te])
results["Delta model (residual target, all features)"] = oof_resid + vina_score

# 3. Direct regression, B-F + vina_score
X_bf_plus = np.column_stack([df[feat_cols_bf].values, vina_score])
results["Direct regression (B-F + vina_score)"] = ridge_oof(X_bf_plus, y_pkd, alpha=50.0)

# 4. Direct regression, B-F only (no vina_score)
X_bf = df[feat_cols_bf].values
results["Direct regression (B-F only, no vina_score)"] = ridge_oof(X_bf, y_pkd, alpha=50.0)

print("=== FOUR-MODEL COMPARISON, identical LOGO-CV split (n=238) ===\n")
print(f"{'Model':45s} {'Pearson R':>10s} {'95% CI':>20s} {'r(score,HAC)':>14s}")
summary_rows = []
for name, pred in results.items():
    r, p = stats.pearsonr(pred, y_pkd)
    ci = bootstrap_ci_pearson(pred, y_pkd)
    r_hac, _ = stats.pearsonr(pred, hac)
    print(f"{name:45s} {r:10.3f} {'['+f'{ci[0]:.3f}'+', '+f'{ci[1]:.3f}'+']':>20s} {r_hac:14.3f}")
    summary_rows.append({"model": name, "pearson_r": r, "ci_low": ci[0], "ci_high": ci[1], "r_hac": r_hac})

pd.DataFrame(summary_rows).to_csv("results/casf2016/stage1_four_model_comparison.csv", index=False)

# CI overlap check: does model 3's CI exclude Vina's point estimate / overlap Vina's CI?
vina_r, _ = stats.pearsonr(results["Vina baseline"], y_pkd)
vina_ci = bootstrap_ci_pearson(results["Vina baseline"], y_pkd)
m3_r, _ = stats.pearsonr(results["Direct regression (B-F + vina_score)"], y_pkd)
m3_ci = bootstrap_ci_pearson(results["Direct regression (B-F + vina_score)"], y_pkd)
print(f"\nVina baseline: {vina_r:.3f} CI{vina_ci}")
print(f"Direct model (B-F+vina): {m3_r:.3f} CI{m3_ci}")
print(f"CIs overlap: {not (m3_ci[0] > vina_ci[1] or vina_ci[0] > m3_ci[1])}")

# ---------------- STAGE 2: decorrelation, on model 3 (best validated model) ----------------
print("\n\n=== STAGE 2: SIZE DECORRELATION (on direct regression, B-F + vina_score) ===\n")

# (i) residualize each feature against HAC (fit on train fold only), then Ridge on residuals
def residualize_against_hac(X, hac_vals, train_idx):
    X_resid = np.zeros_like(X)
    hac_train = hac_vals[train_idx].reshape(-1, 1)
    for j in range(X.shape[1]):
        col = X[:, j]
        # simple linear fit of col ~ hac on train, apply to all
        A = np.vstack([hac_vals[train_idx], np.ones(len(train_idx))]).T
        coef, resid_ss, rank, sv = np.linalg.lstsq(A, col[train_idx], rcond=None)
        pred_all = coef[0]*hac_vals + coef[1]
        X_resid[:, j] = col - pred_all
    return X_resid

oof_decorr_feat = np.zeros(len(df))
for tr, te in logo.split(X_bf_plus, y_pkd, groups):
    X_resid = residualize_against_hac(X_bf_plus, hac, tr)
    scaler = StandardScaler().fit(X_resid[tr])
    Xtr, Xte = scaler.transform(X_resid[tr]), scaler.transform(X_resid[te])
    m = Ridge(alpha=50.0)
    m.fit(Xtr, y_pkd[tr])
    oof_decorr_feat[te] = m.predict(Xte)

r_decorr_feat, _ = stats.pearsonr(oof_decorr_feat, y_pkd)
ci_decorr_feat = bootstrap_ci_pearson(oof_decorr_feat, y_pkd)
r_hac_decorr_feat, _ = stats.pearsonr(oof_decorr_feat, hac)
print(f"(i) Feature-residualized: R={r_decorr_feat:.3f} CI{ci_decorr_feat}  r(score,HAC)={r_hac_decorr_feat:.3f}")

# (ii) target = ligand efficiency (pKd / HAC)
y_le = y_pkd / hac
oof_le = ridge_oof(X_bf_plus, y_le, alpha=50.0)
oof_le_as_pkd = oof_le * hac  # convert back to pKd scale for comparison
r_le, _ = stats.pearsonr(oof_le_as_pkd, y_pkd)
ci_le = bootstrap_ci_pearson(oof_le_as_pkd, y_pkd)
r_hac_le, _ = stats.pearsonr(oof_le_as_pkd, hac)
r_hac_le_native, _ = stats.pearsonr(oof_le, y_le)  # native LE-space check
print(f"(ii) LE-target (converted back to pKd): R={r_le:.3f} CI{ci_le}  r(score,HAC)={r_hac_le:.3f}")

print(f"\nBaseline (direct model, model 3): R={m3_r:.3f} CI{m3_ci}  r(score,HAC)={stats.pearsonr(results['Direct regression (B-F + vina_score)'], hac)[0]:.3f}")

stage2_rows = [
    {"model": "Direct regression (baseline for Stage 2)", "pearson_r": m3_r, "ci_low": m3_ci[0], "ci_high": m3_ci[1], "r_hac": stats.pearsonr(results['Direct regression (B-F + vina_score)'], hac)[0]},
    {"model": "Decorrelated (i): features residualized vs HAC", "pearson_r": r_decorr_feat, "ci_low": ci_decorr_feat[0], "ci_high": ci_decorr_feat[1], "r_hac": r_hac_decorr_feat},
    {"model": "Decorrelated (ii): ligand-efficiency target", "pearson_r": r_le, "ci_low": ci_le[0], "ci_high": ci_le[1], "r_hac": r_hac_le},
]
pd.DataFrame(stage2_rows).to_csv("results/casf2016/stage2_decorrelation.csv", index=False)
print("\nWritten results/casf2016/stage1_four_model_comparison.csv and stage2_decorrelation.csv")

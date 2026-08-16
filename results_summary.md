# Results Summary (Phase 8)

N = 112 candidates generated and docked; 108 with complete crystal/ensemble data, 108 with complete rigid-control data.

## Crystal vs. ensemble (raw, not noise-corrected)

Paired t-test (ensemble_best_affinity vs crystal_affinity), N=108: t=-2.300, p=0.0234

Wilcoxon signed-rank, N=108: W=2246.5, p=0.0458

Mean delta_ensemble_vs_crystal = 0.074 kcal/mol (N=108)

## Noise-corrected comparison (rigid_max_over_seeds vs. ensemble_best_affinity) -- the test that matters

Paired t-test, N=108: t=-2.146, p=0.0341

Wilcoxon signed-rank, N=108: W=2139.5, p=0.0138

Mean noise_corrected_delta = -0.060 kcal/mol (median -0.060), 95% bootstrap CI [-0.114, -0.005] (10000 resamples), N=108

Cohen's d (paired) = -0.207

42/108 candidates (39%) beat the noise-matched rigid baseline.

The 95% CI on the mean excludes zero.

## Crystal beat every ensemble conformer

44/108 candidates (41%), N=108

## Spearman correlation: delta_ensemble_vs_crystal vs. winner_pocket_ca_rmsd

N=108: rho=-0.122, p=0.2101

## Distribution summaries

- Ligand efficiency (N=108): mean=0.282, median=0.271, sd=0.083, range=[0.152, 0.801]
- QED (N=112): mean=0.538, median=0.524, sd=0.150, range=[0.304, 0.901]
- SA score (N=112): mean=4.858, median=4.906, sd=0.752, range=[2.659, 5.985]
- MW (N=112): mean=334.6, median=348.1, sd=81.9, range=[79.1, 504.5]

## Notes

- Every statistic above reports N, the test used, and the p-value where applicable -- none should be described as significant without those.
- Raw affinity improvements are never reported without ligand efficiency alongside them (see the distribution summary above).

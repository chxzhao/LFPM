# Method and source provenance

Target manuscript: `lfpm_final_integration_20260929/LFPM_revised_20260929.tex`.
The associated simulation chapter is `outputs/simulation_20260929_english_JOE`.
Its numerical source is `simulation_ma2024_20260929`, using the preserved
`simulation_five_methods_20260927` data generator and comparison methods.

The code exposed by this R package is a factor-pricing implementation. In
particular, COM is the explicit OR--Bonferroni comparison and max-k is the
specified adaptation, not a claim of exact reproduction of every procedure in
the cited papers. Empirical and simulation comparator implementations must
remain identified separately.

## Sources cited by the manuscript

- Ma, Feng and Wang (2024), arXiv:2410.14308v1:
  https://arxiv.org/abs/2410.14308v1. Candidate grid and adaptive ordered sums.
- Feng, Lan, Liu and Ma (2022), Journal of Econometrics 229(1), 152--175:
  https://doi.org/10.1016/j.jeconom.2021.07.011. Maximum-test comparison and
  factor-pricing simulation context.
- Pesaran and Yamagata (2024): https://doi.org/10.1093/jjfinec/nbad002.
  Sum-type alpha test.
- Xia and Zhang (2024): https://doi.org/10.1080/07350015.2023.2217871.
  Precision-based adaptive benchmark motivating the documented max-k adaptation.
- Liu and Xie (2020), Cauchy combination test:
  https://doi.org/10.1080/01621459.2018.1554485.

## Software attribution

The simulation routines are local implementations using base R and standard R
packages. There is no required third-party R package implementing PY, MAX,
COM or max-k in this repository. In the source project,
`cvCovEst::thresholdingEst` was used only as an optional numerical cross-check.
It is not required to run this package.

The empirical backend is retained in Python to preserve the original NumPy
random streams. Its guide lists the actual Python packages rather than
inventing R-package equivalents. QUIC 1.2 is attributed to Cho-Jui Hsieh and
Matyas A. Sustik and was obtained through `skggm` 0.2.8:
https://github.com/skggm/skggm/tree/master/inverse_covariance/pyquic.
Its source and GPL-3 license accompany the optional backend. No precompiled
library or private R/Python runtime is distributed.

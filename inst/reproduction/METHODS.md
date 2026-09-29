# Method and source provenance

Target manuscript: *Testing factor pricing restrictions in high dimensions without knowing the sparsity of alphas*.
The associated simulation chapter is `outputs/simulation_20260929_english_JOE`.
Its numerical source is `simulation_ma2024_20260929`, using the preserved
`simulation_five_methods_20260927` generator and comparator routines.
The source paths identify the research project and are not runtime dependencies.

Version 0.1.1 distributes reusable R methods and simulation reproduction.
Case-study material is limited to [external data-source links](DATA_SOURCES.md);
case-study code, outputs and market observations are not distributed.

## Statistical definitions and references

The R tests are factor-pricing implementations. COM is the explicit
OR--Bonferroni combination of PY and MAX. The max-k comparator retains the
simulation's four-order grid, thresholded covariance estimate and original
ridge shift; it is a specified adaptation, not a claim of exact reproduction
of every procedure in the cited papers.

- Ma, Feng and Wang (2024), [arXiv:2410.14308v1](https://arxiv.org/abs/2410.14308v1):
  candidate grid and adaptive ordered sums.
- Feng, Lan, Liu and Ma (2022), Journal of Econometrics 229(1), 152--175,
  [DOI](https://doi.org/10.1016/j.jeconom.2021.07.011): maximum-test comparison
  and factor-pricing simulation context.
- Pesaran and Yamagata (2024),
  [DOI](https://doi.org/10.1093/jjfinec/nbad002): sum-type alpha test.
- Xia and Zhang (2024),
  [DOI](https://doi.org/10.1080/07350015.2023.2217871): precision-based adaptive
  benchmark motivating the documented max-k adaptation.
- Liu and Xie (2020),
  [Cauchy combination test](https://doi.org/10.1080/01621459.2018.1554485).

## Generic calibration choices

`adal_test()` supports two explicit component-calibration options on
user-supplied or synthetic return/factor matrices. `calibration="simulation"`
uses the manuscript simulation's empirical tail for order 5 and bootstrap-
moment normal tails for other dyadic orders, with machine-epsilon clipping.
`calibration="empirical"` uses empirical bootstrap tails for every order and
clipping at `1/(N+T)`. Both refit residual variances in each bootstrap draw.
The second option is a reusable statistical method, not a case-study script.
The simulation runner uses its original chapter-specific convention.

## Software attribution

The package uses base R and standard R packages. There is no required
third-party R package implementing PY, MAX, COM or max-k in this repository,
and these routines are not presented as official author-package calls.
In the source project, `cvCovEst::thresholdingEst` was used only as an optional
numerical cross-check. It is not needed to install or run LFPM.

The package is distributed under GPL-3 and has no compiled extension.

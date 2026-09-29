# Method and source provenance

Target manuscript: *Testing factor pricing restrictions in high dimensions without knowing the sparsity of alphas*.
The associated simulation chapter is `outputs/simulation_20260929_english_JOE`.
Its numerical source is `simulation_ma2024_20260929`, using the preserved
`simulation_five_methods_20260927` generator and comparator routines.
The source paths identify the research project and are not runtime dependencies.

Version 0.1.2 distributes reusable R methods, simulation reproduction, and
empirical reproduction code. Market observations and prepared panels are not
included. The [data-source page](DATA_SOURCES.md) documents acquisition and
input requirements; the [empirical guide](EMPIRICAL.md) documents computation.

## Statistical definitions and references

The tests are factor-pricing implementations. COM is the explicit
OR--Bonferroni combination of PY and MAX. The simulation max-k comparator
retains the four-order grid, thresholded covariance estimate and original
ridge shift. The empirical comparator uses QUIC precision estimation,
Gaussian multipliers and its specified candidate orders. These are documented
adaptations, not claims to reproduce every procedure in the cited papers.

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

## Calibration choices and random streams

`adal_test()` supports two component-calibration options on user-supplied
or synthetic matrices. `calibration="simulation"` uses the manuscript
simulation's empirical tail for order 5 and bootstrap-moment normal tails
for other dyadic orders, with machine-epsilon clipping.
`calibration="empirical"` uses empirical bootstrap tails for every order and
clipping at `1/(N+T)`. Both refit residual variances in each bootstrap draw.
The simulation runner preserves its original chapter-specific convention.

`lfpm_empirical_test()` and `reproduce_empirical()` use the Python backend to
retain the empirical implementation's NumPy PCG64 streams, seed derivation
and numerical batching. The generic R option has the same stated statistical
calibration convention but does not reproduce the NumPy draws from an
identical integer seed. Empirical and simulation comparator conventions also
remain distinct.

## Software attribution

The reusable R methods and simulations use base R and standard R packages.
There is no required third-party R package implementing PY, MAX, COM or max-k
in this repository, and these routines are not official author-package calls.
`cvCovEst::thresholdingEst` was used only as an optional numerical cross-check
in the source project; it is not needed to install or run LFPM.

The optional empirical backend uses the following actual Python packages:

| Software | Role |
|---|---|
| NumPy | Arrays, linear algebra and random draws |
| pandas | CSV inputs and result tables |
| SciPy | Statistical distributions and numerical routines |
| statsmodels | Residual diagnostics and multiple-testing adjustment |
| threadpoolctl | Numerical-library thread control |
| matplotlib | Optional empirical figures |

The precision estimator includes QUIC 1.2 source attributed to Cho-Jui Hsieh
and Matyas A. Sustik, obtained through
[skggm 0.2.8](https://github.com/skggm/skggm/tree/master/inverse_covariance/pyquic).
The source, GPL-3 license and upstream notices accompany the backend; no
precompiled solver binary is shipped. `skggm` itself is not a required Python
installation. The solver is compiled at runtime when requested, with the
platform requirements described in the empirical guide. Python source
provenance is recorded in `inst/python/SOURCE_MAP.json`.

## Verification boundary

Short numerical and synthetic integration tests are distinct from a complete
rerun of the simulations and empirical study. This package does not claim that
all published numbers have been newly verified. Current data downloads can
differ from the historical input vintages.

The final archived China max-k penalty-sensitivity endpoint lacks its original
seed and draw-budget metadata. The sensitivity interface performs a fresh
calculation using explicit recorded seeds and flags the missing historical
random-stream parity; it does not copy the paper's endpoint into its outputs.
This limitation is separate from the recovered seed streams for the main tests.

The package is distributed under GPL-3; included third-party source retains
its own upstream attribution and license notices.

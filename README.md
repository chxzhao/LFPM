# LFPM

An R package for *Testing factor pricing restrictions in high dimensions without knowing the sparsity of alphas*.
Version 0.1.1 contains reusable R methods and simulation reproduction.
Case-study material is limited to [external data-source links](inst/reproduction/DATA_SOURCES.md);
case-study analysis code and market data are not included.

Author: **Chenxi Zhao**.
Maintainer: Chenxi Zhao <chxzhao@mail.nankai.edu.cn>.
This is a research repository, not a CRAN release.

## Install

With R 4.2 or later, from this repository's root:

```sh
R CMD build .
R CMD INSTALL LFPM_0.1.1.tar.gz
```

The package uses R's standard packages and requires no compiled extensions.

## Test your own panel

Rows are observations (dates), columns of `Y` are asset excess returns, and
columns of `F` are factors. The functions add the intercept themselves.
This example uses synthetic data and shows two general component-calibration
options for the same reusable AdaL function.

```r
library(LFPM)
dat <- simulate_lfpm(model = 1, N = 40, T_size = 60,
                     distribution = "Normal", m = 0, seed = 42)
fit <- lfpm_fit(dat$Y, dat$F)
l_test(dat$Y, dat$F, k = 5, B = 99, seed = 43)
adal_test(dat$Y, dat$F, B = 99, calibration = "simulation", seed = 43)
adal_test(dat$Y, dat$F, B = 99, calibration = "empirical", seed = 43)
lfpm_test(dat$Y, dat$F, B = 99, seed = 43)
```

`B = 99` is a quick example, not the paper's simulation budget of 999 bootstrap
draws. The `"empirical"` option means empirical bootstrap-tail probabilities;
it is a generic statistical calibration option, not a case-study workflow.
See the help pages for return values, input restrictions and seed handling.

## Reproduce the simulations

```r
simulation_design()
result <- reproduce_simulation(output_dir = "reproduction-output/smoke",
                               preset = "smoke")
plot_simulation(result, file = "reproduction-output/smoke/power.pdf")
```

The `paper` preset uses all 96 size scenarios and 3200 power scenarios,
including every m from 1 through 200. Its budgets are 1000 independent
evaluation repetitions, 1999 independent null-calibration repetitions and
999 bootstrap draws per repetition. This is a large computation, with
checkpoint support. Read `?reproduce_simulation` and
[the simulation guide](inst/reproduction/SIMULATION.md) before starting it.

```r
result <- reproduce_simulation(output_dir = "reproduction-output/paper",
                               preset = "paper")
```

Bundled **simulation summary tables** are cached reference outputs.
`reproduction_reference()` reads them without rerunning the experiments.
The tables and plots distinguish cached references from newly computed runs.

## Calibration options

The reusable AdaL function makes its calibration convention explicit:

| Convention | `calibration = "simulation"` | `calibration = "empirical"` |
|---|---|---|
| Candidate orders | 5 plus the complete Ma dyadic grid | Same grid |
| Order 5 p-value | Add-one empirical bootstrap right tail | Same |
| Other component p-values | Normal tail using bootstrap mean and sample SD | Add-one empirical bootstrap right tail |
| Cauchy endpoint clipping | Machine epsilon | `1 / (N + T)` |
| Bootstrap studentization | OLS residual variance refitted each draw | Same |

The simulation runner preserves the manuscript's simulation convention.
The later all-bootstrap-tail simulation revision is not substituted for those
results. Independent-null calibration in the simulation is a known-DGP
diagnostic comparison. Finite-sample size and power experiments do not
establish an exact fixed-level theorem for the complete growing candidate grid.

## Other methods and software

| Method | Implementation in this package |
|---|---|
| L and AdaL | The manuscript's factor-model implementation, with explicit calibration choices |
| PY | Local implementation of the documented Pesaran--Yamagata correction |
| MAX | Local maximum-test implementation with the simulation's finite-sample scaling |
| COM | OR--Bonferroni combination of PY and MAX |
| max-k | The simulation's precision-transformed benchmark, inspired by Xia--Zhang, with its four-order grid and thresholded covariance estimate |

These are **not calls to official author R packages**. Covariance thresholding
uses base R. `cvCovEst::thresholdingEst` was an optional numerical validation
reference in the original project and is not a runtime dependency.
The benchmark is the documented adaptation, not a claim to reproduce every
published max-k procedure. See [method attribution](inst/reproduction/METHODS.md)
for the statistical references.

## Case-study data sources

The [data-source page](inst/reproduction/DATA_SOURCES.md) links to the sources
used for the paper's case studies. This release does not include case-study
analysis scripts, case-study reproduction/reporting functions, downloaded
files, cleaned market panels or credentials. Users obtain external data under
the providers' access conditions.

## Checks and reproducibility limits

```sh
R CMD build .
R CMD check --no-manual LFPM_0.1.1.tar.gz
```

A successful short smoke run does not mean the full simulation study has been
rerun. Floating-point differences across R/BLAS platforms can affect tail
comparisons even with the same random draws. Reproduction manifests validate
settings, budgets and engine compatibility before checkpoints are reused.

License: GPL-3.

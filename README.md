# LFPM

An R package for *Adaptive L-tests for alpha in high-dimensional linear factor
pricing models*, using the integrated manuscript `LFPM_revised_20260929`.

Author: **Chenxi Zhao**.
Maintainer: Chenxi Zhao <chxzhao@mail.nankai.edu.cn>.
This is a research repository, not a CRAN release.

## Install

With R 4.2 or later, from this repository's root:

```sh
R CMD build .
R CMD INSTALL LFPM_0.1.0.tar.gz
```

The core methods and simulations use R's standard packages. The optional
empirical workflow additionally requires Python scientific packages and,
for its QUIC precision estimator, a C++ compiler with BLAS/LAPACK.
These dependencies are documented in `inst/reproduction/EMPIRICAL.md`.

## Test your own panel

Rows are observations (dates), columns of `Y` are asset excess returns, and
columns of `F` are factors. The functions add the intercept themselves.

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

`B = 99` is a quick example, not the paper's bootstrap budget. The reference
simulation uses 999 bootstrap draws. See the help pages for full return values,
input restrictions and the handling of random seeds.

## Reproduce the paper

```r
simulation_design()
result <- reproduce_simulation(output_dir = "reproduction-output/smoke",
                               preset = "smoke")
```

The `paper` preset uses the complete design and original budgets; it is a large
computation. The runner supports checkpointed reproduction. Run the small
example first and read `?reproduce_simulation` and
`inst/reproduction/SIMULATION.md` before starting it.

```r
result <- reproduce_simulation(output_dir = "reproduction-output/paper",
                               preset = "paper")
```

Bundled **simulation summary tables** are reference outputs, not a substitute
for rerunning the experiments. `reproduction_reference()` accesses them.

For empirical reproduction, obtain the documented external datasets and use
`reproduce_empirical()`; `lfpm_empirical_test()` exposes the empirical test
definitions on user-supplied matrices. See `inst/reproduction/EMPIRICAL.md` and
`inst/reproduction/DATA_SOURCES.md` for exact inputs and commands. No raw or
cleaned market panels, downloaded archives or credential files are bundled.
The China datasets require the relevant CSMAR access. Downloading current
versions of external datasets does not guarantee identical historical files.
Main `reproduce_empirical()` runs automatically assemble tables and the alpha
distribution figure under `output_dir/report`. Use `report_empirical(output_dir)`
to rebuild this report from previously computed empirical outputs.

## Calibration is part of the method

The integrated manuscript reports two implementations. They must not be
silently interchanged:

| Convention | Simulation chapter | Empirical chapter |
|---|---|---|
| Candidate orders | 5 plus the full Ma dyadic grid | Same grid |
| Order 5 p-value | Add-one empirical bootstrap right tail | Same |
| Other component p-values | Normal tail using bootstrap mean and sample SD | Add-one empirical bootstrap right tail |
| Cauchy endpoint clipping | Machine epsilon | `1 / (N + T)` |
| Bootstrap studentization | OLS residual variance refitted each draw | OLS residual variance refitted each draw |

The later all-empirical-tail simulation revision is not substituted for this
manuscript's simulation results. The known-DGP, independent-null calibration
in the simulation is a diagnostic comparison, not an extra empirical-data
calibration procedure. Finite-sample size and power experiments do not establish
an exact fixed-level theorem for the complete growing candidate grid.

## Other methods and software

| Method | Implementation in this repository |
|---|---|
| L and AdaL | The manuscript's factor-model implementation, with explicit calibration choices |
| PY | A local implementation of the documented Pesaran--Yamagata correction |
| MAX | A local maximum-test implementation; finite-sample scaling follows the relevant chapter |
| COM | The paper's OR--Bonferroni combination of PY and MAX |
| max-k | The specified precision-transformed benchmark, inspired by Xia--Zhang; grid and precision estimation differ between simulation and empirical workflows |

These are **not calls to official author R packages**. The simulation's
covariance thresholding is implemented in base R; `cvCovEst::thresholdingEst`
was an optional validation reference in the original project, not a runtime
dependency. The empirical max-k implementation uses QUIC 1.2 source distributed
with `skggm` 0.2.8. Its original authorship, GPL-3 license and source history are
retained with the source. Python dependencies and their use are listed in the
empirical guide. See `inst/reproduction/METHODS.md` for references.

## Checks and reproducibility limits

```sh
R CMD build .
R CMD check --no-manual LFPM_0.1.0.tar.gz
```

The package separates deterministic implementation checks, short workflow
checks and full paper reproduction. A short smoke run is not evidence that the
full study has been rerun. Floating-point differences across R/BLAS platforms
can affect tail comparisons. The empirical backend retains NumPy's random
generator: replacing it by R's generator with the same integer seed does not
produce the same draws.

The archived final China max-k penalty-sensitivity values lack their original
seed and draw-budget metadata. The sensitivity runner performs a fresh,
deterministic calculation and labels this historical matching limitation;
it does not insert the paper's numbers into newly computed results.

GPL-3. See the preserved third-party license notices for the QUIC source.

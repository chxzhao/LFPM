# LFPM

An R package for *Testing factor pricing restrictions in high dimensions without knowing the sparsity of alphas*.
Version 0.1.2 provides reusable R methods, simulation reproduction, and empirical
reproduction code. **Actual case-study data are not included.** Obtain them from
[the external data sources](inst/reproduction/DATA_SOURCES.md), which also document
required files, fields, units and sample-construction rules.

Author: **Chenxi Zhao**.
Maintainer: Chenxi Zhao <chxzhao@mail.nankai.edu.cn>.
This is a research repository, not a CRAN release.

## Install

With R 4.2 or later, from this repository's root:

```sh
R CMD build .
R CMD INSTALL LFPM_0.1.2.tar.gz
```

The reusable R tests and simulations use R's standard packages. Empirical
reproduction uses an optional Python 3.9+ environment with NumPy, pandas,
SciPy, statsmodels and threadpoolctl; plotting additionally uses matplotlib.
The empirical max-k benchmark compiles the included QUIC source at runtime
and requires a C++ compiler with Accelerate on macOS or LP64 BLAS/LAPACK
on Linux. It is not currently supported on Windows; use
`include_maxk = FALSE` there to omit that comparator. Package installation
itself does not compile this optional solver. See the
[empirical guide](inst/reproduction/EMPIRICAL.md) for setup and tested versions.

## Test your own panel

Rows are observations (dates), columns of `Y` are asset excess returns, and
columns of `F` are factors. The functions add the intercept themselves.
This example uses synthetic data and shows the two component-calibration
options for the reusable R AdaL function.

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

`B = 99` is a quick example, not a formal reproduction budget. The
`"empirical"` option means empirical bootstrap-tail probabilities and is also
available as a general R statistical method. The paper's empirical reproduction
interface below preserves its original Python random streams; an identical
integer seed in R and NumPy does not generate identical draws.

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

## Reproduce the empirical analyses

Acquire the files described in [DATA_SOURCES.md](inst/reproduction/DATA_SOURCES.md)
and keep them outside the package and public repository. The China sources
require authorized CSMAR access. Public source vintages can change, so new
downloads do not guarantee the same historical observations.

```r
py <- "/path/to/your/python"

# One main setting, using externally acquired data.
china <- reproduce_empirical(
  data_dir = "/external/market-data",
  output_dir = "/external/results/china",
  analysis = "china", python = py)

# All supported final settings and checks; substantial computation.
# paper <- reproduce_empirical(
#   data_dir = "/external/market-data",
#   output_dir = "/external/results/paper",
#   analysis = "all", python = py)

# Rebuild tables and figures from an already computed output directory.
report <- report_empirical(china$output_dir)
```

`lfpm_empirical_test(Y, F, python = py)` exposes the empirical test definitions
on user-supplied matrices. `reproduce_empirical()` provides the China, U.S. FF3
and U.S. FF5 analyses, diagnostics, fitted-null calibration, break reselection
and max-k sensitivity. Main-setting runs automatically create a report with
joint tests, component probabilities, alpha magnitudes, stock-level BH counts,
available diagnostics, coverage flags and the alpha-distribution figure.
Missing settings and deliberately omitted comparators remain missing; historical
paper values are not inserted. Outputs may contain stock identifiers and
estimated alphas, so keep output directories outside the public repository.
See the [empirical guide](inst/reproduction/EMPIRICAL.md) for options, budgets,
seeds, complete inputs and output files. No raw or cleaned market panels,
downloaded archives, database exports or credentials are bundled.

## Calibration is part of the method

| Convention | Simulation implementation | Empirical implementation |
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
| MAX | Local maximum-test implementation; scaling follows the relevant workflow |
| COM | OR--Bonferroni combination of PY and MAX |
| max-k | The specified precision-transformed benchmark inspired by Xia--Zhang; simulation uses a four-order thresholded-covariance version, while empirical reproduction uses QUIC and Gaussian multipliers |

These are **not calls to official author R packages**. Simulation covariance
thresholding uses base R. `cvCovEst::thresholdingEst` was an optional numerical
validation reference in the source project and is not a runtime dependency.
The empirical solver uses QUIC 1.2 source obtained through `skggm` 0.2.8;
upstream authorship, license and source provenance are retained with it.
The benchmark is the documented adaptation, not a claim to reproduce every
published max-k procedure. See [method attribution](inst/reproduction/METHODS.md).

## Checks and reproducibility limits

```sh
R CMD build .
R CMD check --no-manual LFPM_0.1.2.tar.gz
```

Short numerical, synthetic-data and integration checks do not establish that
the full simulation and real-data study has been rerun or that every published
number has been revalidated. Floating-point differences across R/BLAS and
Python numerical libraries can affect tail comparisons. Simulation manifests
validate settings, budgets and engine compatibility before reusing checkpoints;
empirical runs record source-file hashes, software versions and seeds.

The archived final China max-k penalty-sensitivity endpoint lacks its original
seed and draw-budget metadata. The sensitivity interface performs a fresh,
deterministic calculation and flags `historical_rng_parity_verified = FALSE`.
It cannot certify exact regeneration of that historical Monte Carlo value.
The [empirical guide](inst/reproduction/EMPIRICAL.md) explains this limitation
separately from the recovered seed streams for the main tests.

License: GPL-3. The included QUIC source retains its upstream license notices.

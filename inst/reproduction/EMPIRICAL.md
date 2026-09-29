# Empirical reproduction

These workflows target *Testing factor pricing restrictions in high dimensions without knowing the sparsity of alphas*, including the final complete
Ma candidate grid. They rebuild the specified panels from externally supplied
raw CSVs and never require old result directories or saved bootstrap arrays.
See [DATA_SOURCES.md](DATA_SOURCES.md) for acquisition, fields, access, and vintages.

## Environment

Use Python 3.9+ with `numpy`, `pandas`, `scipy`, `statsmodels`, and `threadpoolctl`.
`matplotlib` is required only when `plots=TRUE`. A separate environment can be
prepared by the reader, for example:

```sh
python3 -m venv lfpm-python
# Activate this environment using your platform's venv instructions.
python -m pip install numpy pandas scipy statsmodels threadpoolctl matplotlib
```

The package does not install Python dependencies automatically. Pass the
environment's executable explicitly in `python=` if the system Python lacks
these packages. The initial portable-backend synthetic tests used Python 3.10.11,
NumPy 1.23.5, pandas 1.5.3, and SciPy 1.10.0. Provenance records the versions
actually used in every run.
The complete tested package versions are listed in
`python/requirements-tested.txt` inside the installed package (or
`inst/python/requirements-tested.txt` in the source repository). These pins
document the test environment; availability of compatible wheels depends on
the Python version and platform.

`include_maxk=TRUE` builds the included GPL-3 QUIC 1.2 source into
`output_dir/compiled` (or a setting subdirectory) on first use. macOS requires
a C++ compiler with Accelerate; Linux requires a C++ compiler and LP64 LAPACK
and BLAS development libraries. Only source, upstream attribution, and license
are distributed. No `.dylib`, `.so`, or other binary is shipped. The empirical
QUIC backend currently does not support Windows; use `include_maxk=FALSE` there
or run the full workflow on macOS/Linux. Omitted max-k results are absent, never
filled with archived numbers. The separate R implementation may use a different
solver; it should not be assumed bit-identical to this empirical backend.

## Call functions

```r
library(LFPM)
py <- "/path/to/lfpm-python/bin/python"

# Run only the main China analysis.
china <- reproduce_empirical(
  data_dir = "/external/market-data",
  output_dir = "/external/results/china",
  analysis = "china", python = py)

# Run main settings, diagnostics, fitted-null checks, break reselection,
# and max-k penalty sensitivity.
paper <- reproduce_empirical(
  data_dir = "/external/market-data",
  output_dir = "/external/results/paper",
  analysis = "all", python = py)

# Other choices: "us_ff3", "us_ff5", "break_reselection", "fitted_null",
# "diagnostics", and "sensitivity".
# Small B/repetitions are useful for checking the pipeline, not paper replication.
quick <- reproduce_empirical(
  "/external/market-data", "/external/results/quick", analysis = "us_ff3",
  python = py, B = 99, B_maxk = 99, include_maxk = FALSE)
```

For a caller's own data:

```r
# Y: T by N matrix of decimal monthly excess returns; N >= 40.
# F: T by p factor matrix, without intercept; matching months and complete data.
answer <- lfpm_empirical_test(Y, F, B = 9999, seed = 123,
  python = py, include_maxk = FALSE)
answer$p.value
answer$results
answer$components

# Fixed loading break: common intercept, each factor receives a post interaction.
F_piecewise <- cbind(F, F * as.numeric(months >= "2017-12"))
answer <- lfpm_empirical_test(Y, F_piecewise, python = py)
```

`lfpm_empirical_test()` does not choose the break. The paper's break-reselection
workflow is separate and repeats the specified score scan within every draw.

## Paper settings and conventions

| Workflow | Draws and seeds | Model/calibration |
|---|---|---|
| China main | 99,999; legacy base 20260903 plus SHA-256 seed key `2021-01_2025-12_top1000alpha1` | FF3; top 1,000 pre-window capitalization, 980 complete assets |
| U.S. FF3 | 19,999; base 20260906, key `US_piecewise_loading_2012_2019alpha1` | FF3 loading interactions from 2017-12; 288 assets |
| U.S. FF5 | 19,999; base 20260906, key `US_FF5_piecewise_L_AdaL` | Five loading interactions from 2015-07; same assets, FF5 RF |
| max-k | China uses 9,999 draws after its separate base-20260906 `CN_main_2021-01_2025-12_top1000alpha1` L stream; U.S. uses 9,999 on separate legacy seed keys | QUIC graphical lasso, Gaussian multiplier, K=5,10,30 |
| Main loading diagnostic | 1,999 draws; base 20260906, key `slope` followed by T and N | Constant-loading FF3, 20%-trimmed mean/max score scan, OR--Bonferroni |
| FF5 loading diagnostic | 4,999 draws; base 20260906, key `US_FF5_loading_break_scan` | Five-factor constant-loading score scan; records the average-score selected month |
| Break reselection | 9,999; seed `20260928 + B + 1000` | 20% trimmed FF3 score scan, minimum segment 12, batch 32; reselect in every draw |
| Fitted-null | 500 outer samples, 999 inner draws | Original outer/inner NumPy streams; paired frozen/refitted denominators |

The integer legacy seed is `(base + int(sha256(key)[:8], 16)) mod 2^32`.
NumPy `default_rng` uses PCG64; R's `set.seed()` does not set these streams.
The main bootstrap preserves the original batch-128 arithmetic; fitted-null
uses batch 256 and the original `bar_h` numerator arithmetic. CPU libraries can
still change the last numerical bits.

The Chinese AdaL and max-k entries in the final manuscript retain different
archived runs. AdaL uses the earlier 99,999-draw L stream; max-k's reported
0.00110 comes from the later 9,999-draw bimarket run. Using the former run's
max-k result (0.00075) would not reproduce the final comparator table.

All empirical AdaL components use `(1 + exceedances)/(B + 1)`, equal weights,
the complete dyadic grid, and endpoint bound `1/(N+T)`. The final main-study
probabilities fall inside that bound, so its raw and bounded formulas agree.
Residual variances are recomputed for each bootstrap draw. T1 is reported
separately. The simulation chapter instead uses a bootstrap-moment normal
tail for the dyadic components; do not silently merge those calibrations.

PY, MAX, and max-k are the study's own adaptations of published methods, not
calls to the original authors' official packages. Numerical services are
provided by SciPy/NumPy; QUIC provides the graphical-lasso solver. `COM` is
`min(1, 2*min(PY_p, MAX_p))`, the study's OR--Bonferroni rule. It is not an
arithmetic mean or a claim to reproduce Feng's literal original COM definition.
QUIC authorship and license are in `python/vendor/quic/SOURCE.md` and `LICENSE`.

## Outputs and verification scope

Main-setting calls (`all`, `china`, `us_ff3`, `us_ff5`) automatically run
`report_empirical(output_dir, plots=plots)`. Its `report/` subdirectory contains
joint and component tables, economic magnitudes, available diagnostics,
source-file checksums, a manuscript mapping, coverage flags, and the combined
China/U.S. alpha figure when requested. Missing settings are labelled rather
than filled with historical results. Rebuild this report with
`report_empirical(output_dir)` without recomputing the bootstrap.

Main runs write per-setting `results.csv`, `components.csv`, `alphas.csv`, and
`metadata.json`, plus combined `main_results.csv` and
`component_probabilities.csv`. Optional residual diagnostics reproduce BP(6),
BP(12), LB(12), and ARCH-LM(6) with separate BH adjustments; ARCH uses the full
design rank as `ddof`. The main settings also write `model_diagnostics.csv`:
the descriptive complex-half-sample eigenvalue-gap ratio and the constant-loading
model's average/maximum score-scan probabilities and selected month. `B` overrides
the scan's 1,999-draw FF3 (4,999-draw FF5) budget too. `economic_magnitudes.csv` contains alpha quantiles,
absolute magnitudes, and BH counts at 5%/10%. Optional PDF figures show estimated
and standardized alpha distributions.
The figures use the recomputed values but do not promise the paper's exact
typographical layout.

Fitted-null runs write raw probabilities, component probabilities, rejection
rates, Monte Carlo standard errors, and Wilson intervals. Break-reselection
runs write paired fixed/reselected probabilities, components, the observed
score scan, bootstrap break distribution, frequency of reselecting the observed
break, and central 95% selected-date quantiles. `selected_stock_paths.csv` and
optional PDF use the eight largest absolute standardized alphas, with raw-return
or adjusted-price indices normalized to 100 in the first month. These paths are
descriptive, not factor-adjusted investment performance. Data-derived outputs may
contain stock identifiers and estimated alphas: store them outside the public
repository. `provenance.json` and `python.log` document every successful run;
failed runs keep their logs. Named output files are replaced when reusing an
output directory.

## Penalty-sensitivity provenance limitation

`analysis="sensitivity"` and full `analysis="all"` recompute max-k at penalty
multipliers 0.75, 1, and 1.25 for each constant-loading FF3 panel, reporting
K=5,10,30 and solver diagnostics. Their CSV and JSON record the explicit
`sensitivity_seed` base and SHA-256-derived keys. These are recomputations,
not copied summary values.

The final archived Chinese endpoint 0.0264 occurs in the original
`main_sensitivity_results.csv`, but available research scripts only consume that
table; its generating command, draw count, and seed were not archived. The older
script with a documented seed reports 0.0282. The package therefore cannot certify
exact regeneration of the 0.0264 Monte Carlo value. It provides a deterministic
rerun of the stated procedure, flags `historical_rng_parity_verified=FALSE`, and
allows a recovered seed convention to be incorporated later. This limitation is
separate from the main tests, whose seed streams were recovered.
`include_maxk=FALSE` explicitly omits solver-dependent sensitivity from `all` and
records that omission in provenance.

The initial package work verifies synthetic-data numerical behavior, source
algorithm preservation, runtime QUIC compilation, and R/Python plumbing. It
does **not** rerun the full real-data calculations or claim to have revalidated
all published probabilities against downloaded market data. Historical
exploratory windows, outdated two/three-component AdaL reports, and abandoned
calibrations are not part of `analysis="all"`. The maintained interface covers the
final main settings, residual/loading/gap diagnostics, economic magnitudes and BH
counts, selected-stock paths, fitted-null checks, FF3 break reselection, and max-k
penalty sensitivity. The sensitivity's missing historical RNG provenance is
disclosed above. It does not reproduce abandoned exploratory pretests or the
manuscript's page layout.

The empirical evidence remains conditional on inspected stock pools, windows,
and model specifications. Fitted-null size diagnostics are not a theorem for
the complete growing AdaL grid.

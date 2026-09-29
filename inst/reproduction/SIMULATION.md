# Chapter 3: standalone simulation reproduction

This package follows *Testing factor pricing restrictions in high dimensions without knowing the sparsity of alphas* and its Ma-rule AdaL simulation,
from `simulation_ma2024_20260929` plus the preserved four comparators in
`simulation_five_methods_20260927`. It does not use the subsequent alternative
`sim_fullgrid_refit` experiment. No original project directory, archive,
precomputed raw checkpoint, Python installation, or third-party R package is
needed for the simulation. It uses base/recommended R packages only.

## Small functional run

```r
library(LFPM)
smoke <- reproduce_simulation("simulation-smoke")
plot_simulation(smoke, "simulation-smoke/power.pdf")
plot_simulation(smoke, "simulation-smoke/size.pdf", phase = "size")
```

The default is explicitly a smoke run: one covariance model, Normal errors,
N=40, T=30, m in 1,5,20,40; R=6, Rcal=39, B=39. These are not manuscript
estimates. Its small null reference is only a software exercise and has coarse
rank p-values.

## Full manuscript computation

```r
paper <- reproduce_simulation("simulation-paper", preset = "paper", workers = 4)
```

The paper preset includes all four models, Normal/t5/MixNorm/Chisq errors,
six size dimensions (T in 100,200 crossed with N in 100,200,400), and power at
T=100, N=200 for **every m=1,...,200**. It uses R=1000 independent evaluation
repetitions, Rcal=1999 independent calibration repetitions, B=999 residual
wild-bootstrap draws, seed=20260927, signal energy=10 and strict p<0.05.
This yields 96 size scenarios and 3200 power scenarios, with five methods and
both native/equal-size calibration. The complete job is expensive and can
run for many hours depending on hardware; do not infer a full reproduction
from a successful smoke run.

To inspect/select designs before committing compute:

```r
design <- simulation_design("paper", models = 1, distributions = "Normal",
                            size_N = 200, size_T = 100, m = c(1, 5, 20, 100, 200))
design
subset <- reproduce_simulation("simulation-subset", design = design)
```

A supplied paper design still defaults to the full repetition/bootstrap
budgets. Explicit `R`, `Rcal`, `B` overrides are allowed, but then this is a
modified-budget experiment. Omitted m under the paper preset always means all
integers through `power_N`. If m is restricted, band summaries use exactly
the selected values and expose their count; they do not reconstruct missing
points. Asset counts must be at least 40 for the Ma candidate set.

## Reproducibility and checkpoints

Rerun the **same call and output directory** to resume. Budgets, grid, seed,
RNG kind, R version, and engine definitions must match the saved manifest.
Incompatible reuse raises an error before sampling. Use a new directory for
different budgets or scientific settings. `workers`, `chunk`, `quiet`, and
`phase` may change without changing random draws. Do not run concurrent jobs
in the same output directory. Windows runs sequentially; macOS/Linux can use
fork workers.

`phase="calibration"` runs only independent null references. `phase="size"`
and `phase="power"` compute their needed references and evaluations; `all`
does both. Summary CSVs are updated after completed scenarios, while RDS
checkpoints are updated after each chunk (default 25 repetitions). An
interruption loses at most the current chunk. Keep the manifest and
checkpoints together. `resume=FALSE` rejects an existing run instead of
overwriting it.

For PY/MAX/COM/max-k, the retained engine is recomputed. Ma-rule AdaL is then
computed separately after restoring the same per-repetition seed. Both see
the same panel and bootstrap draws, reproducing the original combined
experiment without depending on old comparator checkpoints. Seeds use the
original `calibration`, `size`, and `power_curve` namespaces and setting keys.
Power draws a full support permutation **before** bootstrap sampling, then
uses nested supports and common random numbers across all m. Changing the
selected m grid therefore does not change an included point's random panel.
All bootstrap variances are full OLS refits. The previous RNG state of the R
session is restored when reproduction returns or errors.

Bitwise identity across different BLAS libraries, processors, or R versions
is not promised. The original numerical definitions and RNG sequence are
preserved. The factor recursion intentionally retains the manuscript code's
innovation-square recursion; it is not silently replaced by standard GARCH.

## Tables, uncertainty and figures

The runner writes `size_summary.csv`, `power_summary.csv`, their
`*_paired_differences.csv` files, `power_band_summary.csv` and
`power_band_paired_differences.csv`. It also saves `manifest.rds`,
`sessionInfo.txt`, `progress.rds`, and all completed raw chunks in checkpoints.
When the corresponding phases finish it derives the manuscript-format
`size_native_percent.csv`, `size_equal_size_percent.csv`,
`power_selected_m_percent.csv`, and `manuscript_claims.csv` from the new
summaries. Selected-m tables retain only requested points among
1,5,10,20,50,100,200; claims retain the original descriptive comparator rule.

Native results use each method's p-value. Equal-size results use add-one
right-tail ranks against a **separate** null calibration sample, with ties
counted conservatively. Neither size nor power is evaluated on its own
calibration sample. Pointwise rates use binomial MCSE and Wilson 95% intervals.
AdaL-minus-comparator differences are paired by independent repetition. Band
averages (1:10, 11:50, 51:200, 1:200) are formed within each repetition before
computing standard errors, retaining dependence across m. Equal-size
intervals condition on the realized calibration sample. These are pointwise
intervals, not simultaneous or comparator-selection-adjusted intervals.

```r
plot_simulation(paper, "power_main_equal_size.pdf",
                distributions = c("Normal", "t5"))
plot_simulation(paper, "power_supplement_equal_size.pdf",
                distributions = c("MixNorm", "Chisq"))
plot_simulation(paper, "power_main_native.pdf", calibration = "native",
                distributions = c("Normal", "t5"))
plot_simulation(paper, "power_supplement_native.pdf", calibration = "native",
                distributions = c("MixNorm", "Chisq"))
plot_simulation(paper, "power_all_equal_size.pdf")
plot_simulation(paper, "size_all_native.pdf", phase = "size", calibration = "native")
plot_simulation(paper, "size_all_equal_size.pdf", phase = "size")
```

The base-R PDFs reproduce the underlying curves/rates and Wilson intervals.
Their typography is not intended to be byte-identical to the manuscript PDFs.

## Cached manuscript references: zero computation

```r
reference <- reproduction_reference()
reference$validation_summary
reference$power_summary
reference$manuscript_claims
plot_simulation(reference, "reference-power.pdf")
plot_simulation(reference, "reference-size.pdf", phase = "size")
```

This only reads bundled **cached aggregate tables** from the final manuscript
simulation. It does not simulate new draws. The plot titles and R provenance
fields explicitly mark cached results. The reference tables include all
32000 power summary rows and 960 size summary rows, paired uncertainties,
band means, selected-m and size percentage tables, and descriptive manuscript
claims. Source checksums and budgets are in `reference/simulation/PROVENANCE.json`.
No individual panels, raw Monte Carlo arrays or empirical case-study data
are included. Directly selecting the numerically largest comparator in
`manuscript_claims` is descriptive; its intervals are not selection-adjusted.

## Comparator names

PY, MAX, COM and the AT_XZ-style max-k comparator are the supplied paper-code
implementations. No third-party method package is installed or called by this
runner. The max-k precision threshold reproduces the original code's
thresholding operation (previously checked against
`cvCovEst::thresholdingEst`), including its original ridge shift. The package
does not claim that this preserved comparator is an exact implementation of
every published max-k variant. Consult the method documentation and package
references for the statistical definitions and attribution.

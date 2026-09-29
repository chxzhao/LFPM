\name{simulation_design}
\alias{simulation_design}
\alias{reproduce_simulation}
\alias{plot_simulation}
\alias{reproduction_reference}
\title{Reproduce the September 29, 2026 Manuscript Simulation}
\description{
Construct a simulation design, recompute all five methods from fresh random
draws with restartable checkpoints, read cached manuscript summary tables,
or plot results using base R graphics.
}
\usage{
simulation_design(preset = c("smoke", "paper"), models = NULL,
  distributions = NULL, size_N = NULL, size_T = NULL,
  power_N = NULL, power_T = NULL, m = NULL)

reproduce_simulation(output_dir, preset = c("smoke", "paper"),
  design = NULL, R = NULL, Rcal = NULL, B = NULL,
  seed = 20260927L, workers = 1L, chunk = 25L,
  resume = TRUE, phase = c("all", "calibration", "size", "power"),
  quiet = FALSE)

reproduction_reference(table = c("all", "size_summary", "power_summary",
  "size_paired_differences", "power_paired_differences",
  "power_band_summary", "power_band_paired_differences",
  "size_native_percent", "size_equal_size_percent",
  "power_selected_m_percent", "manuscript_claims", "validation_summary"))

plot_simulation(x, file = NULL, phase = c("power", "size"),
  calibration = c("equal_size", "native"), models = NULL,
  distributions = NULL, width = NULL, height = NULL)
}
\arguments{
  \item{preset}{\code{"smoke"} is a small functional check, not a scientific
  reproduction. \code{"paper"} selects the manuscript design and budgets.}
  \item{models}{Covariance model numbers selected from \code{1:4}.}
  \item{distributions}{Innovation distributions selected from
  \code{c("Normal", "t5", "MixNorm", "Chisq")}.}
  \item{size_N, size_T}{Asset counts and time lengths for independent size
  evaluation. All combinations are used. Defaults: paper, \code{c(100,200,400)}
  and \code{c(100,200)}; smoke, 40 and 30.}
  \item{power_N, power_T}{Single asset count and time length for power curves.
  Defaults: paper, 200 and 100; smoke, 40 and 30.}
  \item{m}{Positive nonzero-intercept counts. Paper defaults to every integer
  from 1 through \code{power_N}; smoke defaults to 1, 5, 20, 40, capped at
  \code{power_N}. Asset counts must be at least 40 for the Ma candidate grid.}
  \item{output_dir}{An empty output directory for a new run, or the same
  directory and configuration for resumption. Do not run concurrent jobs in
  the same directory.}
  \item{design}{A \code{simulation_design()} result. Its preset determines
  default budgets. To choose scenarios, supply a restricted design.}
  \item{R, Rcal, B}{Independent evaluation repetitions, independent null
  calibration repetitions, and wild-bootstrap draws per repetition. Defaults
  are 1000, 1999, 999 for paper and 6, 39, 39 for smoke. Each is at least 2.
  Changing these values requires a new output directory.}
  \item{seed}{Nonnegative integer seed, default 20260927.}
  \item{workers}{Fork workers on macOS/Linux. Windows uses sequential
  execution. Random seeds do not depend on worker count.}
  \item{chunk}{Repetitions between checkpoint writes.}
  \item{resume}{Whether to resume a compatible output directory.}
  \item{phase}{For reproduction: all, calibration only, size with calibration,
  or power with calibration. For plotting: power or size.}
  \item{quiet}{Suppress progress messages.}
  \item{table}{Cached aggregate reference table to load, or all tables.}
  \item{x}{A reproduction result, a reference-table list, or a summary data
  frame with the corresponding columns.}
  \item{file}{Optional PDF filename. \code{NULL} plots on the current device.}
  \item{calibration}{Native method p-values or independent-null rank
  calibration (\code{"equal_size"}).}
  \item{width, height}{Optional PDF dimensions in inches.}
}
\details{
The full design has 96 size scenarios and 3200 power scenarios (four covariance
models, four distributions, six size dimensions, and all 200 sparsities at
T=100 and N=200). Rejection uses the strict rule \eqn{p < 0.05}.
The five methods are PY, MAX, COM, the retained AT\_XZ-style max-k comparator,
and the Ma-rule AdaL used in \emph{LFPM\_revised\_20260929}. This is the
chapter 3 experiment, not the later alternative full-grid simulation.

The four earlier comparator statistics and Ma-rule AdaL are recomputed with
the original phase/scenario/repetition seed keys. Both engines restart from
the same seed, preserving identical panels, support draws, and bootstrap
multipliers. Power uses common random numbers across m with nested uniform
supports, while calibration, size, and power use separate seed namespaces.
All bootstrap variances use full OLS refits. Signal energy is fixed at 10,
with all nonzero intercepts positive. No external checkpoint is required.

Checkpoints are written after each chunk and include repetition identifiers.
The output manifest validates the complete budgets, design, seed, RNG kind,
R version, and engine definitions before any saved results are reused.
Workers, chunk size, and requested phase may change on resume; scientific
configuration may not. R's previous RNG state is restored after the call.
Exact floating-point identity across R versions or BLAS implementations is
not guaranteed, even with the same random draws.

Pointwise intervals are Wilson 95 percent intervals and MCSE is the binomial
standard error. AdaL differences are paired by repetition. Power-band means
are first averaged within repetition, preserving dependence across m; only
selected m values in each named band contribute, and \code{m_count} records
that count. Equal-size uncertainty is conditional on the fixed independent
calibration reference; it does not include all uncertainty from estimating
that reference. Intervals are pointwise, with no multiple-comparison or
comparator-selection adjustment.

\code{reproduction_reference()} only reads bundled cached CSV summaries and
does not rerun any simulation. Its provenance field/attribute and plot titles
identify this distinction. These are aggregate simulation results, not raw
simulation arrays or empirical financial data. Full reproduction is a large
computation; the smoke preset and examples use much smaller budgets.
}
\value{
\code{simulation_design()} returns an \code{lfpm_simulation_design} list with
size and power setting data frames and the sparsity grid. A reproduction
returns an \code{lfpm_simulation} list of summary, paired-difference and
band-summary data frames, derived manuscript percentage/claim tables,
the output directory, manifest and provenance.
The same tables are saved as CSV along with RDS checkpoints and session info.
\code{reproduction_reference()} returns a data frame or named list.
\code{plot_simulation()} invisibly returns the plotted data frame.
}
\examples{
design <- simulation_design("paper", models = 1, distributions = "Normal")
nrow(design$size)
length(design$m)
reference <- reproduction_reference("size_summary")
head(reference)
\dontrun{
# Small end-to-end functional check; not paper results.
smoke <- reproduce_simulation(file.path(tempdir(), "lfpm-smoke"))
plot_simulation(smoke, file = file.path(tempdir(), "smoke-power.pdf"))

# All manuscript settings and exact manuscript budgets; substantial runtime.
paper <- reproduce_simulation("simulation-paper", preset = "paper", workers = 4)
# Resume with the identical scientific configuration.
paper <- reproduce_simulation("simulation-paper", preset = "paper", workers = 2)

# A selected paper scenario; still defaults to R=1000, Rcal=1999, B=999.
one <- simulation_design("paper", models = 1, distributions = "Normal",
                         size_N = 200, size_T = 100)
selected <- reproduce_simulation("simulation-selected", design = one)

# Plot cached manuscript tables without running simulations.
ref <- reproduction_reference()
plot_simulation(ref, file = "paper-power.pdf")
plot_simulation(ref, file = "paper-size.pdf", phase = "size")
}
}

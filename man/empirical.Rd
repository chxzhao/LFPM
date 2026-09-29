\name{lfpm_empirical_test}
\alias{lfpm_empirical_test}
\alias{reproduce_empirical}
\title{Empirical Factor-Model Tests and Paper Reproduction}
\description{Runs the portable Python implementation of the empirical analyses
in the final LFPM manuscript. Market data must be supplied separately.}
\usage{
lfpm_empirical_test(Y, F, B = 9999L, seed = 20260929,
  python = Sys.which("python3"), include_maxk = FALSE, B_maxk = B,
  output_dir = tempfile("lfpm-empirical-"))

reproduce_empirical(data_dir, output_dir,
  analysis = c("all", "china", "us_ff3", "us_ff5",
    "break_reselection", "fitted_null", "diagnostics", "sensitivity"),
  python = Sys.which("python3"), B = NULL, B_maxk = 9999L,
  repetitions = 500L, inner_B = 999L, sensitivity_seed = 20260906,
  include_maxk = TRUE,
  diagnostics = TRUE, plots = TRUE)
}
\arguments{
\item{Y}{Finite numeric matrix: rows are months, columns are assets, in decimal
excess-return units. At least 40 assets are required.}
\item{F}{Finite numeric matrix of factors, aligned with the rows of \code{Y}.
Do not include an intercept. For fixed loading breaks, append factor-by-post-break
interactions; one common intercept is added internally.}
\item{B}{Number of bootstrap draws. For reproduction, \code{NULL} uses 99999 for
China, 19999 for U.S. main analyses, and 9999 for break reselection. A non-NULL
override is for shorter runs and does not reproduce the paper's Monte Carlo budget.}
\item{seed}{Nonnegative NumPy PCG64 seed for a custom test. R's random seed is not used.}
\item{python}{Python 3.9 or later executable with numpy, pandas, scipy, statsmodels,
and threadpoolctl. matplotlib is additionally required for plots.}
\item{include_maxk}{Whether to compute the QUIC precision-transformed max-k benchmark.
Requires a C++ compiler and BLAS/LAPACK on macOS or Linux. The empirical QUIC
backend does not currently support Windows.}
\item{B_maxk}{Draws for max-k and its penalty sensitivity. China preserves the
final comparator's separate bimarket stream, including its preceding L draws,
rather than the main AdaL stream. The paper uses 9999.}
\item{output_dir}{Directory for tables, provenance, logs, and generated figures.
The first max-k call also compiles QUIC into a subdirectory here. Existing files
with the same output names are replaced. Keep this directory outside the package.}
\item{data_dir}{External directory with \code{China/CSMAR} and/or \code{America}
subdirectories. See the installed \file{reproduction/DATA_SOURCES.md}.}
\item{analysis}{Paper setting or workflow. \code{"all"} includes three main settings,
paired fitted-null diagnostics, U.S. break reselection, and max-k penalty
sensitivity when \code{include_maxk=TRUE}. The \code{"diagnostics"} choice runs
the two main samples' residual and specification diagnostics only.}
\item{repetitions}{Outer fitted-null repetitions; the paper uses 500.}
\item{inner_B}{Inner fitted-null draws; the paper uses 999.}
\item{sensitivity_seed}{Nonnegative seed base for the documented penalty-sensitivity
rerun. Original penalty-sensitivity seeds were not archived, so this is an
explicit new convention and not a claim of identical Monte Carlo draws.}
\item{diagnostics}{Write the four residual-diagnostic families, BH counts, and
the main settings' eigenvalue-gap and loading-stability score diagnostics.}
\item{plots}{Write alpha-distribution PDF figures.}
}
\details{Every empirical AdaL component uses an add-one empirical upper tail,
clipped to \eqn{[1/(N+T),1-1/(N+T)]}, and equal-weight Cauchy aggregation.
Residual variances are re-estimated in each bootstrap draw. This differs from
the simulation chapter's bootstrap-moment normal calibration for dyadic components.
The maximum statistic T1 is reported separately from the AdaL candidate grid.
COM is this paper's OR--Bonferroni combination of PY and MAX.

The manuscript's fixed U.S. loading breaks are 2017-12 for FF3 and 2015-07 for
FF5. The break-reselection workflow repeats the trimmed FF3 score scan in every
null sample. It conditions on one common loading break, rather than selecting
between zero and multiple breaks. Fitted-null diagnostics pair frozen and
refitted denominators on the original NumPy seed streams.

No market observations, prepared panels, or Python binary libraries are
distributed. CSMAR requires appropriate database access. Public data revisions,
NumPy versions, and numerical libraries can affect exact numeric replication;
the output records input hashes and software versions. The full market-data
workflows have not been rerun as part of the initial package assembly.}
\value{\code{lfpm_empirical_test} returns a list containing the AdaL statistic and
p-value, all requested method results, component probabilities, asset-level
alpha estimates, and output paths. \code{reproduce_empirical} returns the analysis
name, output directory, generated file list, provenance path, and a report
assembled by \code{\link{report_empirical}} for main-setting runs. Other
analyses return \code{report=NULL} and keep their dedicated backend outputs.}
\examples{
\dontrun{
set.seed(1)
F <- matrix(rnorm(60 * 3), 60, 3)
Y <- F %*% matrix(rnorm(3 * 40), 3, 40) + matrix(rnorm(60 * 40), 60, 40)
answer <- lfpm_empirical_test(Y, F, B = 99)
answer$p.value

reproduce_empirical("/external/market-data", "/external/lfpm-results",
  analysis = "all", python = "/path/to/python")
}
}

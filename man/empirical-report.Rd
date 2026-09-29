\name{report_empirical}
\alias{report_empirical}
\title{Build Manuscript Tables and Alpha Histograms from Empirical Output}
\description{
Read locally computed empirical backend results and create joint-test,
ordered-component, economic-magnitude and stock-level summaries, together
with the manuscript's four-panel alpha-distribution figure. No financial
input data or historical manuscript values are bundled or substituted.
}
\usage{
report_empirical(output_dir, plots = TRUE)
}
\arguments{
  \item{output_dir}{The existing output directory created by
  \code{reproduce_empirical()}. Recognized subdirectories are \code{china},
  \code{us_ff3}, and \code{us_ff5}.}
  \item{plots}{Create a PDF with China and US FF3 monthly-alpha and
  standardized-alpha histograms.}
}
\details{
Each setting is reported only if its \file{results.csv},
\file{components.csv}, and \file{alphas.csv} files all exist. Incomplete or
absent settings are recorded explicitly in \file{coverage.csv}; they are
not filled from the manuscript. Malformed present files cause an error.
Skipped max-k computations and unavailable selected component orders are
also recorded. Filenames do not prove that a run used the final manuscript's
sample, model or budget; retain and inspect backend metadata and provenance.

The joint-test table selects PY, MAX, COM, max-k with 30 candidate orders,
and AdaL. The selected component table includes T1 from the separately
bootstrapped maximum and T5, ceiling(0.25N), and ceiling(0.50N) from the
saved component table. A full-component table preserves all computed
orders, empirical probabilities, clipped probabilities and weights.
US FF5 results are distinguished from the China/US FF3 main comparison.
The backend records method-specific bootstrap budgets. Missing B is allowed
for the analytical PY, MAX and COM methods and remains missing in the report;
AdaL, T1 and max-k require finite positive integer budgets. If an analytical
row includes a budget from an older output format, it must also be valid.

Alpha magnitudes include mean, sample standard deviation, median, fifth
and ninety-fifth percentiles, mean absolute value, maximum absolute value,
and the share with absolute alpha at least 50 basis points. All magnitude
fields are in monthly basis points; the share is a fraction in [0,1].
Quantiles use R type 7 interpolation, matching pandas' default linear rule.
Benjamini--Hochberg discovery counts use q-values less than or equal to
0.05 and 0.10. Saved q-values are checked against a fresh BH adjustment of
the individual p-values. Decimal and basis-point alphas must agree.

The PDF has two panels per main market: monthly basis-point alphas and
standardized alphas. Dashed lines at -1.96 and 1.96 appear only on standardized
panels and are pointwise normal references, not multiplicity-adjusted
critical values. If one main setting is missing its panels say so.
No plot is generated if both main-market settings are incomplete or missing.

The generated \file{manuscript_mapping.csv} maps files and columns to
\code{emp:tab:joint}, \code{emp:tab:magnitudes}, \code{emp:fig:alphas}, the
stock-level evidence, full-grid component probabilities, and FF5 robustness
paragraphs in \emph{LFPM\_revised\_20260929}. Residual diagnostics, break
reselection, fitted-null calibration, and max-k sensitivity are not recomputed
by this reporting function. When optional \file{model_diagnostics.csv} and
\file{diagnostic_summary.csv} exist, they are summarized in
\file{diagnostic_table.csv}, mapping to \code{emp:tab:diagnostics}; coverage
flags explicitly identify missing diagnostics. The static gap ratio is
descriptive. Loading stability uses the constant-loading specification,
whereas residual diagnostic counts use the fitted main specification.
}
\value{
An \code{lfpm_empirical_report} list containing data frames named
\code{coverage}, \code{joint_tests}, \code{selected_components},
\code{full_components}, \code{alpha_magnitudes}, \code{source_files}, and
\code{manuscript_mapping}, optional \code{model_diagnostics},
\code{residual_diagnostic_summary}, and \code{diagnostic_table};
the report directory, figure path, a
\code{main_complete} indicator, and provenance. CSV tables and
\file{REPORT.txt} are saved under \file{output_dir/report/}. Source-file MD5
checksums bind the report to the CSVs used. The optional PDF is
\file{report/figures/Figure_4_2_alpha_histograms.pdf}.
}
\seealso{\code{\link{reproduce_empirical}}}
\examples{
\dontrun{
# All numerical inputs here are newly produced backend outputs.
run <- reproduce_empirical("external-data", "empirical-run", analysis = "china")
report <- report_empirical(run$output_dir)
report$coverage
report$alpha_magnitudes
report$selected_components
}
}

\name{lfpm_methods}
\alias{lfpm_methods}
\alias{lfpm_fit}
\alias{simulate_lfpm}
\alias{lfpm_candidate_grid}
\alias{l_test}
\alias{adal_test}
\alias{lfpm_test}
\title{Studentized Intercept Tests for Linear Factor Pricing Models}
\description{
Fit a factor pricing model, simulate the paper's three-factor designs, and
compute fixed-order or adaptive L-statistics. The five-method interface uses
the simulation definitions accompanying the integrated LFPM manuscript.
}
\usage{
lfpm_fit(Y, F)

simulate_lfpm(model, N, T_size, distribution = "Normal", m = 0L,
              signal_energy = 10, seed = NULL)

lfpm_candidate_grid(N)

l_test(Y, F, k, B = 999L, seed = NULL, multipliers = NULL)

adal_test(Y, F, B = 999L, calibration = c("simulation", "empirical"),
          seed = NULL, multipliers = NULL)

lfpm_test(Y, F, B = 999L, seed = NULL)
}
\arguments{
\item{Y}{Finite numeric matrix or numeric data frame, with observations in rows
and assets in columns. Missing values are not omitted.}
\item{F}{Finite numeric factor matrix or numeric data frame with the same rows
as \code{Y}, excluding the intercept. A matrix with zero columns specifies an
intercept-only model. The design including its intercept must have full rank
and positive residual degrees of freedom.}
\item{model}{Covariance model, an integer from 1 to 4.}
\item{N}{Number of assets. Simulation requires at least 3; the adaptive candidate
grid requires at least 40. The fixed-order test supports one or more assets.}
\item{T_size}{Number of observations; at least 5 for the three-factor simulator.}
\item{distribution}{Innovation law: \code{"Normal"}, \code{"t5"},
\code{"MixNorm"}, or \code{"Chisq"}. Each is standardized to unit variance
before the covariance transformation.}
\item{m}{Number of nonzero intercepts, an integer between 0 and \code{N}.}
\item{signal_energy}{Finite nonnegative signal constant. For positive \code{m},
each active intercept is \eqn{\sqrt{\mathrm{signal\_energy}\log(N)/(T m)}}.
Supports are selected uniformly without replacement.}
\item{seed}{Optional nonnegative integer. Supplying a seed preserves the
caller's complete prior \code{.Random.seed}, including its absence, even on
error. With \code{NULL}, random draws use and advance the current RNG stream.}
\item{k}{Fixed aggregation order, one integer from 1 to the number of assets.}
\item{B}{Number of residual wild-bootstrap draws, an integer of at least 2.}
\item{multipliers}{Optional numeric matrix with \code{T_size} rows and \code{B}
columns, containing only -1 and 1. Each column is a bootstrap draw; the same
time multiplier is applied to every asset. When supplied, no multiplier draws
are generated.}
\item{calibration}{Marginal calibration for AdaL. \code{"simulation"} uses the
paper's Ma-rule simulation implementation. \code{"empirical"} uses empirical
bootstrap tails for every component. This argument does not request independent
null simulation or size-adjusted diagnostic calibration.}
}
\details{
The coordinate statistic is the fitted intercept times \eqn{\sqrt{q}} divided
by its residual standard deviation, where \eqn{q=1'M_F1}. Variances use
\eqn{T-p-1} residual degrees of freedom. Every wild-bootstrap draw refits the
intercept and all factors and recomputes the residual variances, using the
equivalent projection identity. For \code{l_test} and simulation AdaL,
nonpositive bootstrap residual sums of squares raise an error and no variance
floor is applied.

The fixed-order statistic sums the largest \code{k} squared coordinate
statistics. Its p-value is \eqn{(1+\#\{L_k^*\geq L_k\})/(B+1)}.

The adaptive grid is \eqn{\{5\}\cup\{\lceil N/2^j\rceil:
j=1,\ldots,\lfloor\log_2(N/20)\rfloor\}}. Its returned order is 5 followed
by the proportional orders from largest to smallest. The simulation AdaL
uses an add-one empirical upper tail for order 5 and normal upper tails
centered and scaled with bootstrap means and sample standard deviations for
the proportional orders. Before equal-weight Cauchy combination, probabilities
are clipped only at double-precision machine epsilon. This is the manuscript's
Ma-rule simulation, not the later all-empirical fullgrid/refit experiment.

Empirical AdaL instead uses add-one empirical upper tails, including ties,
for every candidate, and clips at \eqn{1/(N+T)} and \eqn{1-1/(N+T)}. In both
cases the analytic Cauchy upper-tail probability is evaluated using
\code{pcauchy} for numerical stability. P-values are returned without applying
a rejection threshold.

The \code{calibration="empirical"} option applies a
\eqn{10^{-20}} lower bound on refitted bootstrap variances, including zero
or negative values from numerical subtraction. Its \code{variance_audit}
reports the smallest raw variance and the number of entries at or below this
bound. This safeguard is a numerical convention of this calibration; it is
not the zero-statistic convention discussed in the mathematical theory, and
it is not applied to the simulation branch. Nonfinite numerical results raise
an error. With very small data scales this absolute floor can affect results,
so callers should inspect the audit fields.

\code{lfpm_test} returns PY, MAX, COM, max-k, and simulation AdaL, sharing one
multiplier matrix. PY includes the original correlation correction, MAX its
finite-sample scaling, COM is OR--Bonferroni, and max-k retains the thresholded
precision estimate and four-point grid from the simulation code. PY requires
residual degrees of freedom greater than four. These are the paper's simulation
benchmarks and are not claimed to duplicate every published algorithm.
These are general-purpose methods for caller-supplied matrices. Case-study
analysis and reproduction scripts are not included; only data-source links
are provided in \file{reproduction/DATA_SOURCES.md}.

The simulator retains the frozen factor recursion, innovation order, covariance
construction, 50-observation burn-in, and random loading and support draws. It
returns one panel, not a Monte Carlo study.
}
\value{
\code{lfpm_fit} returns an \code{lfpm_fit} list containing intercepts
\code{alpha}, residuals \code{U}, residual variances \code{s2}, residual
degrees of freedom \code{df}, \code{q}, projection quantities \code{P} and
\code{Q}, dimensions \code{T_size} and \code{N}, and coordinate statistics
\code{tstat}.

\code{simulate_lfpm} returns matrices \code{Y} and \code{F} in a list.
\code{lfpm_candidate_grid} returns an integer vector.

\code{l_test} and \code{adal_test} return \code{htest} objects with
\code{statistic}, \code{p.value}, parameters, and calibration metadata.
The fixed-order result includes \code{bootstrap.statistics}; the adaptive
result includes a \code{components} data frame with candidate order, observed
statistic, raw and clipped p-values, bootstrap moments, and calibration rule,
plus \code{variance_audit} describing its numerical variance convention.

\code{lfpm_test} returns named vectors \code{scores} and \code{pvalues},
an \code{AdaL} list of component details, \code{parameter}, and
\code{method_version}. COM's score is the negative of its p-value so larger
scores always indicate stronger evidence against the joint null.
}
\examples{
dat <- simulate_lfpm(1, N = 40, T_size = 30, seed = 17)
lfpm_candidate_grid(ncol(dat$Y))
l_test(dat$Y, dat$F, k = 5, B = 19, seed = 23)
adal_test(dat$Y, dat$F, B = 19, seed = 23)
adal_test(dat$Y, dat$F, B = 19, calibration = "empirical", seed = 23)
lfpm_test(dat$Y, dat$F, B = 19, seed = 23)$pvalues
}

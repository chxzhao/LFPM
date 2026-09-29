# Public interfaces. The frozen engines in core.R, fast.R and ma_engine.R
# retain their original definitions; validation and reproducibility live here.

.lfpm_integer <- function(x, name, minimum = 0L, maximum = .Machine$integer.max) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > maximum) {
    stop(name, " must be one integer between ", minimum, " and ", maximum,
         call. = FALSE)
  }
  as.integer(x)
}

.lfpm_matrix <- function(x, name, allow_empty_columns = FALSE) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || !nrow(x) ||
      (!allow_empty_columns && !ncol(x)) || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix (or numeric data frame)",
         call. = FALSE)
  }
  x
}

.lfpm_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  seed <- .lfpm_integer(seed, "seed")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.lfpm_multipliers <- function(multipliers, T_size, B) {
  if (is.null(multipliers)) {
    return(matrix(sample(c(-1, 1), T_size * B, replace = TRUE), T_size, B))
  }
  if (!is.matrix(multipliers) || !is.numeric(multipliers) || is.complex(multipliers) ||
      !identical(dim(multipliers), c(as.integer(T_size), as.integer(B))) ||
      any(!is.finite(multipliers)) || any(!multipliers %in% c(-1, 1))) {
    stop("multipliers must be a T_size by B numeric matrix containing only -1 and 1",
         call. = FALSE)
  }
  multipliers
}

lfpm_fit <- function(Y, F) {
  Y <- .lfpm_matrix(Y, "Y")
  F <- .lfpm_matrix(F, "F", allow_empty_columns = TRUE)
  if (nrow(Y) != nrow(F)) stop("Y and F must have the same number of rows", call. = FALSE)
  if (nrow(Y) <= ncol(F) + 1L) {
    stop("There must be positive residual degrees of freedom: T_size > ncol(F) + 1",
         call. = FALSE)
  }
  fit <- fit_panel(Y, F)
  if (any(!is.finite(c(fit$alpha, fit$s2, fit$q))) || any(fit$s2 <= 0)) {
    stop("Nonfinite or degenerate fitted panel", call. = FALSE)
  }
  fit$tstat <- sqrt(fit$q) * fit$alpha / sqrt(fit$s2)
  if (any(!is.finite(fit$tstat))) stop("Nonfinite studentized intercept", call. = FALSE)
  class(fit) <- "lfpm_fit"
  fit
}

simulate_lfpm <- function(model, N, T_size, distribution = "Normal", m = 0L,
                          signal_energy = 10, seed = NULL) {
  model <- .lfpm_integer(model, "model", 1L, 4L)
  N <- .lfpm_integer(N, "N", 3L)
  T_size <- .lfpm_integer(T_size, "T_size", 5L)
  m <- .lfpm_integer(m, "m", 0L, N)
  if (!is.character(distribution) || length(distribution) != 1L ||
      is.na(distribution) || !distribution %in% c("Normal", "t5", "MixNorm", "Chisq")) {
    stop("distribution must be one of Normal, t5, MixNorm or Chisq", call. = FALSE)
  }
  if (!is.numeric(signal_energy) || is.complex(signal_energy) || length(signal_energy) != 1L ||
      !is.finite(signal_energy) || signal_energy < 0) {
    stop("signal_energy must be one finite nonnegative number", call. = FALSE)
  }
  .lfpm_seed(seed, make_data(model, N, T_size, distribution, m, signal_energy))
}

lfpm_candidate_grid <- function(N) {
  N <- .lfpm_integer(N, "N", 40L)
  ma_candidate_k(N)
}

l_test <- function(Y, F, k, B = 999L, seed = NULL, multipliers = NULL) {
  data_name <- paste(deparse(substitute(Y)), "and", deparse(substitute(F)))
  fit <- lfpm_fit(Y, F)
  k <- .lfpm_integer(k, "k", 1L, fit$N)
  B <- .lfpm_integer(B, "B", 2L)
  .lfpm_seed(seed, {
    weights <- .lfpm_multipliers(multipliers, fit$T_size, B)
    # The frozen row_top_sums assumes at least two candidate columns. Duplicate
    # this single rank only to preserve its matrix shape, including when N=1.
    boot <- ma_bootstrap_L(fit, K = c(k, k), B = B,
                           multipliers = weights)[, 1L]
    observed <- top_sums(fit$alpha^2 * fit$q / fit$s2, k)
    p_value <- (1 + sum(boot >= observed)) / (B + 1)
    structure(list(statistic = setNames(observed, paste0("L", k)),
                   parameter = c(k = k, B = B, N = fit$N, T_size = fit$T_size),
                   p.value = p_value, alternative = "greater",
                   method = "Studentized LFPM L-test with refitted residual wild bootstrap",
                   data.name = data_name, bootstrap.statistics = boot,
                   calibration = "empirical", bootstrap_variance = "refit"),
              class = "htest")
  })
}

.lfpm_empirical_bootstrap_L <- function(fit, K, B, multipliers) {
  # Preserve the audited safeguard in original_stream_bootstrap() from
  # scripts/adal_fullgrid_empirical_oldstream_20260929.py. This historical
  # empirical rule is deliberately separate from the Ma simulation engine.
  alpha_boot <- crossprod(multipliers * fit$P, fit$U)
  projected_ss <- matrix(0, B, fit$N)
  for (j in seq_len(ncol(fit$Q))) {
    projected_ss <- projected_ss + crossprod(multipliers * fit$Q[, j], fit$U)^2
  }
  variance <- (matrix(colSums(fit$U^2), B, fit$N, byrow = TRUE) - projected_ss) / fit$df
  if (any(!is.finite(variance))) stop("Nonfinite empirical bootstrap variance", call. = FALSE)
  audit <- list(variance_floor = 1e-20,
    minimum_unclipped_refitted_variance = min(variance),
    legacy_variance_safeguard_count = sum(variance <= 1e-20))
  variance <- pmax(variance, 1e-20)
  L <- row_top_sums(alpha_boot^2 * fit$q / variance, K)
  if (any(!is.finite(L))) stop("Nonfinite empirical bootstrap statistic", call. = FALSE)
  list(L = L, audit = audit)
}

.lfpm_adal <- function(fit, B, calibration, multipliers) {
  K <- lfpm_candidate_grid(fit$N)
  if (calibration == "simulation") {
    prepared <- ma_prepare(fit, B = B, multipliers = multipliers)
    audited <- ma_score_alphas(fit$alpha, prepared, audit = TRUE)
    observed <- as.numeric(audited$Lobs)
    p_component <- as.numeric(audited$p_L)
    epsilon <- MA_P_ENDPOINT_EPS
    statistic <- unname(audited$values[1L, "score_AdaL"])
    p_value <- unname(audited$values[1L, "p_AdaL"])
    boot_mean <- audited$boot_mean
    boot_sd <- audited$boot_sd
    rule <- c("add-one empirical upper tail", rep("bootstrap-moment normal upper tail", length(K) - 1L))
    variance_audit <- list(variance_floor = NULL)
  } else {
    empirical_boot <- .lfpm_empirical_bootstrap_L(fit, K, B, multipliers)
    boot <- empirical_boot$L
    variance_audit <- empirical_boot$audit
    observed <- top_sums(fit$alpha^2 * fit$q / fit$s2, K)
    p_component <- (1 + colSums(sweep(boot, 2L, observed, ">="))) / (B + 1)
    epsilon <- 1 / (fit$N + fit$T_size)
    clipped <- pmin(pmax(p_component, epsilon), 1 - epsilon)
    statistic <- mean(tan(pi * (.5 - clipped)))
    p_value <- pcauchy(statistic, lower.tail = FALSE)
    boot_mean <- colMeans(boot)
    boot_sd <- apply(boot, 2L, sd)
    rule <- rep("add-one empirical upper tail", length(K))
  }
  list(statistic = statistic, p.value = p_value,
       components = data.frame(k = K, statistic = observed,
         p.value = p_component,
         p.combined = pmin(pmax(p_component, epsilon), 1 - epsilon),
         bootstrap.mean = as.numeric(boot_mean), bootstrap.sd = as.numeric(boot_sd),
         rule = rule, row.names = NULL, stringsAsFactors = FALSE),
       endpoint_epsilon = epsilon, calibration = calibration,
       bootstrap_variance = "refit", variance_audit = variance_audit)
}

adal_test <- function(Y, F, B = 999L, calibration = c("simulation", "empirical"),
                      seed = NULL, multipliers = NULL) {
  data_name <- paste(deparse(substitute(Y)), "and", deparse(substitute(F)))
  fit <- lfpm_fit(Y, F)
  lfpm_candidate_grid(fit$N)
  B <- .lfpm_integer(B, "B", 2L)
  calibration <- match.arg(calibration)
  answer <- .lfpm_seed(seed, {
    weights <- .lfpm_multipliers(multipliers, fit$T_size, B)
    .lfpm_adal(fit, B, calibration, weights)
  })
  answer$statistic <- c(AdaL = answer$statistic)
  answer$parameter <- c(B = B, N = fit$N, T_size = fit$T_size)
  answer$alternative <- "greater"
  answer$method <- paste("Adaptive LFPM L-test with", calibration, "component calibration")
  answer$data.name <- data_name
  class(answer) <- "htest"
  answer
}

# Four benchmark formulas and operation order match the frozen test_panel.
# Passing a single multiplier matrix avoids changing the comparator stream
# when the Ma candidate grid differs from the original four-point grid.
.lfpm_comparators <- function(fit, B, multipliers) {
  N <- fit$N
  T_size <- fit$T_size
  v <- fit$df
  t2 <- fit$alpha^2 * fit$q / fit$s2
  theta <- qnorm(1 - .1 / (N - 1) / 2)^2
  R <- cor(fit$U)
  off_diag <- row(R) != col(R)
  rho2 <- sum(R[off_diag & R^2 * v > theta]^2) / (N * (N - 1))
  denom <- (v / (v - 2)) * sqrt(2 * (v - 1) / (v - 4) * (1 + N * rho2))
  PY <- (sum(t2) - N * v / (v - 2)) / sqrt(N) / denom
  MAX <- (max(t2) - 2 * log(N) + log(log(N))) / (1 + 1 / log(T_size * sqrt(N)))
  p_py <- pnorm(PY, lower.tail = FALSE)
  p_max <- -expm1(-exp(-MAX / 2) / sqrt(pi))
  p_com <- min(1, 2 * min(p_py, p_max))
  K <- unique(as.integer(c(1L, 5L, floor(N * .25), floor(N * .5))))
  Gamma <- precision_estimate(fit, theta)
  boot <- bootstrap_components(fit, Gamma, K, B, "refit", multipliers)
  Zobs <- as.vector(Gamma %*% (fit$alpha / sqrt(fit$s2)))
  Gobs <- top_sums(fit$q * Zobs^2 / diag(Gamma), K)
  EG <- colMeans(boot$G)
  SDG <- apply(boot$G, 2L, sd) + 1e-8
  maxk <- max((Gobs - EG) / SDG)
  maxk_boot <- apply(sweep(sweep(boot$G, 2L, EG, "-"), 2L, SDG, "/"), 1L, max)
  p_maxk <- (1 + sum(maxk_boot >= maxk)) / (B + 1)
  list(scores = setNames(c(PY, MAX, -p_com, maxk), METHODS[1:4]),
       pvalues = setNames(c(p_py, p_max, p_com, p_maxk), METHODS[1:4]))
}

lfpm_test <- function(Y, F, B = 999L, seed = NULL) {
  fit <- lfpm_fit(Y, F)
  lfpm_candidate_grid(fit$N)
  if (fit$df <= 4L) stop("The PY benchmark requires residual degrees of freedom > 4", call. = FALSE)
  B <- .lfpm_integer(B, "B", 2L)
  .lfpm_seed(seed, {
    weights <- .lfpm_multipliers(NULL, fit$T_size, B)
    answer <- .lfpm_comparators(fit, B, weights)
    adal <- .lfpm_adal(fit, B, "simulation", weights)
    answer$scores <- c(answer$scores, AdaL = adal$statistic)
    answer$pvalues <- c(answer$pvalues, AdaL = adal$p.value)
    answer$AdaL <- adal
    answer$parameter <- c(B = B, N = fit$N, T_size = fit$T_size)
    answer$method_version <- MA_ENGINE_VERSION
    if (any(!is.finite(c(answer$scores, answer$pvalues)))) {
      stop("Nonfinite test result", call. = FALSE)
    }
    answer
  })
}

# Ma (2024)-definition AdaL within the unchanged LFPM residual-wild-bootstrap
# experiment. Source the preserved lfpm_core.R before this file. No execution
# takes place at source time. All random draws follow the old phase/seed order.
MA_ENGINE_VERSION <- "ma2024-adal-lfpm-20260929-v1"
MA_P_ENDPOINT_EPS <- .Machine$double.eps
MA_BOOTSTRAP_VARIANCE <- "refit"
MA_COMPONENT_P_RULE <- "k5:add-one >= empirical tail; proportional k: normal right-tail with bootstrap mean and sample sd"

ma_candidate_k <- function(N) {
  if (length(N) != 1L || !is.finite(N) || N != floor(N) || N < 40L) {
    stop("Ma candidate grid requires one integer N >= 40")
  }
  K <- as.integer(floor(log(N / 20) / log(2)))
  candidates <- as.integer(c(5L, ceiling(N / 2^seq_len(K))))
  if (anyDuplicated(candidates) || length(candidates) != K + 1L ||
      any(candidates < 1L | candidates > N)) stop("Invalid Ma candidate set")
  candidates
}

ma_bootstrap_L <- function(fit, K = ma_candidate_k(fit$N), B = 999L,
                           multipliers = NULL) {
  if (length(B) != 1L || !is.finite(B) || B < 2L || B != floor(B)) stop("Invalid B")
  B <- as.integer(B)
  T_size <- as.integer(fit$T_size)
  if (is.null(multipliers)) {
    # Same month-level Rademacher draws, column ordering and RNG consumption as
    # bootstrap_components() in the preserved five-method LFPM engine.
    multipliers <- matrix(sample(c(-1, 1), T_size * B, replace = TRUE), T_size, B)
  }
  if (!identical(dim(multipliers), c(T_size, B)) ||
      any(!multipliers %in% c(-1, 1))) stop("Invalid bootstrap multipliers")
  alpha_boot <- crossprod(multipliers * fit$P, fit$U)
  rss <- matrix(colSums(fit$U^2), B, fit$N, byrow = TRUE)
  # Full OLS-refit RSS on (1,F), exactly matching the old refit engine.
  for (j in seq_len(ncol(fit$Q))) {
    rss <- rss - crossprod(multipliers * fit$Q[, j], fit$U)^2
  }
  if (any(!is.finite(rss)) || any(rss <= 0)) stop("Invalid bootstrap RSS")
  L <- row_top_sums(alpha_boot^2 * fit$q / (rss / fit$df), K)
  colnames(L) <- paste0("k", K)
  if (any(!is.finite(L))) stop("Nonfinite bootstrap L statistic")
  L
}

ma_prepare <- function(fit, B = 999L, multipliers = NULL, keep_bootstrap = FALSE) {
  K <- ma_candidate_k(fit$N)
  boot_L <- ma_bootstrap_L(fit, K, B, multipliers)
  boot_mean <- colMeans(boot_L)
  boot_sd <- apply(boot_L, 2L, sd)
  # No statistical sd regularization or old .001/.999 truncation is used.
  if (any(!is.finite(boot_sd)) || any(boot_sd <= 0)) stop("Degenerate bootstrap component")
  result <- list(fit = fit, B = as.integer(B), K = K,
                 boot_mean = boot_mean, boot_sd = boot_sd,
                 fixed5_sorted = sort(boot_L[, 1L]))
  if (keep_bootstrap) result$boot_L <- boot_L
  result
}

ma_combine_p <- function(p_L) {
  if (is.null(dim(p_L))) p_L <- matrix(p_L, nrow = 1L)
  if (!is.matrix(p_L) || !ncol(p_L) || any(!is.finite(p_L)) ||
      any(p_L < 0 | p_L > 1)) stop("Invalid component p values")
  p_numeric <- pmin(pmax(p_L, MA_P_ENDPOINT_EPS), 1 - MA_P_ENDPOINT_EPS)
  # There are K+1 distinct candidates, hence rowMeans is precisely 1/(K+1).
  statistic <- rowMeans(tan(pi * (.5 - p_numeric)))
  # Stable numerical evaluation of 1/2 - atan(T_C)/pi, including large tails.
  combined_p <- pcauchy(statistic, lower.tail = FALSE)
  out <- cbind(score_AdaL = statistic, p_AdaL = combined_p)
  if (any(!is.finite(out))) stop("Nonfinite AdaL statistic")
  out
}

ma_score_alphas <- function(alpha, prepared, audit = FALSE) {
  if (is.null(dim(alpha))) alpha <- matrix(alpha, nrow = 1L)
  fit <- prepared$fit
  if (!is.matrix(alpha) || ncol(alpha) != fit$N ||
      nrow(alpha) < 1L || any(!is.finite(alpha))) stop("Invalid alpha matrix")
  t2 <- sweep(alpha^2 * fit$q, 2L, fit$s2, "/")
  Lobs <- row_top_sums(t2, prepared$K)
  colnames(Lobs) <- paste0("k", prepared$K)
  p_L <- matrix(NA_real_, nrow(alpha), length(prepared$K), dimnames = dimnames(Lobs))
  p_L[, 1L] <- (1 + prepared$B - findInterval(Lobs[, 1L], prepared$fixed5_sorted,
                                            left.open = TRUE)) / (prepared$B + 1)
  normal_idx <- seq.int(2L, length(prepared$K))
  z <- sweep(sweep(Lobs[, normal_idx, drop = FALSE], 2L,
                   prepared$boot_mean[normal_idx], "-"), 2L,
             prepared$boot_sd[normal_idx], "/")
  p_L[, normal_idx] <- pnorm(z, lower.tail = FALSE)
  values <- ma_combine_p(p_L)
  if (!audit) return(values)
  list(values = values, K = prepared$K, Lobs = Lobs, p_L = p_L,
       boot_mean = prepared$boot_mean, boot_sd = prepared$boot_sd,
       endpoint_epsilon = MA_P_ENDPOINT_EPS)
}

ma_validate_cfg <- function(cfg) {
  if (!is.null(cfg$bootstrap_variance) && cfg$bootstrap_variance != MA_BOOTSTRAP_VARIANCE) {
    stop("Ma engine requires full OLS-refit bootstrap variance")
  }
  if (is.null(cfg$B) || is.null(cfg$seed) || is.null(cfg$signal_energy)) {
    stop("cfg must provide B, seed, signal_energy")
  }
  invisible(TRUE)
}

ma_one_repetition <- function(index, setting, cfg, phase, audit = FALSE) {
  ma_validate_cfg(cfg)
  if (!phase %in% c("calibration", "size")) stop("Use calibration or size phase")
  if (!identical(as.integer(setting$m), 0L)) stop("Calibration/size must have m=0")
  key <- paste(phase, setting$model, setting$distribution, setting$T_size,
               setting$N, setting$m, sep = "_")
  set.seed(seed_for(cfg$seed, key, index))
  dat <- make_data(setting$model, setting$N, setting$T_size, setting$distribution,
                   setting$m, cfg$signal_energy)
  fit <- fit_panel(dat$Y, dat$F)
  prepared <- ma_prepare(fit, cfg$B)
  ans <- ma_score_alphas(fit$alpha, prepared, audit)
  if (audit) {
    ans$values <- ans$values[1L, ]
    return(ans)
  }
  ans[1L, ]
}

ma_one_power_curve <- function(index, setting, cfg, m_grid, audit = FALSE) {
  ma_validate_cfg(cfg)
  N <- setting$N
  if (!length(m_grid) || any(!is.finite(m_grid)) || any(m_grid != floor(m_grid)) ||
      any(m_grid < 0L | m_grid > N) || anyDuplicated(m_grid)) stop("Invalid power m grid")
  key <- paste("power_curve", setting$model, setting$distribution, setting$T_size,
               N, sep = "_")
  set.seed(seed_for(cfg$seed, key, index))
  dat <- make_data(setting$model, N, setting$T_size, setting$distribution,
                   0L, cfg$signal_energy)
  # Draw support BEFORE bootstrap just as the previous CRN engine did.
  support <- sample.int(N)
  fit <- fit_panel(dat$Y, dat$F)
  prepared <- ma_prepare(fit, cfg$B)
  alpha <- matrix(fit$alpha, length(m_grid), N, byrow = TRUE)
  for (j in seq_along(m_grid)) {
    m <- m_grid[j]
    if (m > 0L) {
      selected <- support[seq_len(m)]
      alpha[j, selected] <- alpha[j, selected] +
        sqrt(cfg$signal_energy * log(N) / (setting$T_size * m))
    }
  }
  ans <- ma_score_alphas(alpha, prepared, audit)
  if (audit) {
    rownames(ans$values) <- rownames(ans$Lobs) <- rownames(ans$p_L) <- as.character(m_grid)
    return(ans)
  }
  rownames(ans) <- as.character(m_grid)
  ans
}

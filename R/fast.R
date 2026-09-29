# Exact cross-m reuse of the original five-method statistics.
# Source lfpm_core.R first. Each m retains its original marginal DGP; a curve
# shares factors/errors, a nested uniform support and bootstrap multipliers.
# Calibration and independent size evaluation must keep their original phases.
FAST_ENGINE_VERSION <- "five-methods-fast-curve-20260928-v1"

prepare_tests <- function(fit, B = 999L, bootstrap_variance = "refit", adal_include_l1 = FALSE) {
  N <- fit$N
  v <- fit$df
  theta <- qnorm(1 - .1 / (N - 1) / 2)^2
  R <- cor(fit$U)
  off_diag <- row(R) != col(R)
  rho2 <- sum(R[off_diag & R^2 * v > theta]^2) / (N * (N - 1))
  py_denom <- (v / (v - 2)) * sqrt(2 * (v - 1) / (v - 4) * (1 + N * rho2))
  K <- unique(as.integer(c(1L, 5L, floor(N * .25), floor(N * .5))))
  Gamma <- precision_estimate(fit, theta)
  boot <- bootstrap_components(fit, Gamma, K, B, bootstrap_variance)
  EG <- colMeans(boot$G)
  SDG <- apply(boot$G, 2L, sd) + 1e-8
  maxk_boot <- apply(sweep(sweep(boot$G, 2L, EG, "-"), 2L, SDG, "/"), 1L, max)
  L_sorted <- vapply(seq_along(K), function(j) sort(boot$L[, j]), numeric(B))
  list(fit = fit, B = B, K = K, Gamma = Gamma, EG = EG, SDG = SDG,
       py_center = N * v / (v - 2), py_denom = py_denom,
       max_center = 2 * log(N) - log(log(N)),
       max_scale = 1 + 1 / log(fit$T_size * sqrt(N)),
       adal_idx = if (adal_include_l1) seq_along(K) else which(K != 1L),
       L_sorted = L_sorted, maxk_sorted = sort(maxk_boot))
}

right_tail_from_sorted <- function(observed, null_sorted) {
  # left.open=TRUE counts null values strictly below observed, including ties
  # exactly as (1 + sum(null >= observed)) / (B + 1) in the original engine.
  (1 + length(null_sorted) - findInterval(observed, null_sorted, left.open = TRUE)) /
    (length(null_sorted) + 1)
}

score_alphas <- function(alpha, prepared) {
  # Rows are observed intercept vectors; columns are assets.
  if (is.null(dim(alpha))) alpha <- matrix(alpha, nrow = 1L)
  fit <- prepared$fit
  if (!is.matrix(alpha) || ncol(alpha) != fit$N || any(!is.finite(alpha))) stop("Invalid alpha matrix")
  t2 <- sweep(alpha^2 * fit$q, 2L, fit$s2, "/")
  PY <- (rowSums(t2) - prepared$py_center) / sqrt(fit$N) / prepared$py_denom
  MAX <- (apply(t2, 1L, max) - prepared$max_center) / prepared$max_scale
  p_py <- pnorm(PY, lower.tail = FALSE)
  p_max <- -expm1(-exp(-MAX / 2) / sqrt(pi))
  p_com <- pmin(1, 2 * pmin(p_py, p_max))

  Lobs <- row_top_sums(t2, prepared$K)
  p_L <- vapply(seq_along(prepared$K), function(j) {
    right_tail_from_sorted(Lobs[, j], prepared$L_sorted[, j])
  }, numeric(nrow(alpha)))
  p_L <- matrix(p_L, nrow(alpha), length(prepared$K))
  p_s <- pmax(pmin(p_L[, prepared$adal_idx, drop = FALSE], .999), .001)
  AdaL <- rowMeans(tan(pi * (.5 - p_s)))
  p_adal <- pcauchy(AdaL, lower.tail = FALSE)

  # Transpose preserves the original observed Gamma %*% alpha expression,
  # including any insignificant floating-point asymmetry in its inverse.
  Z <- sweep(alpha, 2L, sqrt(fit$s2), "/") %*% t(prepared$Gamma)
  Gobs <- row_top_sums(sweep(Z^2 * fit$q, 2L, diag(prepared$Gamma), "/"), prepared$K)
  maxk <- apply(sweep(sweep(Gobs, 2L, prepared$EG, "-"), 2L, prepared$SDG, "/"), 1L, max)
  p_maxk <- right_tail_from_sorted(maxk, prepared$maxk_sorted)
  out <- cbind(PY, MAX, -p_com, maxk, AdaL, p_py, p_max, p_com, p_maxk, p_adal)
  colnames(out) <- c(paste0("score_", METHODS), paste0("p_", METHODS))
  if (any(!is.finite(out))) stop("Nonfinite statistic")
  out
}

score_alpha <- function(alpha, prepared) {
  if (length(alpha) != prepared$fit$N || !is.null(dim(alpha))) stop("Expected one intercept vector")
  result <- score_alphas(alpha, prepared)[1L, ]
  list(scores = setNames(result[1:5], METHODS), pvalues = setNames(result[6:10], METHODS))
}

one_power_curve <- function(index, setting, cfg, m_grid) {
  N <- setting$N
  if (!length(m_grid) || any(!is.finite(m_grid)) || any(m_grid != floor(m_grid)) ||
      any(m_grid < 0L | m_grid > N) || anyDuplicated(m_grid)) stop("Invalid power m grid")
  key <- paste("power_curve", setting$model, setting$distribution, setting$T_size, N, sep = "_")
  set.seed(seed_for(cfg$seed, key, index))
  dat <- make_data(setting$model, N, setting$T_size, setting$distribution, 0L, cfg$signal_energy)
  support <- sample.int(N)
  fit <- fit_panel(dat$Y, dat$F)
  prepared <- prepare_tests(fit, cfg$B, cfg$bootstrap_variance, cfg$adal_include_l1)
  alpha <- matrix(fit$alpha, length(m_grid), N, byrow = TRUE)
  for (j in seq_along(m_grid)) {
    m <- m_grid[j]
    if (m > 0L) {
      selected <- support[seq_len(m)]
      alpha[j, selected] <- alpha[j, selected] + sqrt(cfg$signal_energy * log(N) / (setting$T_size * m))
    }
  }
  out <- score_alphas(alpha, prepared)
  rownames(out) <- as.character(m_grid)
  out
}

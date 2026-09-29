# Five-method LFPM simulation. Base R only; source without running simulations.
# max-k retains the original AT_XZ implementation, not a claim of exact
# reproduction of a particular published max-k algorithm.
METHODS <- c("PY", "MAX", "COM", "max-k", "AdaL")
ENGINE_VERSION <- "five-methods-20260927-v1"

make_covariance <- function(model, N) {
  if (model == 1L) return(0.7 ^ abs(outer(seq_len(N), seq_len(N), "-")))
  if (model == 2L) {
    d <- runif(N, 1, 2)
    b <- numeric(N)
    s <- floor(N ^ 0.3)
    b[sample.int(N, s)] <- runif(s, 0.7, 0.9)
    R <- diag(N) + tcrossprod(b) - diag(b ^ 2)
    return(R * sqrt(outer(d, d)))
  }
  if (model == 3L) return(solve(0.6 ^ abs(outer(seq_len(N), seq_len(N), "-"))))
  if (model == 4L) {
    W <- matrix(0, N, N)
    for (i in 2:(N - 1L)) W[i, c(i - 1L, i + 1L)] <- 0.5
    W[1L, 2L] <- W[N, N - 1L] <- 1
    A <- solve(diag(N) - 0.5 * W)
    b <- numeric(N)
    b[seq_len(floor(N ^ 0.3))] <- runif(floor(N ^ 0.3), 0.7, 0.9)
    return(tcrossprod(b) + tcrossprod(A))
  }
  stop("Unknown covariance model")
}

make_factors <- function(T_size) {
  pars <- rbind(c(.06, .85, .89, .11, .53),
                c(.19, .74, .62, .19, .19),
                c(.05, .76, .80, .15, .19))
  F <- matrix(0, T_size, 3L)
  for (j in 1:3) {
    f <- numeric(T_size + 50L)
    h <- rep(1, T_size + 50L)
    z <- rnorm(T_size + 50L)
    for (i in seq_len(T_size + 49L)) {
      # Preserve original recursion; z^2 is not h*z^2 (not standard GARCH).
      h[i + 1L] <- pars[j, 2L] * h[i] + pars[j, 3L] + pars[j, 4L] * z[i]^2
      f[i + 1L] <- pars[j, 1L] * f[i] + sqrt(h[i + 1L]) * z[i + 1L] + pars[j, 5L]
    }
    F[, j] <- f[51L:(T_size + 50L)]
  }
  F
}

make_innovations <- function(T_size, N, distribution) {
  n <- T_size * N
  z <- switch(distribution,
    Normal = rnorm(n),
    t5 = rt(n, 5) / sqrt(5 / 3),
    MixNorm = {
      a <- rbinom(n, 1, .1)
      ((1 - a) * rnorm(n) + a * rnorm(n, sd = 3)) / sqrt(1.8)
    },
    Chisq = (rchisq(n, 10) - 10) / sqrt(20),
    stop("Unknown distribution"))
  matrix(z, T_size, N)
}

make_data <- function(model, N, T_size, distribution, m, signal_energy = 10) {
  Sigma <- make_covariance(model, N)
  U <- make_innovations(T_size, N, distribution) %*% chol(Sigma)
  alpha <- numeric(N)
  if (m > 0L) alpha[sample.int(N, m)] <- sqrt(signal_energy * log(N) / (T_size * m))
  F <- make_factors(T_size)
  beta <- matrix(runif(3L * N, -1, 1), 3L, N)
  Y <- sweep(F %*% beta + U, 2L, alpha, "+")
  list(Y = Y, F = F)
}

fit_panel <- function(Y, F) {
  T_size <- nrow(Y)
  design <- cbind(1, F)
  design_qr <- qr(design)
  if (design_qr$rank != ncol(design)) stop("Rank deficient factor design")
  one_resid <- qr.resid(qr(F), rep(1, T_size))
  q <- sum(one_resid)
  residuals <- qr.resid(design_qr, Y)
  alpha <- as.numeric(qr.coef(design_qr, Y)[1L, ])
  df <- T_size - ncol(design)
  s2 <- colSums(residuals^2) / df
  if (q <= 0 || any(s2 <= 0)) stop("Degenerate regression")
  list(alpha = alpha, U = residuals, s2 = s2, df = df, q = q,
       P = one_resid / q, Q = qr.Q(design_qr), T_size = T_size, N = ncol(Y))
}

top_sums <- function(x, K) cumsum(sort(x, decreasing = TRUE))[K]
row_top_sums <- function(x, K) t(vapply(seq_len(nrow(x)), function(i) top_sums(x[i, ], K), numeric(length(K))))

precision_estimate <- function(fit, theta) {
  sut <- sweep(fit$U, 2L, sqrt(fit$s2), "/")
  # Equivalent to cvCovEst::thresholdingEst(dat, gamma), checked locally.
  S <- cov(sut)
  S[abs(S) < 1.5 * sqrt(theta / fit$df)] <- 0
  # Preserve original AT_XZ ridge, including its shift when S is positive definite.
  shift <- abs(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)) + .05
  S <- S + diag(shift, fit$N)
  solve(S)
}

bootstrap_components <- function(fit, Gamma, K, B, variance_mode = "refit", multipliers = NULL) {
  if (!variance_mode %in% c("refit", "legacy")) stop("Invalid bootstrap variance mode")
  T_size <- fit$T_size
  N <- fit$N
  if (is.null(multipliers)) multipliers <- matrix(sample(c(-1, 1), T_size * B, replace = TRUE), T_size, B)
  stopifnot(identical(dim(multipliers), c(as.integer(T_size), as.integer(B))))
  # Each month uses the same multiplier for every asset.
  ab <- crossprod(multipliers * fit$P, fit$U)
  rss <- matrix(colSums(fit$U^2), B, N, byrow = TRUE)
  if (variance_mode == "refit") {
    # RSS of a complete OLS refit on (1,F): ||u*||^2 - ||Q'u*||^2.
    for (j in seq_len(ncol(fit$Q))) rss <- rss - crossprod(multipliers * fit$Q[, j], fit$U)^2
  } else {
    # Exact algebraic form of the old code's sum((u* - alpha*)^2).
    rss <- rss - 2 * ab * crossprod(multipliers, fit$U) + T_size * ab^2
  }
  if (any(rss <= 0)) stop("Nonpositive bootstrap residual sum of squares")
  L <- row_top_sums(ab^2 * fit$q / (rss / fit$df), K)
  Z <- sweep(ab, 2L, sqrt(fit$s2), "/") %*% Gamma
  G <- row_top_sums(sweep(Z^2 * fit$q, 2L, diag(Gamma), "/"), K)
  list(L = L, G = G)
}

test_panel <- function(Y, F, B = 999L, bootstrap_variance = "refit", adal_include_l1 = FALSE) {
  fit <- fit_panel(Y, F)
  N <- fit$N
  T_size <- fit$T_size
  v <- fit$df
  t2 <- fit$alpha^2 * fit$q / fit$s2
  theta <- qnorm(1 - .1 / (N - 1) / 2)^2
  R <- cor(fit$U)
  # Explicitly exclude the diagonal. Equivalent on the original grid, and
  # avoids a negative variance estimate for very short custom time series.
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
  boot <- bootstrap_components(fit, Gamma, K, B, bootstrap_variance)
  Lobs <- top_sums(t2, K)
  p_L <- (1 + colSums(sweep(boot$L, 2L, Lobs, ">="))) / (B + 1)
  adal_idx <- if (adal_include_l1) seq_along(K) else which(K != 1L)
  # Keep the original truncation explicit to preserve the original AdaL definition.
  p_s <- pmax(pmin(p_L[adal_idx], .999), .001)
  AdaL <- mean(tan(pi * (.5 - p_s)))
  p_adal <- pcauchy(AdaL, lower.tail = FALSE)

  Zobs <- as.vector(Gamma %*% (fit$alpha / sqrt(fit$s2)))
  Gobs <- top_sums(fit$q * Zobs^2 / diag(Gamma), K)
  EG <- colMeans(boot$G)
  SDG <- apply(boot$G, 2L, sd) + 1e-8
  maxk <- max((Gobs - EG) / SDG)
  maxk_boot <- apply(sweep(sweep(boot$G, 2L, EG, "-"), 2L, SDG, "/"), 1L, max)
  p_maxk <- (1 + sum(maxk_boot >= maxk)) / (B + 1)
  scores <- setNames(c(PY, MAX, -p_com, maxk, AdaL), METHODS)
  pvalues <- setNames(c(p_py, p_max, p_com, p_maxk, p_adal), METHODS)
  if (any(!is.finite(c(scores, pvalues)))) stop("Nonfinite statistic")
  list(scores = scores, pvalues = pvalues)
}

# A seed depends on the phase, scenario and repetition, not scheduling or resume order.
seed_for <- function(base_seed, key, repetition) {
  modulus <- 2147483646
  h <- as.double(base_seed)
  for (v in utf8ToInt(paste(key, repetition, sep = ":"))) h <- (h * 131 + v) %% modulus
  as.integer(h + 1)
}

one_repetition <- function(index, setting, cfg, phase) {
  key <- paste(phase, setting$model, setting$distribution, setting$T_size, setting$N, setting$m, sep = "_")
  set.seed(seed_for(cfg$seed, key, index))
  dat <- make_data(setting$model, setting$N, setting$T_size, setting$distribution, setting$m, cfg$signal_energy)
  result <- test_panel(dat$Y, dat$F, cfg$B, cfg$bootstrap_variance, cfg$adal_include_l1)
  c(setNames(result$scores, paste0("score_", METHODS)), setNames(result$pvalues, paste0("p_", METHODS)))
}

calibrated_pvalues <- function(scores, calibration_scores) {
  # Independent-null rank p values. >= makes ties conservative.
  out <- vapply(seq_len(ncol(scores)), function(j) {
    vapply(scores[, j], function(s) (1 + sum(calibration_scores[, j] >= s)) /
             (nrow(calibration_scores) + 1), numeric(1))
  }, numeric(nrow(scores)))
  labels <- if (ncol(scores) == length(METHODS)) METHODS else colnames(scores)
  matrix(out, nrow(scores), ncol(scores), dimnames = list(NULL, labels))
}

binomial_summary <- function(rejections) {
  n <- length(rejections)
  p <- mean(rejections)
  z <- qnorm(.975)
  den <- 1 + z^2 / n
  mid <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(rate = p, mcse = sqrt(p * (1 - p) / n), lower95 = mid - half, upper95 = mid + half, repetitions = n)
}

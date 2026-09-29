library(LFPM)

# Optional integration check: base R package checks do not install a Python stack.
# Set LFPM_PYTHON to a prepared Python executable to exercise this backend.
empirical_python <- Sys.getenv("LFPM_PYTHON", unset = "")
if (nzchar(empirical_python)) {
  set.seed(67)
  empirical_f <- matrix(rnorm(60 * 3), 60, 3)
  empirical_y <- empirical_f %*% matrix(rnorm(3 * 40), 3, 40) +
    matrix(rnorm(60 * 40), 60, 40)
  empirical_out <- tempfile("lfpm-empirical-test-")
  empirical_first <- lfpm_empirical_test(empirical_y, empirical_f, B = 31,
    seed = 762, python = empirical_python, output_dir = empirical_out)
  empirical_second <- lfpm_empirical_test(empirical_y, empirical_f, B = 31,
    seed = 762, python = empirical_python)
  stopifnot(is.finite(empirical_first$p.value), empirical_first$p.value >= 0,
    empirical_first$p.value <= 1,
    identical(empirical_first$components, empirical_second$components),
    all(empirical_first$components$p_value > 0),
    all(empirical_first$components$p_value <= 1),
    all(c("AdaL", "PY", "MAX", "COM", "T1") %in% empirical_first$results$method),
    abs(sum(empirical_first$components$weight) - 1) < 1e-14,
    nrow(empirical_first$alphas) == 40,
    file.exists(empirical_first$metadata))
  empirical_bad <- try(lfpm_empirical_test(empirical_y,
    cbind(1, empirical_f), B = 3, python = empirical_python), silent = TRUE)
  stopifnot(inherits(empirical_bad, "try-error"))
  unlink(c(empirical_out, empirical_second$output_dir), recursive = TRUE)
}

library(LFPM)

# Deliberately synthetic CSVs only: no financial panels or original outputs.
report_fun <- getFromNamespace("report_empirical", "LFPM")
out <- tempfile("synthetic-empirical-report-")
dir.create(out)
write_synthetic <- function(key, shift = 0, maxk = TRUE) {
  path <- file.path(out, key); dir.create(path)
  N <- 80L
  bp <- rep(c(-100, -50, 0, 25, 50, 100, 150, 200), 10) + shift
  p <- c(.0001, .0005, .0015, .004, rep(.8, N - 4))
  alphas <- data.frame(asset = paste0("synthetic-", seq_len(N)),
    alpha_monthly = bp / 10000, alpha_bp_monthly = bp,
    standard_error_bp = rep(25, N), standardized_alpha = bp / 25,
    individual_p = p, BH_q = p.adjust(p, "BH"))
  methods <- c("AdaL", "T1", "PY", "MAX", "COM", if (maxk) "max-k (K=30)")
  results <- data.frame(setting = key, method = methods,
    statistic = seq_along(methods), p_value = seq_along(methods) / 100,
    B = ifelse(methods %in% c("PY", "MAX", "COM"), NA_integer_,
               ifelse(startsWith(methods, "max-k"), 49L, 99L)), N = N, T = 60L)
  results$statistic[results$method == "COM"] <- NA_real_
  components <- data.frame(k = c(5, 20, 40), statistic = c(50, 100, 150),
    exceedances = c(1, 2, 3), p_value = c(.02, .03, .04),
    clipped_p = c(.02, .03, .04), weight = rep(1/3, 3))
  write.csv(results, file.path(path, "results.csv"), row.names = FALSE)
  write.csv(components, file.path(path, "components.csv"), row.names = FALSE)
  write.csv(alphas, file.path(path, "alphas.csv"), row.names = FALSE)
  alphas
}
china <- write_synthetic("china", maxk = FALSE)
partial <- report_fun(out, plots = TRUE)
stopifnot(!partial$main_complete,
          identical(partial$coverage$status, c("reported", "not_computed", "not_computed")),
          partial$coverage$missing_joint_methods[1] == "max-k (K=30)",
          nrow(partial$joint_tests) == 4L,
          nrow(partial$selected_components) == 4L,
          identical(partial$selected_components$k, c(1L, 5L, 20L, 40L)),
          isTRUE(all.equal(partial$selected_components$p_value, c(.02, .02, .03, .04))),
          file.exists(partial$figure), file.info(partial$figure)$size > 1000)
stopifnot(all(is.na(partial$joint_tests$B[partial$joint_tests$method %in% c("PY", "MAX", "COM")])),
          partial$joint_tests$B[partial$joint_tests$method == "AdaL"] == 99L)
# The same allowance must not hide missing/invalid budgets on bootstrap rows.
result_path <- file.path(out, "china", "results.csv")
saved_result <- read.csv(result_path)
for (invalid in c(NA_real_, 0, 1.5, Inf)) {
  bad_result <- saved_result
  bad_result$B[bad_result$method == "AdaL"] <- invalid
  write.csv(bad_result, result_path, row.names = FALSE)
  stopifnot(inherits(try(report_fun(out, FALSE), silent = TRUE), "try-error"))
}
write.csv(saved_result, result_path, row.names = FALSE)
z <- partial$alpha_magnitudes
bp <- china$alpha_bp_monthly
stopifnot(z$mean_bp == mean(bp), z$sd_bp == sd(bp), z$median_bp == median(bp),
          z$q05_bp == unname(quantile(bp, .05)), z$q95_bp == unname(quantile(bp, .95)),
          z$mean_abs_bp == mean(abs(bp)), z$max_abs_bp == max(abs(bp)),
          z$share_abs_ge_50bp == mean(abs(bp) >= 50),
          z$BH_005 == 3L, z$BH_010 == 4L)
us <- write_synthetic("us_ff3", 12)
ff5 <- write_synthetic("us_ff5", -5)
write.csv(data.frame(setting = "china", static_complex_gap_ratio = 2.5,
  loading_average_p = .2, loading_maximum_p = .3, loading_bonferroni_p = .4,
  loading_B = 99L, loading_seed = 100L, loading_average_break = "synthetic-date",
  gap_ratio_is_descriptive = TRUE), file.path(out, "china", "model_diagnostics.csv"), row.names = FALSE)
write.csv(data.frame(setting = "china", diagnostic = c("BP6", "BP12", "LB12", "ARCH6"),
  raw_count = c(4, 5, 6, 7), BH_count = c(1, 2, 3, 4)),
  file.path(out, "china", "diagnostic_summary.csv"), row.names = FALSE)
full <- report_fun(out, plots = TRUE)
stopifnot(full$main_complete, nrow(full$alpha_magnitudes) == 3L,
          nrow(full$joint_tests) == 14L, nrow(full$source_files) == 11L,
          full$alpha_magnitudes$mean_bp[2L] == mean(us$alpha_bp_monthly),
          full$alpha_magnitudes$role[3L] == "factor-set robustness",
          file.exists(full$figure), file.info(full$figure)$size > 1000)
stopifnot(full$coverage$model_diagnostics_available[1L],
          !full$coverage$model_diagnostics_available[2L],
          full$model_diagnostics$static_complex_gap_ratio == 2.5,
          full$diagnostic_table$value[full$diagnostic_table$diagnostic == "Loading-stability p-value"] == .4,
          full$diagnostic_table$value[full$diagnostic_table$diagnostic == "ARCH6 BH discoveries"] == 4L)
head <- readBin(full$figure, "raw", n = 4)
stopifnot(identical(rawToChar(head), "%PDF"))
without_plot <- report_fun(out, plots = FALSE)
stopifnot(is.null(without_plot$figure), any(grepl("plots=FALSE", without_plot$manuscript_mapping$artifact, fixed = TRUE)))

# Missing component orders are declared, never replaced with paper values.
component_path <- file.path(out, "us_ff5", "components.csv")
component <- read.csv(component_path)
write.csv(component[component$k != 20, ], component_path, row.names = FALSE)
missing <- report_fun(out, plots = FALSE)
stopifnot(grepl("T_quarter", missing$coverage$missing_selected_orders[3L]),
          nrow(missing$selected_components[missing$selected_components$setting == "us_ff5", ]) == 3L)

# Partial writes are reported but not treated as complete settings.
unlink(file.path(out, "us_ff5", "alphas.csv"))
incomplete <- report_fun(out, plots = FALSE)
stopifnot(incomplete$coverage$status[3L] == "incomplete",
          nrow(incomplete$alpha_magnitudes) == 2L)

# Mismatched units and forged/inconsistent BH values are refused.
alpha_path <- file.path(out, "china", "alphas.csv")
bad <- china; bad$alpha_monthly <- bad$alpha_monthly * 100
write.csv(bad, alpha_path, row.names = FALSE)
stopifnot(inherits(try(report_fun(out, FALSE), silent = TRUE), "try-error"))
bad <- china; bad$BH_q[1L] <- .99
write.csv(bad, alpha_path, row.names = FALSE)
stopifnot(inherits(try(report_fun(out, FALSE), silent = TRUE), "try-error"))
unlink(out, recursive = TRUE)

# Opt-in end-to-end integration uses actual installed Python backend CSVs,
# rather than assuming a synthetic schema. CI without optional Python packages
# can omit LFPM_PYTHON; local verification sets it to the configured Python.
python <- Sys.getenv("LFPM_PYTHON", unset = "")
if (nzchar(python)) {
  backend_out <- tempfile("actual-backend-synthetic-report-")
  dir.create(backend_out)
  set.seed(4103)
  F <- matrix(rnorm(28 * 3), 28, 3)
  Y <- F %*% matrix(rnorm(3 * 80, sd = .2), 3, 80) +
    matrix(rnorm(28 * 80), 28, 80)
  fit <- lfpm_empirical_test(Y, F, B = 19L, seed = 730L, python = python,
    include_maxk = FALSE, output_dir = file.path(backend_out, "china"))
  # Only relabel the generic matrix input's setting for the report interface.
  actual <- read.csv(file.path(backend_out, "china", "results.csv"))
  stopifnot(all(is.na(actual$B[actual$method %in% c("PY", "MAX", "COM")])),
            all(actual$B[actual$method %in% c("AdaL", "T1")] == 19L))
  actual$setting <- "china"
  write.csv(actual, file.path(backend_out, "china", "results.csv"), row.names = FALSE)
  actual_report <- report_fun(backend_out, plots = TRUE)
  stopifnot(nrow(actual_report$joint_tests) == 4L,
            nrow(actual_report$selected_components) == 4L,
            actual_report$alpha_magnitudes$N == 80L,
            isTRUE(all.equal(actual_report$joint_tests$p_value[actual_report$joint_tests$method == "AdaL"], fit$p.value)),
            isTRUE(all.equal(actual_report$alpha_magnitudes$mean_bp, mean(fit$alphas$alpha_bp_monthly))),
            file.exists(actual_report$figure), file.info(actual_report$figure)$size > 1000)
  unlink(backend_out, recursive = TRUE)
}

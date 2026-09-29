# Reports derived only from caller-supplied backend output, never paper values.
.empirical_report_read <- function(path, required) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(x) || !all(required %in% names(x)))
    stop("Missing rows or required columns in ", path, call. = FALSE)
  x
}

.empirical_report_numeric <- function(x, columns, path, probability = FALSE) {
  for (nm in columns) {
    if (!is.numeric(x[[nm]]) || any(!is.finite(x[[nm]])) ||
        (probability && any(x[[nm]] < 0 | x[[nm]] > 1)))
      stop("Invalid ", nm, " in ", path, call. = FALSE)
  }
  invisible(TRUE)
}

.empirical_report_csv <- function(x, path) {
  temporary <- tempfile(".report-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(x, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Cannot save report table: ", path)
}

.empirical_report_histograms <- function(panels, path) {
  grDevices::pdf(path, width = 10, height = 7.8, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4.1, 4.1, 2.7, 1),
                oma = c(2.4, 0, 1.9, 0), mgp = c(2.5, .8, 0))
  for (setting in c("china", "us_ff3")) {
    market <- if (setting == "china") "China" else "United States (FF3)"
    d <- panels[[setting]]
    for (standardized in c(FALSE, TRUE)) {
      title <- paste(market, if (standardized) "standardized alphas" else "monthly alphas", sep = ": ")
      if (is.null(d)) {
        graphics::plot.new()
        graphics::title(main = title, cex.main = .95)
        graphics::text(.5, .5, "Not computed: complete backend output unavailable", cex = .8)
        next
      }
      values <- if (standardized) d$standardized_alpha else d$alpha_bp_monthly
      xr <- range(c(values, if (standardized) c(-1.96, 1.96)))
      if (diff(xr) == 0) xr <- xr + c(-.5, .5)
      graphics::hist(values, breaks = 35, col = "#2879A8", border = "white",
        main = title, cex.main = .95, xlim = xr,
        xlab = if (standardized) "Standardized alpha" else "Monthly alpha (basis points)",
        ylab = "Number of stocks")
      if (standardized) graphics::abline(v = c(-1.96, 1.96), lty = 2, lwd = 1.1, col = "#A23B30")
    }
  }
  graphics::mtext("Estimated and standardized alphas: supplied backend results",
                 outer = TRUE, side = 3, line = .4, cex = 1)
  graphics::mtext("Dashed lines: pointwise +/-1.96 normal references, not multiplicity-adjusted thresholds.",
                 outer = TRUE, side = 1, line = .6, cex = .75)
  invisible(path)
}

#' Assemble empirical manuscript tables and the alpha distribution figure
#' @export
report_empirical <- function(output_dir, plots = TRUE) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !dir.exists(output_dir))
    stop("output_dir must be an existing empirical backend output directory", call. = FALSE)
  if (!is.logical(plots) || length(plots) != 1L || is.na(plots))
    stop("plots must be TRUE or FALSE", call. = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  settings <- c("china", "us_ff3", "us_ff5")
  labels <- c(china = "China", us_ff3 = "United States (FF3)", us_ff5 = "United States (FF5)")
  roles <- c(china = "main", us_ff3 = "main", us_ff5 = "factor-set robustness")
  coverage <- joint <- selected <- full <- magnitudes <- panels <- sources <- list()
  model_diagnostics <- residual_diagnostics <- diagnostic_table <- list()
  files <- c("results.csv", "components.csv", "alphas.csv")
  methods <- c("PY", "MAX", "COM", "max-k (K=30)", "AdaL")
  for (key in settings) {
    paths <- file.path(output_dir, key, files)
    present <- file.exists(paths)
    if (!all(present)) {
      coverage[[key]] <- data.frame(setting = key, market = unname(labels[key]),
        role = unname(roles[key]), status = if (any(present)) "incomplete" else "not_computed",
        missing_files = paste(files[!present], collapse = "; "),
        missing_joint_methods = "", missing_selected_orders = "", N = NA_integer_, T = NA_integer_)
      next
    }
    r <- .empirical_report_read(paths[1L], c("setting", "method", "statistic", "p_value", "B", "N", "T"))
    c <- .empirical_report_read(paths[2L], c("k", "statistic", "exceedances", "p_value", "clipped_p", "weight"))
    a <- .empirical_report_read(paths[3L], c("asset", "alpha_monthly", "alpha_bp_monthly",
      "standard_error_bp", "standardized_alpha", "individual_p", "BH_q"))
    if (anyNA(r$setting) || any(r$setting != key) || anyNA(r$method) || anyDuplicated(r$method))
      stop("Invalid setting or duplicate method rows in ", paths[1L])
    .empirical_report_numeric(r, "p_value", paths[1L], TRUE)
    .empirical_report_numeric(r, c("N", "T"), paths[1L])
    if (length(unique(r$N)) != 1L || length(unique(r$T)) != 1L ||
        any(r$N < 1 | r$T < 1) ||
        any(r$N != floor(r$N) | r$T != floor(r$T)))
      stop("Inconsistent N or T in ", paths[1L])
    # The backend records method-specific budgets: PY/MAX/COM do not use
    # bootstrap draws and therefore have blank B fields. AdaL/T1/max-k must
    # still have a finite positive integer budget. Preserve analytical NA.
    if (is.logical(r$B) && all(is.na(r$B))) r$B <- as.numeric(r$B)
    analytical <- r$method %in% c("PY", "MAX", "COM")
    observed_B <- r$B[!is.na(r$B)]
    if (!is.numeric(r$B) || anyNA(r$B[!analytical]) ||
        any(!is.finite(observed_B)) || any(observed_B < 1 | observed_B != floor(observed_B)))
      stop("Invalid or missing bootstrap-method B in ", paths[1L])
    N <- r$N[1L]; T_size <- r$T[1L]
    .empirical_report_numeric(c, c("k", "statistic", "exceedances", "weight"), paths[2L])
    .empirical_report_numeric(c, c("p_value", "clipped_p"), paths[2L], TRUE)
    if (anyDuplicated(c$k) || any(c$k != floor(c$k) | c$k < 1 | c$k > N) ||
        any(c$exceedances < 0 | c$exceedances != floor(c$exceedances)) || any(c$weight < 0))
      stop("Invalid component orders, counts or weights in ", paths[2L])
    .empirical_report_numeric(a, c("alpha_monthly", "alpha_bp_monthly", "standard_error_bp",
                                   "standardized_alpha"), paths[3L])
    .empirical_report_numeric(a, c("individual_p", "BH_q"), paths[3L], TRUE)
    if (nrow(a) != N || anyNA(a$asset) || anyDuplicated(a$asset) || any(a$standard_error_bp <= 0))
      stop("Asset identifiers/counts or standard errors are invalid in ", paths[3L])
    if (!isTRUE(all.equal(a$alpha_bp_monthly, 10000 * a$alpha_monthly,
                         tolerance = 1e-10, check.attributes = FALSE)))
      stop("Monthly decimal and basis-point alphas disagree in ", paths[3L])
    bh <- stats::p.adjust(a$individual_p, method = "BH")
    if (any(abs(bh - a$BH_q) > 1e-10)) stop("Saved BH_q disagrees with individual_p in ", paths[3L])
    meta <- data.frame(setting = key, market = unname(labels[key]), role = unname(roles[key]), N = N, T = T_size)
    existing <- methods[methods %in% r$method]
    if (length(existing)) {
      rows <- r[match(existing, r$method), c("method", "statistic", "p_value", "B"), drop = FALSE]
      joint[[key]] <- cbind(meta[rep(1L, nrow(rows)), , drop = FALSE], rows)
      rownames(joint[[key]]) <- NULL
    }
    orders <- c(T1 = 1L, T5 = 5L, T_quarter = as.integer(ceiling(.25 * N)),
                T_half = as.integer(ceiling(.50 * N)))
    selected_rows <- list(); missing_orders <- character()
    for (order in names(orders)) {
      k <- orders[[order]]
      if (order == "T1") {
        row <- r[r$method == "T1", , drop = FALSE]
        if (!nrow(row)) { missing_orders <- c(missing_orders, "T1"); next }
        z <- data.frame(component = order, k = k, statistic = row$statistic,
                        p_value = row$p_value, calibration = "empirical residual bootstrap")
      } else {
        row <- c[c$k == k, , drop = FALSE]
        if (!nrow(row)) { missing_orders <- c(missing_orders, paste0(order, " (k=", k, ")")); next }
        z <- data.frame(component = order, k = k, statistic = row$statistic,
                        p_value = row$p_value, calibration = "empirical residual bootstrap")
      }
      selected_rows[[order]] <- cbind(meta, z)
    }
    if (length(selected_rows)) selected[[key]] <- do.call(rbind, selected_rows)
    full[[key]] <- cbind(meta[rep(1L, nrow(c)), , drop = FALSE], c[, c("k", "statistic", "exceedances", "p_value", "clipped_p", "weight")])
    bp <- a$alpha_bp_monthly
    magnitudes[[key]] <- cbind(meta, data.frame(mean_bp = mean(bp), sd_bp = stats::sd(bp),
      median_bp = stats::median(bp), q05_bp = unname(stats::quantile(bp, .05, type = 7)),
      q95_bp = unname(stats::quantile(bp, .95, type = 7)), mean_abs_bp = mean(abs(bp)),
      max_abs_bp = max(abs(bp)), share_abs_ge_50bp = mean(abs(bp) >= 50),
      BH_005 = sum(bh <= .05), BH_010 = sum(bh <= .10)))
    coverage[[key]] <- data.frame(setting = key, market = unname(labels[key]), role = unname(roles[key]),
      status = "reported", missing_files = "", missing_joint_methods = paste(setdiff(methods, r$method), collapse = "; "),
      missing_selected_orders = paste(missing_orders, collapse = "; "), N = N, T = T_size)
    panels[[key]] <- a
    sources[[key]] <- data.frame(setting = key, file = file.path(key, files),
      md5 = unname(tools::md5sum(paths)), stringsAsFactors = FALSE)
    diagnostic_rows <- data.frame(setting = key, market = unname(labels[key]),
      diagnostic = c("T", "N"), value = c(T_size, N), interpretation = "sample dimension")
    model_path <- file.path(output_dir, key, "model_diagnostics.csv")
    residual_path <- file.path(output_dir, key, "diagnostic_summary.csv")
    for (optional in c(model_path, residual_path)) {
      if (file.exists(optional)) sources[[key]] <- rbind(sources[[key]],
        data.frame(setting = key, file = file.path(key, basename(optional)),
          md5 = unname(tools::md5sum(optional))))
    }
    if (file.exists(model_path)) {
      d <- .empirical_report_read(model_path, c("setting", "static_complex_gap_ratio",
        "loading_average_p", "loading_maximum_p", "loading_bonferroni_p", "loading_B",
        "loading_seed", "loading_average_break", "gap_ratio_is_descriptive"))
      if (nrow(d) != 1L || anyNA(d$setting) || d$setting != key)
        stop("Invalid setting in ", model_path)
      .empirical_report_numeric(d, c("static_complex_gap_ratio", "loading_B", "loading_seed"), model_path)
      .empirical_report_numeric(d, c("loading_average_p", "loading_maximum_p", "loading_bonferroni_p"), model_path, TRUE)
      model_diagnostics[[key]] <- cbind(meta, d[, setdiff(names(d), "setting"), drop = FALSE])
      diagnostic_rows <- rbind(diagnostic_rows, data.frame(setting = key,
        market = unname(labels[key]), diagnostic = c("Static eigenvalue-gap ratio", "Loading-stability p-value"),
        value = c(d$static_complex_gap_ratio, d$loading_bonferroni_p),
        interpretation = c("descriptive; no 7.50 critical-value comparison", "OR-Bonferroni scan on constant-loading specification")))
    }
    if (file.exists(residual_path)) {
      d <- .empirical_report_read(residual_path, c("setting", "diagnostic", "raw_count", "BH_count"))
      if (anyNA(d$setting) || any(d$setting != key) || anyNA(d$diagnostic) || anyDuplicated(d$diagnostic))
        stop("Invalid setting or diagnostic families in ", residual_path)
      .empirical_report_numeric(d, c("raw_count", "BH_count"), residual_path)
      if (any(d$raw_count < 0 | d$raw_count > N | d$raw_count != floor(d$raw_count) |
              d$BH_count < 0 | d$BH_count > N | d$BH_count != floor(d$BH_count)))
        stop("Invalid discovery counts in ", residual_path)
      residual_diagnostics[[key]] <- cbind(meta[rep(1L, nrow(d)), , drop = FALSE],
        d[, c("diagnostic", "raw_count", "BH_count"), drop = FALSE])
      diagnostic_rows <- rbind(diagnostic_rows, data.frame(setting = key, market = unname(labels[key]),
        diagnostic = paste0(d$diagnostic, " BH discoveries"), value = d$BH_count,
        interpretation = "BH at 5% within diagnostic family; fitted main specification"))
    }
    diagnostic_table[[key]] <- diagnostic_rows
  }
  bind <- function(x) {
    if (!length(x)) return(data.frame())
    z <- do.call(rbind, x); rownames(z) <- NULL; z
  }
  result <- list(coverage = bind(coverage), joint_tests = bind(joint),
    selected_components = bind(selected), full_components = bind(full),
    alpha_magnitudes = bind(magnitudes), source_files = bind(sources),
    model_diagnostics = bind(model_diagnostics), residual_diagnostic_summary = bind(residual_diagnostics),
    diagnostic_table = bind(diagnostic_table))
  result$coverage$model_diagnostics_available <- result$coverage$setting %in% names(model_diagnostics)
  result$coverage$residual_diagnostics_available <- result$coverage$setting %in% names(residual_diagnostics)
  main_complete <- all(result$coverage$status[match(c("china", "us_ff3"), result$coverage$setting)] == "reported")
  result$manuscript_mapping <- data.frame(
    manuscript_label = c("emp:tab:joint Panel A", "emp:tab:joint Panel B", "Full-grid component probabilities in text",
      "emp:tab:magnitudes", "Stock-level evidence paragraph", "emp:fig:alphas", "FF5 robustness paragraph"),
    artifact = c("joint_tests.csv", "selected_components.csv", "full_components.csv", "alpha_magnitudes.csv",
      "alpha_magnitudes.csv", if (plots) "figures/Figure_4_2_alpha_histograms.pdf" else "not generated (plots=FALSE)",
      "joint_tests.csv; full_components.csv"),
    scope = c(rep("china and us_ff3", 6), "us_ff5"),
    fields = c("PY, MAX, COM, max-k (K=30), AdaL: p_value",
      "T1, T5, T_quarter=ceil(0.25N), T_half=ceil(0.50N): p_value",
      "k, statistic, exceedances, p_value, clipped_p, weight",
      "median_bp, q05_bp, q95_bp, mean_abs_bp, max_abs_bp, share_abs_ge_50bp; plus mean_bp and sd_bp",
      "BH_005 and BH_010 from individual two-sided p-values",
      "monthly alpha in basis points and standardized_alpha; pointwise +/-1.96 on standardized panels",
      "AdaL and full-grid components under FF5"), stringsAsFactors = FALSE)
  result$manuscript_mapping <- rbind(result$manuscript_mapping, data.frame(
    manuscript_label = "emp:tab:diagnostics", artifact = "diagnostic_table.csv; model_diagnostics.csv; residual_diagnostic_summary.csv",
    scope = "china and us_ff3; optional diagnostics reported only if computed",
    fields = "T, N, BP6/BP12/LB12/ARCH6 BH discoveries, descriptive static eigenvalue-gap ratio, OR-Bonferroni loading-stability p-value"))
  report_dir <- file.path(output_dir, "report")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(result)) .empirical_report_csv(result[[nm]], file.path(report_dir, paste0(nm, ".csv")))
  figure <- NULL
  has_main <- any(c("china", "us_ff3") %in% names(panels))
  if (plots && has_main) {
    dir.create(file.path(report_dir, "figures"), showWarnings = FALSE)
    figure <- file.path(report_dir, "figures", "Figure_4_2_alpha_histograms.pdf")
    .empirical_report_histograms(panels, figure)
  }
  if (plots && !has_main)
    result$manuscript_mapping$artifact[result$manuscript_mapping$manuscript_label == "emp:fig:alphas"] <- "not generated (no complete main-market setting)"
  .empirical_report_csv(result$manuscript_mapping, file.path(report_dir, "manuscript_mapping.csv"))
  notes <- c("LFPM empirical report from supplied backend CSV outputs",
    paste0("Main-market coverage: ", if (main_complete) "China and US FF3 available" else "PARTIAL: consult coverage.csv"),
    paste0("Reported settings: ", if (length(panels)) paste(names(panels), collapse = ", ") else "none"),
    "No historic manuscript estimates are inserted. Missing settings/methods/orders remain missing.",
    "Outputs summarize the supplied run; verify its backend metadata/provenance for the selected sample, model and bootstrap budget.",
    "The main-paper mapping assumes china and us_ff3 use the intended main designs; filenames alone do not verify that scientific specification.",
    "us_ff5 is reported separately as a factor-set robustness result.",
    "B is method-specific: PY/MAX/COM have no bootstrap budget and retain missing B; AdaL/T1 and max-k retain their respective recorded budgets.",
    "Monthly returns are decimals; 1 basis point equals 0.0001 monthly return. Share_abs_ge_50bp is a fraction, not a percent.",
    "SD uses the sample denominator N-1; 5th/95th percentiles use linear interpolation (R type 7 / pandas default).",
    "BH counts use <= 0.05 and <= 0.10. Saved BH_q is checked against BH adjustment of individual_p.",
    "Histogram +/-1.96 reference lines are pointwise normal references and are not multiplicity-adjusted thresholds.",
    "This report covers the joint tests, component profile, stock-level counts, economic magnitudes and alpha distribution figure.",
    "Optional model_diagnostics.csv and diagnostic_summary.csv feed diagnostic_table.csv when available; coverage flags show omissions.",
    "Residual diagnostics use the fitted main specification; loading stability uses the constant-loading specification. The static gap ratio is descriptive.",
    "Fitted-null calibration, break reselection and max-k sensitivity retain their dedicated backend outputs.")
  writeLines(notes, file.path(report_dir, "REPORT.txt"))
  result$output_dir <- report_dir
  result$figure <- figure
  result$main_complete <- main_complete
  result$provenance <- "Derived solely from supplied backend results; no historical manuscript values inserted"
  structure(result, class = "lfpm_empirical_report")
}

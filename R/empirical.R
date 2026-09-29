.lfpm_python_backend <- function(python) {
  if (length(python) != 1L || is.na(python) || !nzchar(python)) {
    stop("Choose python= pointing to Python 3.9+ with numpy, pandas, scipy, statsmodels and threadpoolctl installed.", call. = FALSE)
  }
  resolved <- if (file.exists(python)) normalizePath(python) else Sys.which(python)
  if (!nzchar(resolved)) stop("Python executable was not found: ", python, call. = FALSE)
  script <- system.file("python", "empirical_backend.py", package = "LFPM")
  if (!nzchar(script)) stop("Install LFPM before calling the Python empirical backend.", call. = FALSE)
  list(python = resolved, script = script)
}

.lfpm_empirical_count <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) || x < 1 || x != floor(x)) {
    stop(name, " must be a positive whole number.", call. = FALSE)
  }
  format(x, scientific = FALSE, trim = TRUE)
}

.lfpm_empirical_run <- function(python, arguments, output_dir) {
  backend <- .lfpm_python_backend(python)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  log <- file.path(output_dir, "python.log")
  status <- system2(backend$python,
    args = vapply(c(backend$script, arguments, "--output-dir", output_dir), shQuote, character(1)),
    stdout = log, stderr = log, env = "PYTHONDONTWRITEBYTECODE=1")
  if (!identical(as.integer(status), 0L)) {
    details <- if (file.exists(log)) tail(readLines(log, warn = FALSE), 15L) else "No Python log was created."
    stop("Empirical backend failed. Log: ", log, "\n", paste(details, collapse = "\n"), call. = FALSE)
  }
  output_dir
}

#' Test a complete factor-model panel with the empirical Python implementation
#'
#' Uses the empirical full-grid calibration from the final paper. The Python
#' backend preserves NumPy PCG64 streams and the original bootstrap arithmetic.
#' See the installed reproduction documents for dependencies and data sources.
#' @export
lfpm_empirical_test <- function(Y, F, B = 9999L, seed = 20260929,
                                python = Sys.which("python3"),
                                include_maxk = FALSE, B_maxk = B,
                                output_dir = tempfile("lfpm-empirical-")) {
  if (!is.matrix(Y) || !is.numeric(Y) || !is.matrix(F) || !is.numeric(F) ||
      nrow(Y) != nrow(F) || ncol(Y) < 40L || any(!is.finite(Y)) || any(!is.finite(F))) {
    stop("Y and F must be finite numeric matrices with matching rows and at least 40 assets in Y.", call. = FALSE)
  }
  if (ncol(F) < 1L || nrow(Y) <= ncol(F) + 5L) {
    stop("F must have at least one column and residual degrees of freedom must exceed four.", call. = FALSE)
  }
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) || seed < 0 || seed != floor(seed)) {
    stop("seed must be a nonnegative whole number.", call. = FALSE)
  }
  if (!is.logical(include_maxk) || length(include_maxk) != 1L || is.na(include_maxk)) {
    stop("include_maxk must be TRUE or FALSE.", call. = FALSE)
  }
  # Input matrices are transient caller data, never part of the installed package.
  temporary <- tempfile("lfpm-matrices-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  y_file <- file.path(temporary, "y.csv")
  f_file <- file.path(temporary, "f.csv")
  utils::write.table(Y, y_file, sep = ",", row.names = FALSE, col.names = FALSE)
  utils::write.table(F, f_file, sep = ",", row.names = FALSE, col.names = FALSE)
  arguments <- c("test", "--y", y_file, "--f", f_file,
    "--B", .lfpm_empirical_count(B, "B"),
    "--seed", format(seed, scientific = FALSE, trim = TRUE),
    "--B-maxk", .lfpm_empirical_count(B_maxk, "B_maxk"),
    if (!include_maxk) "--skip-maxk")
  output_dir <- .lfpm_empirical_run(python, arguments, output_dir)
  results <- utils::read.csv(file.path(output_dir, "results.csv"), stringsAsFactors = FALSE)
  answer <- list(statistic = stats::setNames(results$statistic[results$method == "AdaL"], "Cauchy"),
    p.value = results$p_value[results$method == "AdaL"],
    method = "Empirical full-grid AdaL with refitted residual wild bootstrap",
    parameter = c(B = B, N = ncol(Y), T = nrow(Y)),
    results = results,
    components = utils::read.csv(file.path(output_dir, "components.csv")),
    alphas = utils::read.csv(file.path(output_dir, "alphas.csv")),
    output_dir = output_dir,
    metadata = file.path(output_dir, "metadata.json"))
  answer
}

#' Reproduce the paper's empirical applications using externally supplied data
#' @export
reproduce_empirical <- function(data_dir, output_dir,
    analysis = c("all", "china", "us_ff3", "us_ff5", "break_reselection", "fitted_null", "diagnostics", "sensitivity"),
    python = Sys.which("python3"), B = NULL, B_maxk = 9999L,
    repetitions = 500L, inner_B = 999L, sensitivity_seed = 20260906, include_maxk = TRUE,
    diagnostics = TRUE, plots = TRUE) {
  analysis <- match.arg(analysis)
  if (length(sensitivity_seed) != 1L || !is.numeric(sensitivity_seed) || !is.finite(sensitivity_seed) || sensitivity_seed < 0 || sensitivity_seed != floor(sensitivity_seed)) {
    stop("sensitivity_seed must be a nonnegative whole number.", call. = FALSE)
  }
  if (length(data_dir) != 1L || !dir.exists(data_dir)) stop("data_dir must be an existing external data directory.", call. = FALSE)
  if (missing(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    stop("Choose an output_dir outside the installed package.", call. = FALSE)
  }
  for (name in c("include_maxk", "diagnostics", "plots")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  arguments <- c("reproduce", "--data-dir", normalizePath(data_dir), "--analysis", analysis,
    if (!is.null(B)) c("--B", .lfpm_empirical_count(B, "B")),
    "--B-maxk", .lfpm_empirical_count(B_maxk, "B_maxk"),
    "--repetitions", .lfpm_empirical_count(repetitions, "repetitions"),
    "--inner-B", .lfpm_empirical_count(inner_B, "inner_B"),
    "--sensitivity-seed", format(sensitivity_seed, scientific = FALSE, trim = TRUE),
    if (!include_maxk) "--skip-maxk", if (diagnostics) "--diagnostics", if (plots) "--plots")
  output_dir <- .lfpm_empirical_run(python, arguments, output_dir)
  report <- if (analysis %in% c("all", "china", "us_ff3", "us_ff5")) {
    report_empirical(output_dir, plots = plots)
  } else NULL
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  structure(list(analysis = analysis, output_dir = output_dir, files = files,
    report = report, provenance = file.path(output_dir, "provenance.json")),
    class = "lfpm_empirical_reproduction")
}

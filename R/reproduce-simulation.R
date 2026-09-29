# Standalone reconstruction of the 2026-09-29 manuscript simulation.
# Internal engine functions deliberately retain their audited random draw order.

.simulation_integer <- function(x, name, minimum = 1L, scalar = FALSE) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x != floor(x)) || any(x < minimum | x > .Machine$integer.max) ||
      (scalar && length(x) != 1L) || anyDuplicated(x)) {
    stop(name, " must contain distinct integers >= ", minimum, call. = FALSE)
  }
  as.integer(x)
}

#' Define a manuscript simulation grid
#' @export
simulation_design <- function(preset = c("smoke", "paper"), models = NULL,
                              distributions = NULL, size_N = NULL, size_T = NULL,
                              power_N = NULL, power_T = NULL, m = NULL) {
  preset <- match.arg(preset)
  paper <- preset == "paper"
  if (is.null(models)) models <- if (paper) 1:4 else 1L
  if (is.null(distributions)) distributions <- if (paper)
    c("Normal", "t5", "MixNorm", "Chisq") else "Normal"
  if (is.null(size_N)) size_N <- if (paper) c(100L, 200L, 400L) else 40L
  if (is.null(size_T)) size_T <- if (paper) c(100L, 200L) else 30L
  if (is.null(power_N)) power_N <- if (paper) 200L else 40L
  if (is.null(power_T)) power_T <- if (paper) 100L else 30L
  models <- .simulation_integer(models, "models")
  if (any(!models %in% 1:4)) stop("models must be selected from 1:4")
  if (!is.character(distributions) || !length(distributions) ||
      anyNA(distributions) || anyDuplicated(distributions) ||
      any(!distributions %in% c("Normal", "t5", "MixNorm", "Chisq"))) {
    stop("Invalid innovation distributions")
  }
  size_N <- .simulation_integer(size_N, "size_N", 40L)
  size_T <- .simulation_integer(size_T, "size_T", 10L)
  power_N <- .simulation_integer(power_N, "power_N", 40L, TRUE)
  power_T <- .simulation_integer(power_T, "power_T", 10L, TRUE)
  if (is.null(m)) m <- if (paper) seq_len(power_N) else
    unique(pmin(c(1L, 5L, 20L, 40L), power_N))
  m <- sort(.simulation_integer(m, "m"))
  if (any(m > power_N)) stop("m must not exceed power_N")
  grid <- function(N, T_size) {
    x <- expand.grid(model = models, distribution = distributions,
                     T_size = T_size, N = N, stringsAsFactors = FALSE)
    x$m <- 0L
    x
  }
  structure(list(preset = preset, size = grid(size_N, size_T),
                 power = grid(power_N, power_T), m = m),
            class = "lfpm_simulation_design")
}

.simulation_check_design <- function(design) {
  if (!inherits(design, "lfpm_simulation_design") ||
      !all(c("preset", "size", "power", "m") %in% names(design))) {
    stop("design must be returned by simulation_design()")
  }
  for (phase in c("size", "power")) {
    x <- design[[phase]]
    if (!is.data.frame(x) || !nrow(x) ||
        !identical(names(x), c("model", "distribution", "T_size", "N", "m")) ||
        anyNA(x) || anyDuplicated(x) || any(!x$model %in% 1:4) ||
        any(!x$distribution %in% c("Normal", "t5", "MixNorm", "Chisq")) ||
        any(x$m != 0L)) stop("Invalid ", phase, " design")
    .simulation_integer(unique(x$N), "N", 40L)
    .simulation_integer(unique(x$T_size), "T_size", 10L)
  }
  .simulation_integer(design$m, "m")
  if (length(unique(design$power$N)) != 1L ||
      length(unique(design$power$T_size)) != 1L ||
      any(design$m > min(design$power$N))) stop("Invalid power grid")
  invisible(TRUE)
}

.simulation_columns <- function() c(paste0("score_", METHODS), paste0("p_", METHODS))

.simulation_one <- function(i, setting, cfg, phase, m) {
  # Each call resets the seed before generating the panel/support/multipliers.
  # This recreates the previous four comparators and Ma-rule AdaL on exactly
  # the same random panel, without either engine changing the other's draws.
  if (phase == "power_curve") {
    answer <- one_power_curve(i, setting, cfg, m)
    ma <- ma_one_power_curve(i, setting, cfg, m)
    answer[, c("score_AdaL", "p_AdaL")] <- ma
  } else {
    answer <- one_repetition(i, setting, cfg, phase)
    answer[c("score_AdaL", "p_AdaL")] <-
      ma_one_repetition(i, setting, cfg, phase)
  }
  answer
}

.simulation_atomic_rds <- function(object, path) {
  temporary <- tempfile(".lfpm-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop("Cannot save checkpoint: ", path)
}

.simulation_atomic_csv <- function(object, path) {
  temporary <- tempfile(".lfpm-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop("Cannot save table: ", path)
}

.simulation_key <- function(s) paste0("M", s$model, "_", s$distribution,
                                      "_T", s$T_size, "_N", s$N, "_m", s$m)

.simulation_rank_p <- function(scores, reference) {
  answer <- vapply(seq_len(5L), function(j) {
    (1 + nrow(reference) - findInterval(scores[, j], sort(reference[, j]),
                                       left.open = TRUE)) / (nrow(reference) + 1)
  }, numeric(nrow(scores)))
  matrix(answer, nrow(scores), 5L, dimnames = list(NULL, METHODS))
}

.simulation_summary <- function(samples, reference, setting, phase) {
  summaries <- paired <- list()
  for (calibration in c("native", "equal_size")) {
    p <- if (calibration == "native") samples[, 6:10, drop = FALSE] else
      .simulation_rank_p(samples[, 1:5, drop = FALSE], reference[, 1:5, drop = FALSE])
    rejects <- p < 0.05
    for (j in seq_along(METHODS)) {
      summaries[[length(summaries) + 1L]] <- cbind(as.data.frame(setting),
        phase = phase, calibration = calibration, method = METHODS[j],
        as.data.frame(as.list(binomial_summary(rejects[, j]))))
    }
    for (j in 1:4) {
      delta <- as.numeric(rejects[, 5L]) - as.numeric(rejects[, j])
      se <- stats::sd(delta) / sqrt(length(delta))
      paired[[length(paired) + 1L]] <- cbind(as.data.frame(setting),
        phase = phase, calibration = calibration,
        comparison = paste("AdaL -", METHODS[j]), difference = mean(delta),
        mcse = se, lower95 = mean(delta) - 1.96 * se,
        upper95 = mean(delta) + 1.96 * se)
    }
  }
  list(summary = do.call(rbind, summaries), paired = do.call(rbind, paired))
}

.simulation_bands <- function(curves, reference, setting, m) {
  bands <- list(m001_010 = 1:10, m011_050 = 11:50,
                m051_200 = 51:200, m001_200 = 1:200)
  result <- paired <- list()
  for (calibration in c("native", "equal_size")) {
    rejects <- lapply(curves, function(x) {
      p <- if (calibration == "native") x[, 6:10, drop = FALSE] else
        .simulation_rank_p(x[, 1:5, drop = FALSE], reference[, 1:5, drop = FALSE])
      p < 0.05
    })
    for (band in names(bands)) {
      ids <- which(m %in% bands[[band]])
      if (!length(ids)) next
      # m values are correlated; average within independent repetition first.
      means <- t(vapply(rejects, function(x) colMeans(x[ids, , drop = FALSE]), numeric(5L)))
      n <- nrow(means)
      for (j in 1:5) {
        rate <- mean(means[, j]); se <- stats::sd(means[, j]) / sqrt(n)
        result[[length(result) + 1L]] <- data.frame(model = setting$model,
          distribution = setting$distribution, T_size = setting$T_size, N = setting$N,
          calibration = calibration, band = band, m_min = min(m[ids]),
          m_max = max(m[ids]), m_count = length(ids), method = METHODS[j],
          mean_power = rate, mcse = se, lower95 = rate - 1.96 * se,
          upper95 = rate + 1.96 * se, independent_repetitions = n)
      }
      for (j in 1:4) {
        delta <- means[, 5L] - means[, j]; se <- stats::sd(delta) / sqrt(n)
        paired[[length(paired) + 1L]] <- data.frame(model = setting$model,
          distribution = setting$distribution, calibration = calibration,
          band = band, comparison = paste("AdaL -", METHODS[j]),
          difference = mean(delta), mcse = se, lower95 = mean(delta) - 1.96 * se,
          upper95 = mean(delta) + 1.96 * se, independent_repetitions = n)
      }
    }
  }
  list(summary = do.call(rbind, result), paired = do.call(rbind, paired))
}

.simulation_presentation_tables <- function(result) {
  wide <- function(dat, id) {
    x <- unique(dat[, id, drop = FALSE])
    x <- x[do.call(order, x), , drop = FALSE]
    key <- function(z) do.call(paste, c(z[, id, drop = FALSE], sep = "|"))
    for (method in METHODS) {
      d <- dat[dat$method == method, , drop = FALSE]
      x[[method]] <- 100 * d$rate[match(key(x), key(d))]
    }
    x
  }
  tables <- list()
  if (nrow(result$size_summary)) {
    for (cal in c("native", "equal_size")) {
      d <- result$size_summary[result$size_summary$calibration == cal, , drop = FALSE]
      tables[[paste0("size_", cal, "_percent")]] <-
        wide(d, c("model", "distribution", "T_size", "N"))
    }
  }
  if (nrow(result$power_summary)) {
    d <- result$power_summary[result$power_summary$m %in% c(1, 5, 10, 20, 50, 100, 200), , drop = FALSE]
    if (nrow(d)) tables$power_selected_m_percent <-
      wide(d, c("model", "distribution", "T_size", "N", "m", "calibration"))
    bands <- result$power_band_summary
    pairs <- result$power_band_paired_differences
    claims <- list()
    grid <- unique(bands[bands$calibration == "equal_size" & bands$band != "m001_200",
                         c("model", "distribution", "band"), drop = FALSE])
    for (i in seq_len(nrow(grid))) {
      s <- grid[i, ]
      d <- bands[bands$model == s$model & bands$distribution == s$distribution &
                   bands$band == s$band & bands$calibration == "equal_size", , drop = FALSE]
      adal <- d[d$method == "AdaL", , drop = FALSE]
      others <- d[d$method != "AdaL", , drop = FALSE]
      best <- others[which.max(others$mean_power), , drop = FALSE]
      delta <- pairs[pairs$model == s$model & pairs$distribution == s$distribution &
        pairs$band == s$band & pairs$calibration == "equal_size" &
        pairs$comparison == paste("AdaL -", best$method), , drop = FALSE]
      claims[[i]] <- data.frame(model = s$model, distribution = s$distribution,
        band = s$band, comparator = best$method, adal_mean = adal$mean_power,
        comparator_mean = best$mean_power, paired_difference = delta$difference,
        lower95 = delta$lower95, upper95 = delta$upper95,
        selection = "descriptive maximum among four fixed comparators; intervals unadjusted")
    }
    if (length(claims)) tables$manuscript_claims <- do.call(rbind, claims)
  }
  tables
}

.simulation_engine_signature <- function() {
  functions <- c("make_covariance", "make_innovations", "make_factors", "make_data",
    "fit_panel", "top_sums", "row_top_sums", "precision_estimate", "bootstrap_components",
    "test_panel", "seed_for", "one_repetition", "prepare_tests", "score_alphas",
    "right_tail_from_sorted", "one_power_curve", "ma_candidate_k", "ma_bootstrap_L",
    "ma_prepare", "ma_combine_p", "ma_score_alphas", "ma_one_repetition",
    "ma_one_power_curve", ".simulation_one", ".simulation_rank_p",
    ".simulation_summary", ".simulation_bands", ".simulation_presentation_tables",
    "binomial_summary")
  list(functions = setNames(lapply(functions, function(nm)
    paste(deparse(get(nm, mode = "function")), collapse = "\n")), functions),
    constants = list(methods = METHODS, endpoint = MA_P_ENDPOINT_EPS,
                     bootstrap_variance = MA_BOOTSTRAP_VARIANCE))
}

#' Recompute the chapter 3 simulation with restartable checkpoints
#' @export
reproduce_simulation <- function(output_dir, preset = c("smoke", "paper"),
                                 design = NULL, R = NULL, Rcal = NULL, B = NULL,
                                 seed = 20260927L, workers = 1L, chunk = 25L,
                                 resume = TRUE, phase = c("all", "calibration", "size", "power"),
                                 quiet = FALSE) {
  preset <- match.arg(preset)
  if (is.null(design)) design <- simulation_design(preset)
  .simulation_check_design(design)
  # A supplied design determines budget defaults, so a paper design cannot
  # silently inherit the smoke budget from the function's default argument.
  if (!missing(preset) && preset != design$preset) stop("preset and design$preset disagree")
  preset <- design$preset
  phase <- match.arg(phase)
  if (is.null(R)) R <- if (preset == "paper") 1000L else 6L
  if (is.null(Rcal)) Rcal <- if (preset == "paper") 1999L else 39L
  if (is.null(B)) B <- if (preset == "paper") 999L else 39L
  R <- .simulation_integer(R, "R", 2L, TRUE)
  Rcal <- .simulation_integer(Rcal, "Rcal", 2L, TRUE)
  B <- .simulation_integer(B, "B", 2L, TRUE)
  seed <- .simulation_integer(seed, "seed", 0L, TRUE)
  workers <- .simulation_integer(workers, "workers", 1L, TRUE)
  chunk <- .simulation_integer(chunk, "chunk", 1L, TRUE)
  for (flag in c("resume", "quiet")) {
    value <- get(flag)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) stop(flag, " must be TRUE or FALSE")
  }
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))
    stop("output_dir must be one directory path")
  if (preset == "smoke" && !quiet) message("SMOKE RUN: reduced dimensions/budgets; not manuscript results.")
  cfg <- list(engine = "LFPM-manuscript-20260929-v1", sim = R,
              calibration = Rcal, B = B, seed = seed, signal_energy = 10,
              bootstrap_variance = "refit", adal_include_l1 = FALSE,
              significance = 0.05, rejection_rule = "p < 0.05")
  manifest <- list(config = cfg, design = design,
    rng = c("Mersenne-Twister", "Inversion", "Rejection"),
    R_version = as.character(getRversion()), engine = .simulation_engine_signature())
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  manifest_path <- file.path(output_dir, "manifest.rds")
  if (file.exists(manifest_path)) {
    if (!resume) stop("Output exists; use resume=TRUE or a new output directory")
    if (!identical(readRDS(manifest_path), manifest))
      stop("Output directory has an incompatible configuration, design, R version or engine. Use a new output directory.")
  } else {
    if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
      stop("New output directory must be empty")
    .simulation_atomic_rds(manifest, manifest_path)
    utils::capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  }
  checkpoint_dir <- file.path(output_dir, "checkpoints")
  dir.create(checkpoint_dir, showWarnings = FALSE)
  old_rng <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_rng))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  do.call(RNGkind, as.list(manifest$rng))
  if (.Platform$OS.type == "windows" && workers > 1L) {
    if (!quiet) message("Windows uses sequential execution; checkpoints remain compatible.")
    workers <- 1L
  }
  run_samples <- function(s, what) {
    curve <- what == "power_curve"
    total <- if (what == "calibration") Rcal else R
    path <- file.path(checkpoint_dir, paste0(what, "_", .simulation_key(s), ".rds"))
    expected <- list(config = cfg, setting = s, phase = what,
                     m = if (curve) design$m else NULL)
    if (file.exists(path)) {
      state <- readRDS(path)
      if (!identical(state$metadata, expected) || !is.list(state$samples) ||
          length(state$samples) > total ||
          !identical(state$rep_ids, seq_along(state$samples))) stop("Checkpoint mismatch: ", path)
    } else state <- list(metadata = expected, rep_ids = integer(), samples = list())
    validate <- function(x) {
      shape <- if (curve) identical(dim(x), c(length(design$m), 10L)) &&
        identical(rownames(x), as.character(design$m)) &&
        identical(colnames(x), .simulation_columns()) else
        is.numeric(x) && is.null(dim(x)) && identical(names(x), .simulation_columns())
      isTRUE(shape) && all(is.finite(x))
    }
    if (!all(vapply(state$samples, validate, logical(1L)))) stop("Invalid stored sample: ", path)
    done <- length(state$samples)
    if (done < total) {
      if (!quiet) message(what, " ", .simulation_key(s), ": ", done, "/", total, " cached")
      for (start in seq.int(done + 1L, total, by = chunk)) {
        ids <- seq.int(start, min(start + chunk - 1L, total))
        work <- function(i) .simulation_one(i, s, cfg, what, design$m)
        added <- if (workers > 1L) parallel::mclapply(ids, work,
          mc.cores = min(workers, length(ids)), mc.set.seed = FALSE) else lapply(ids, work)
        if (!all(vapply(added, validate, logical(1L))))
          stop("Simulation failed; previous checkpoint retained. No repetitions were omitted.")
        state$samples <- c(state$samples, added)
        state$rep_ids <- seq_along(state$samples)
        .simulation_atomic_rds(state, path)
        .simulation_atomic_rds(list(phase = what, setting = s,
          done = length(state$samples), total = total, time = Sys.time()),
          file.path(output_dir, "progress.rds"))
        if (!quiet) message("  ", length(state$samples), "/", total)
      }
    }
    if (curve) state$samples else do.call(rbind, state$samples)
  }
  setting <- function(grid, i) as.list(grid[i, , drop = FALSE])
  result <- list(size_summary = list(), size_paired_differences = list(),
    power_summary = list(), power_paired_differences = list(),
    power_band_summary = list(), power_band_paired_differences = list())
  write_tables <- function(prefix) {
    for (nm in names(result)[startsWith(names(result), prefix)]) {
      if (length(result[[nm]])) .simulation_atomic_csv(do.call(rbind, result[[nm]]),
        file.path(output_dir, paste0(nm, ".csv")))
    }
  }
  if (phase == "calibration") {
    grid <- unique(rbind(design$size, design$power))
    for (i in seq_len(nrow(grid))) run_samples(setting(grid, i), "calibration")
  }
  if (phase %in% c("all", "size")) {
    for (i in seq_len(nrow(design$size))) {
      s <- setting(design$size, i)
      cal <- run_samples(s, "calibration")
      sam <- run_samples(s, "size")
      z <- .simulation_summary(sam, cal, s, "size")
      result$size_summary[[i]] <- z$summary
      result$size_paired_differences[[i]] <- z$paired
      write_tables("size_")
    }
  }
  if (phase %in% c("all", "power")) {
    for (i in seq_len(nrow(design$power))) {
      s <- setting(design$power, i)
      cal <- run_samples(s, "calibration")
      curves <- run_samples(s, "power_curve")
      for (j in seq_along(design$m)) {
        sam <- t(vapply(curves, function(x) x[j, ], numeric(10L)))
        s$m <- design$m[j]
        z <- .simulation_summary(sam, cal, s, "power")
        key <- .simulation_key(s)
        result$power_summary[[key]] <- z$summary
        result$power_paired_differences[[key]] <- z$paired
      }
      z <- .simulation_bands(curves, cal, s, design$m)
      result$power_band_summary[[i]] <- z$summary
      result$power_band_paired_differences[[i]] <- z$paired
      write_tables("power_")
    }
  }
  result <- lapply(result, function(x) if (length(x)) do.call(rbind, x) else data.frame())
  presentation <- .simulation_presentation_tables(result)
  for (nm in names(presentation)) .simulation_atomic_csv(presentation[[nm]],
    file.path(output_dir, paste0(nm, ".csv")))
  result <- c(result, presentation)
  result$output_dir <- output_dir
  result$manifest <- manifest
  result$provenance <- "recomputed from random draws, not bundled reference tables"
  .simulation_atomic_rds(list(phase = "complete", requested_phase = phase,
                              time = Sys.time()), file.path(output_dir, "progress.rds"))
  structure(result, class = "lfpm_simulation")
}

#' Read previously computed manuscript simulation reference tables
#' @export
reproduction_reference <- function(table = c("all", "size_summary", "power_summary",
                                    "size_paired_differences", "power_paired_differences",
                                    "power_band_summary", "power_band_paired_differences",
                                    "size_native_percent", "size_equal_size_percent",
                                    "power_selected_m_percent", "manuscript_claims",
                                    "validation_summary")) {
  table <- match.arg(table)
  path <- system.file("reference", "simulation", package = "LFPM", mustWork = TRUE)
  tables <- sub("\\.csv\\.gz$", "", list.files(path, pattern = "\\.csv\\.gz$"))
  read_table <- function(nm) {
    con <- gzfile(file.path(path, paste0(nm, ".csv.gz")), "rt")
    on.exit(close(con), add = TRUE)
    x <- utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE)
    attr(x, "provenance") <- "Bundled cached manuscript summary; no simulation was run"
    x
  }
  if (table != "all") return(read_table(table))
  x <- setNames(lapply(tables, read_table), tables)
  x$provenance <- "Bundled cached manuscript summaries; no simulation was run"
  structure(x, class = "lfpm_simulation_reference")
}

#' Plot simulation results or bundled reference tables
#' @export
plot_simulation <- function(x, file = NULL, phase = c("power", "size"),
                            calibration = c("equal_size", "native"),
                            models = NULL, distributions = NULL,
                            width = NULL, height = NULL) {
  phase <- match.arg(phase)
  calibration <- match.arg(calibration)
  if (is.data.frame(x)) {
    dat <- x
    provenance <- attr(x, "provenance")
  } else if (is.list(x)) {
    dat <- x[[paste0(phase, "_summary")]]
    provenance <- x$provenance
  } else stop("x must be a simulation result, reference list, or summary data frame")
  needed <- c("model", "distribution", "T_size", "N", "m", "calibration", "method", "rate", "lower95", "upper95")
  if (!is.data.frame(dat) || !all(needed %in% names(dat))) stop("No valid ", phase, " summary in x")
  dat <- dat[dat$calibration == calibration, , drop = FALSE]
  if (!is.null(models)) dat <- dat[dat$model %in% models, , drop = FALSE]
  if (!is.null(distributions)) dat <- dat[dat$distribution %in% distributions, , drop = FALSE]
  if (!nrow(dat)) stop("No scenarios selected")
  panels <- unique(dat[, if (phase == "power") c("model", "distribution", "T_size", "N") else
    c("model", "distribution"), drop = FALSE])
  ncol <- min(4L, length(unique(panels$distribution)))
  nrow <- ceiling(nrow(panels) / ncol)
  if (is.null(width)) width <- max(6, 3.2 * ncol)
  if (is.null(height)) height <- max(4.5, 2.8 * nrow + 1)
  if (!is.null(file)) {
    if (!is.character(file) || length(file) != 1L || !grepl("\\.pdf$", file, ignore.case = TRUE))
      stop("file must be a PDF filename, or NULL for the current device")
    grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE, after = FALSE)
  graphics::par(mfrow = c(nrow, ncol), mar = c(3.5, 3.6, 2.2, 1),
                oma = c(2.1, 0, 2.5, 0), mgp = c(2.1, .65, 0))
  colours <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00")
  types <- c(1, 2, 3, 4, 1)
  for (i in seq_len(nrow(panels))) {
    panel <- panels[i, ]
    keep <- rep(TRUE, nrow(dat))
    for (nm in names(panels)) keep <- keep & dat[[nm]] == panel[[nm]]
    d <- dat[keep, , drop = FALSE]
    title <- paste0("M", panel$model, " | ", panel$distribution)
    if (phase == "power") {
      title <- paste0(title, " | T=", panel$T_size, ", N=", panel$N)
      xrange <- range(d$m)
      if (diff(xrange) == 0) xrange <- xrange + c(-0.5, 0.5)
      graphics::plot(xrange, c(0, 1), type = "n", xlab = "Nonzero intercepts (m)",
                     ylab = "Rejection probability", main = title, cex.main = .85)
      graphics::abline(h = c(.25, .5, .75), col = "grey93")
      for (j in seq_along(METHODS)) {
        z <- d[d$method == METHODS[j], ]; z <- z[order(z$m), ]
        graphics::lines(z$m, z$rate, col = colours[j], lty = types[j],
                        lwd = if (j == 5L) 1.8 else 1.2, type = if (nrow(z) == 1L) "p" else "l")
      }
    } else {
      settings <- unique(d[, c("T_size", "N")]); settings <- settings[order(settings$T_size, settings$N), ]
      keys <- paste(settings$T_size, settings$N)
      ymax <- min(1, max(.1, d$upper95 * 1.04))
      graphics::plot(c(.5, nrow(settings) + .5), c(0, ymax), type = "n",
        xlab = "(T, N)", ylab = "Null rejection probability", xaxt = "n",
        main = title, cex.main = .85)
      graphics::axis(1, seq_along(keys), paste0(settings$T_size, ",", settings$N), cex.axis = .7)
      graphics::abline(h = .05, lty = 2, col = "grey40")
      for (j in seq_along(METHODS)) {
        z <- d[d$method == METHODS[j], ]; z <- z[match(keys, paste(z$T_size, z$N)), ]
        xp <- seq_along(keys) + (j - 3) * .11
        graphics::segments(xp, z$lower95, xp, z$upper95, col = colours[j])
        graphics::points(xp, z$rate, col = colours[j], pch = j, cex = .75)
      }
    }
    if (i == 1L) graphics::legend("topright", METHODS, col = colours,
      lty = if (phase == "power") types else NULL,
      pch = if (phase == "size") 1:5 else NULL, cex = .65, bty = "n")
  }
  cached <- is.character(provenance) && grepl("[Cc]ached|[Bb]undled", provenance)
  budget <- if (is.list(x) && !is.null(x$manifest)) {
    cfg <- x$manifest$config
    paste0(" | R=", cfg$sim, ", Rcal=", cfg$calibration, ", B=", cfg$B)
  } else if (cached) " | R=1000, Rcal=1999, B=999" else ""
  graphics::mtext(paste0(if (cached) "Cached manuscript tables" else "Simulation results",
    " | ", phase, " | ", calibration, budget), outer = TRUE, side = 3, line = 1, cex = .85)
  graphics::mtext(if (phase == "size")
    "Pointwise Wilson 95% intervals; equal-size uncertainty is conditional on calibration."
    else "Pointwise and paired uncertainty are in the CSV tables; no simultaneous bands are implied.",
    outer = TRUE, side = 1, line = .8, cex = .65)
  invisible(dat)
}

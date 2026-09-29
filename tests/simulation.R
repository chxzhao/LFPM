library(LFPM)

internal <- function(name) getFromNamespace(name, "LFPM")
expect_error <- function(expr) stopifnot(inherits(try(expr, silent = TRUE), "try-error"))
paper_design <- simulation_design("paper")
stopifnot(nrow(paper_design$size) == 96L,
          nrow(paper_design$power) * length(paper_design$m) == 3200L,
          identical(paper_design$m, 1:200))
expect_error(simulation_design(power_N = 30L))
expect_error(simulation_design(m = c(1, 1)))

out <- tempfile("lfpm-simulation-")
design <- simulation_design(m = c(1L, 5L, 40L))
set.seed(987L)
before <- .Random.seed
run <- reproduce_simulation(out, design = design, R = 4, Rcal = 7, B = 19,
                             chunk = 2, quiet = TRUE)
stopifnot(identical(before, .Random.seed), nrow(run$size_summary) == 10L,
          nrow(run$power_summary) == 30L)

# Every empirical rate and Wilson interval is computed from the actual saved
# independent-repetition indicators, including p==0.05's strict non-rejection.
size_file <- file.path(out, "checkpoints", "size_M1_Normal_T30_N40_m0.rds")
cal_file <- file.path(out, "checkpoints", "calibration_M1_Normal_T30_N40_m0.rds")
size <- readRDS(size_file)
cal <- readRDS(cal_file)
size_mat <- do.call(rbind, size$samples)
cal_mat <- do.call(rbind, cal$samples)
stopifnot(identical(size$rep_ids, 1:4), identical(cal$rep_ids, 1:7))
for (j in seq_len(5L)) {
  native <- size_mat[, 5L + j] < .05
  equal <- vapply(size_mat[, j], function(v)
    (1 + sum(cal_mat[, j] >= v)) / (nrow(cal_mat) + 1) < .05, logical(1L))
  for (mode in c("native", "equal_size")) {
    rejected <- if (mode == "native") native else equal
    # Avoid depending on row order in the result table.
    method <- sub("^score_", "", colnames(size_mat)[j])
    row <- run$size_summary[run$size_summary$method == method &
                             run$size_summary$calibration == mode, ]
    p <- mean(rejected); n <- length(rejected); z <- qnorm(.975)
    center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    stopifnot(row$rate == p, row$mcse == sqrt(p * (1 - p) / n),
              isTRUE(all.equal(row$lower95, center - half)),
              isTRUE(all.equal(row$upper95, center + half)))
  }
}

# The merged implementation must reproduce the retained four comparators
# and Ma implementation separately from the same original per-scenario seed.
s <- as.list(design$size[1L, ])
cfg <- run$manifest$config
original <- internal("one_repetition")(1L, s, cfg, "size")
ma <- internal("ma_one_repetition")(1L, s, cfg, "size")
keep <- c(1:4, 6:9)
stopifnot(identical(size$samples[[1L]][keep], original[keep]),
          identical(size$samples[[1L]][c(5L, 10L)], ma),
          internal("seed_for")(cfg$seed, "size_1_Normal_30_40_0", 1) !=
            internal("seed_for")(cfg$seed, "calibration_1_Normal_30_40_0", 1))
full <- internal(".simulation_one")(2L, s, cfg, "power_curve", c(1L, 5L, 40L))
subset <- internal(".simulation_one")(2L, s, cfg, "power_curve", c(1L, 40L))
stopifnot(identical(full[c(1L, 3L), , drop = FALSE], subset))

# Truncating a valid checkpoint emulates interruption after a completed
# chunk; resumption must reconstruct the original repetitions and summaries.
complete_size <- size
size$samples <- size$samples[1:2]; size$rep_ids <- 1:2
saveRDS(size, size_file)
resumed <- reproduce_simulation(out, design = design, R = 4, Rcal = 7, B = 19,
                                chunk = 3, quiet = TRUE)
stopifnot(identical(readRDS(size_file), complete_size),
          identical(run$size_summary, resumed$size_summary),
          identical(run$power_summary, resumed$power_summary))
expect_error(reproduce_simulation(out, design = design, R = 5, Rcal = 7, B = 19, quiet = TRUE))
expect_error(reproduce_simulation(out, design = design, R = 4, Rcal = 8, B = 19, quiet = TRUE))
expect_error(reproduce_simulation(out, design = design, R = 4, Rcal = 7, B = 20, quiet = TRUE))
expect_error(reproduce_simulation(out, design = design, R = 4, Rcal = 7, B = 19,
                                  seed = 123, quiet = TRUE))
expect_error(reproduce_simulation(out, design = design, R = 4, Rcal = 7, B = 19,
                                  resume = FALSE, quiet = TRUE))

# A perfectly dependent pair of m values has the same SE as one value;
# treating the two m points as independent would incorrectly reduce it.
curves <- list(matrix(rep(c(rep(0, 5), rep(.01, 5)), each = 2), 2, 10),
               matrix(rep(c(rep(0, 5), rep(.99, 5)), each = 2), 2, 10),
               matrix(rep(c(rep(0, 5), rep(.99, 5)), each = 2), 2, 10))
bands <- internal(".simulation_bands")(curves, cal_mat, s, c(1L, 2L))$summary
native_bands <- bands[bands$calibration == "native", ]
stopifnot(all(abs(native_bands$mcse - sd(c(1, 0, 0)) / sqrt(3)) < 1e-12))

reference <- reproduction_reference()
stopifnot(nrow(reference$size_summary) == 960L,
          nrow(reference$power_summary) == 32000L,
          all(reference$power_summary$repetitions == 1000L),
          grepl("cached", reference$provenance))
derived <- internal(".simulation_presentation_tables")(reference)
for (name in names(derived)) {
  a <- derived[[name]]; b <- reference[[name]]
  if (name == "manuscript_claims") {
    a <- a[order(a$model, a$distribution, a$band), ]
    b <- b[order(b$model, b$distribution, b$band), ]
  }
  stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE)))
}
pdf <- file.path(out, "smoke.pdf")
plot_simulation(run, pdf)
stopifnot(file.exists(pdf), file.info(pdf)$size > 0)
unlink(out, recursive = TRUE)

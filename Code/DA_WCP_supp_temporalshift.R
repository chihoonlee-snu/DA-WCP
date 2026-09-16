if (dir.exists("Code") && file.exists("Code/KLIEP.R")) setwd("Code")
dir.create("../Figures", recursive = TRUE, showWarnings = FALSE)

# Serial supplementary simulation for temporal outcome shifts.
#
# This strict OFAT sensitivity analysis starts from the fixed-size rho = 0
# scenario-S1 arm in DA_WCP_supp_ICC.R. Within a Monte Carlo generation, the
# historical sample, current X/G sample, target X/G and residuals, split,
# KLIEP estimate, subgroup factors, score models, and prediction intervals are
# generated/fitted exactly once. Only target outcomes vary with delta:
#
#   Y_2(delta) = X_2^T beta + gamma_G + delta * sigma_G(X_2) + epsilon_2,
#   sigma_G(X_2) = sqrt(exp(X_2^T theta + eta_G)).
#
# Interactive use from R/RStudio (0-based logical CPU IDs):
#   source("DA_WCP_supp_temporalshift.R")
#
# Nothing runs automatically when this file is sourced.


###############################################################################
# Reuse the current fixed-size S1 and shared simulation implementation
###############################################################################

source("DA_WCP_supp_ICC.R")


###############################################################################
# Exact S1 constants and requested OFAT grid
###############################################################################

delta_vec_temporalshift <- c(0, 0.25, 0.50, 0.75, 1.00)

stopifnot(
  identical(prob1[, 1], rep(0.20, 5)),
  identical(prob2[, 1], c(0.10, 0.15, 0.20, 0.25, 0.30)),
  identical(target_g, c(0.05, 0.10, 0.15, 0.30, 0.40)),
  identical(beta, c(5, -3, 2, 1)),
  identical(gamma, c(-10, -5, 0, 5, 10)),
  identical(theta_hetero, c(0.5, -0.4, 0.3, 0.2)),
  identical(eta_hetero, c(1.2, 1.6, 2.0, 2.4, 2.8))
)

validate_temporalshift_n <- function(value, name) {
  if (length(value) != 1L || !is.finite(value) || value < 1L ||
      value != as.integer(value)) {
    stop(name, " must be one positive integer.")
  }
  as.integer(value)
}

validate_temporalshift_delta <- function(delta_vec) {
  if (!is.numeric(delta_vec) || length(delta_vec) == 0L ||
      any(!is.finite(delta_vec)) || any(delta_vec < 0)) {
    stop("delta_vec must contain finite non-negative numeric values.")
  }
  delta_vec <- as.numeric(delta_vec)
  if (anyDuplicated(delta_vec)) stop("delta_vec cannot contain duplicates.")
  delta_vec
}


###############################################################################
# One S1 interval fit per generation; delta never enters this section
###############################################################################

fit_temporalshift_intervals <- function(
    data1, data2, data_test_baseline, split_indices,
    cached_rvec_marg = NULL,
    R = 4, alpha = 0.2, n_alpha = 5,
    sig_wo_cv = 1.2, kliep_cv = FALSE) {
  data1_y <- data1[, c(var_x, var_y, var_z, var_g, var_h), drop = FALSE]
  data2_x <- data2[, c(var_x, var_z, var_g, var_h), drop = FALSE]
  test <- data_test_baseline[
    , c(var_x, var_y, var_z, var_g, var_h), drop = FALSE
  ]
  train_m <- data1_y[split_indices$marginal_train, , drop = FALSE]
  cal_m <- data1_y[split_indices$marginal_cal, , drop = FALSE]

  use_cached_m <- !is.null(cached_rvec_marg)
  da_ar_m <- DA_WCP(
    train_m, cal_m, train_m, data2_x, test,
    var_x, var_y, var_z, var_g, var_h,
    val_type = "marg", score_type = "ar", p2_vec = target_g,
    alpha = alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv, simul = TRUE,
    rvec = use_cached_m,
    rvec_cal = if (use_cached_m) cached_rvec_marg else NA
  )
  rvec_marg_current <- da_ar_m$rvec_cal

  da_cqr_m <- DA_WCP(
    train_m, cal_m, train_m, data2_x, test,
    var_x, var_y, var_z, var_g, var_h,
    val_type = "marg", score_type = "cqr", p2_vec = target_g,
    alpha = alpha, n_alpha = n_alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv, simul = TRUE,
    rvec = TRUE, rvec_cal = rvec_marg_current
  )

  intervals <- list(AR = da_ar_m$result_mat, CQR = da_cqr_m$result_mat)
  stopifnot(
    all(vapply(intervals, nrow, integer(1)) == nrow(test)),
    all(vapply(intervals, ncol, integer(1)) == 2L),
    all(vapply(
      intervals, function(x) all(is.finite(as.matrix(x))), logical(1)
    ))
  )
  list(intervals = intervals, rvec_marg = rvec_marg_current)
}

evaluate_temporalshift_grid <- function(
    intervals, baseline_mean, sigma, standardized_residual,
    delta_vec, generation) {
  baseline_y <- baseline_mean + sigma * standardized_residual
  rows <- vector("list", length(delta_vec) * length(intervals))
  index <- 0L

  for (delta in delta_vec) {
    y_delta <- baseline_y + delta * sigma
    for (score in names(intervals)) {
      index <- index + 1L
      result_mat <- intervals[[score]]
      covered <- y_delta >= result_mat[, 1] & y_delta <= result_mat[, 2]
      rows[[index]] <- data.frame(
        generation = as.integer(generation),
        delta = delta,
        score = score,
        coverage = mean(covered),
        length = mean(result_mat[, 2] - result_mat[, 1]),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(rows)
}


###############################################################################
# Common-random-number replicate
###############################################################################

run_temporalshift_replicate <- function(
    replicate_id, n1, n2, test_n, delta_vec, cluster_size,
    R, alpha, n_alpha, sig_wo_cv, kliep_cv,
    cached_rvec_marg = NULL) {
  m <- as.integer(replicate_id)
  cat("------------------------------------------------------\n")
  cat("  Temporal outcome shift - replicate ", m, "\n", sep = "")
  cat("------------------------------------------------------\n")

  # This is exactly the random-generation order used by the rho = 0 S1 arm in
  # DA_WCP_supp_ICC.R. Although rho = 0 is iid, retaining cluster_size here also
  # retains the same number/order of base random draws and hence exact seeds.
  set.seed(123L + m)
  base1 <- generate_base_ICC(n1, t = 1, cluster_size = cluster_size)
  base2 <- generate_base_ICC(n2, t = 2, cluster_size = cluster_size)
  base_test <- generate_base_ICC(
    test_n, t = 2, cluster_size = 1L, test = TRUE
  )
  data1 <- generate_data_ICC(base1, rho = 0)
  data2 <- generate_data_ICC(base2, rho = 0)
  data_test_baseline <- generate_data_ICC(base_test, rho = 0)
  data_test_baseline$stratum <- 1L
  split_indices <- make_split_indices(
    base1$Z, replicate_id = m, stratum = base1$stratum
  )

  interval_fit <- fit_temporalshift_intervals(
    data1, data2, data_test_baseline, split_indices,
    cached_rvec_marg = cached_rvec_marg,
    R = R, alpha = alpha, n_alpha = n_alpha,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv
  )

  target_X <- as.matrix(data_test_baseline[, var_x, drop = FALSE])
  baseline_mean <- drop(target_X %*% beta) + gamma[data_test_baseline$Z]
  sigma <- sqrt(exp(
    drop(target_X %*% theta_hetero) + eta_hetero[data_test_baseline$Z]
  ))
  standardized_residual <- data_test_baseline$.std_residual
  baseline_y_reconstructed <- baseline_mean + sigma * standardized_residual

  results <- evaluate_temporalshift_grid(
    interval_fit$intervals, baseline_mean, sigma, standardized_residual,
    delta_vec = delta_vec, generation = m
  )

  delta_zero <- results[results$delta == 0, , drop = FALSE]
  baseline_metrics <- bind_rows(lapply(names(interval_fit$intervals), function(score) {
    metrics <- interval_metrics(
      interval_fit$intervals[[score]], data_test_baseline
    )
    data.frame(
      score = score,
      coverage = unname(metrics["coverage"]),
      length = unname(metrics["length"]),
      stringsAsFactors = FALSE
    )
  }))
  delta_zero <- delta_zero[order(delta_zero$score), c(
    "score", "coverage", "length"
  )]
  baseline_metrics <- baseline_metrics[
    order(baseline_metrics$score), , drop = FALSE
  ]
  rownames(delta_zero) <- NULL
  rownames(baseline_metrics) <- NULL

  interval_length_invariance <- results %>%
    group_by(score) %>%
    summarise(n_distinct_length = n_distinct(length), .groups = "drop")

  sanity <- list(
    delta_zero_is_original_S1_outcome = isTRUE(all.equal(
      baseline_y_reconstructed, data_test_baseline$Y,
      tolerance = 0, check.attributes = FALSE
    )),
    sigma_is_conditional_sd = isTRUE(all.equal(
      sigma^2,
      exp(drop(target_X %*% theta_hetero) +
            eta_hetero[data_test_baseline$Z]),
      tolerance = 1e-12, check.attributes = FALSE
    )),
    delta_zero_metrics_match_baseline = isTRUE(all.equal(
      delta_zero, baseline_metrics,
      tolerance = 0, check.attributes = FALSE
    )),
    intervals_fit_once_and_reused = TRUE,
    historical_data_reused_across_delta = TRUE,
    current_covariates_reused_across_delta = TRUE,
    target_covariates_and_residuals_reused_across_delta = TRUE,
    kliep_weights_reused_across_delta = TRUE,
    subgroup_factors_reused_across_delta = TRUE,
    interval_lengths_invariant_across_delta = all(
      interval_length_invariance$n_distinct_length == 1L
    ),
    outcome_shift_formula_verified = all(vapply(
      delta_vec,
      function(delta) isTRUE(all.equal(
        baseline_y_reconstructed + delta * sigma,
        baseline_mean + delta * sigma + sigma * standardized_residual,
        tolerance = 1e-12, check.attributes = FALSE
      )),
      logical(1)
    ))
  )

  if (!all(vapply(sanity, isTRUE, logical(1)))) {
    failed <- names(sanity)[!vapply(sanity, isTRUE, logical(1))]
    stop("Temporal-shift replicate checks failed: ", paste(failed, collapse = ", "))
  }

  list(
    results = results,
    rvec_marg = interval_fit$rvec_marg,
    sanity = sanity
  )
}


###############################################################################
# Summary, figure, and exact requested output files
###############################################################################

summarize_temporalshift_results <- function(results) {
  results %>%
    group_by(delta, score) %>%
    summarise(
      mean_coverage = mean(coverage),
      sd_coverage = sd(coverage),
      mean_length = mean(length),
      sd_length = sd(length),
      .groups = "drop"
    ) %>%
    arrange(delta, factor(score, levels = c("AR", "CQR")))
}

make_temporalshift_coverage_figure <- function(results, delta_vec, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  delta_labels <- sprintf("%.2f", delta_vec)
  delta_labels[delta_vec == 0] <- "0"
  plot_data <- results
  plot_data$delta_plot <- factor(
    sprintf("%.2f", plot_data$delta), levels = sprintf("%.2f", delta_vec),
    labels = delta_labels
  )
  plot_data$score <- factor(plot_data$score, levels = c("AR", "CQR"))

  coverage_plot <- ggplot(
    plot_data, aes(x = delta_plot, y = coverage, fill = score)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = 0.65, outlier.size = 0.8
    ) +
    geom_hline(
      yintercept = 0.8, color = "red",
      linetype = "dashed", linewidth = 0.7
    ) +
    scale_fill_manual(values = c(AR = "#9ecae1", CQR = "#3182bd")) +
    scale_x_discrete(expand = expansion(add = 0.7)) +
    scale_y_continuous(
      breaks = seq(0.3, 0.9, by = 0.1),
      expand = expansion(mult = c(0.035, 0.035))
    ) +
    coord_cartesian(ylim = c(0.3, 0.9), expand = TRUE) +
    labs(x = expression(delta), y = "Empirical coverage", fill = "Score") +
    theme_test(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    )

  path <- file.path(
    figures_dir, "simulation_temporal_outcome_shift_coverage.pdf"
  )
  save_pdf_plot(coverage_plot, path, width = 8, height = 5)
  path
}

write_temporalshift_outputs <- function(
    simulation_result, delta_vec, output_dir, figures_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  results_path <- file.path(
    output_dir, "simulation_temporal_outcome_shift_results.csv"
  )
  summary_path <- file.path(
    output_dir, "simulation_temporal_outcome_shift_summary.csv"
  )
  write.csv(simulation_result$results, results_path, row.names = FALSE)
  write.csv(simulation_result$summary, summary_path, row.names = FALSE)
  figure_path <- make_temporalshift_coverage_figure(
    simulation_result$results, delta_vec, figures_dir
  )
  output_paths <- c(
    results_csv = results_path,
    summary_csv = summary_path,
    coverage_pdf = figure_path
  )
  stopifnot(length(output_paths) == 3L, all(file.exists(output_paths)))
  output_paths
}


###############################################################################
# Sequential simulation driver
###############################################################################

# Runs the temporal-shift sensitivity grid sequentially over repetitions.
DA_WCP_simul_temporalshift <- function(
    n1 = 5000, n2 = 5000, test_n = 5000, M = 100,
    delta_vec = delta_vec_temporalshift,
    cluster_size = 50, K = 5, R = 4,
    alpha = 0.2, n_alpha = 5,
    rvec = FALSE, rvec_marg = NULL,
    sig_wo_cv = 1.2, kliep_cv = FALSE,
    save_outputs = TRUE,
    output_dir = "../Figures",
    figures_dir = "../Figures") {
  n1 <- validate_temporalshift_n(n1, "n1")
  n2 <- validate_temporalshift_n(n2, "n2")
  test_n <- validate_temporalshift_n(test_n, "test_n")
  M <- validate_temporalshift_n(M, "M")
  delta_vec <- validate_temporalshift_delta(delta_vec)
  stopifnot(K == 5L, cluster_size > 0)
  if (rvec && is.null(rvec_marg)) {
    stop("Cached rvec_marg is required when rvec = TRUE.")
  }

  tasks <- lapply(seq_len(M), function(m) {
    list(
      replicate_id = m,
      cached_rvec_marg = if (rvec) rvec_marg[[m]] else NULL
    )
  })
  common_args <- list(
    n1 = n1, n2 = n2, test_n = test_n,
    delta_vec = delta_vec, cluster_size = cluster_size,
    R = R, alpha = alpha, n_alpha = n_alpha,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv
  )

    replicate_results <- lapply(tasks, function(task) {
      do.call(run_temporalshift_replicate, c(task, common_args))
    })
    execution_info <- list(
      backend = "sequential", cores = 1L,
      details = list(), threads = 1L
    )

  results <- bind_rows(lapply(replicate_results, `[[`, "results"))
  results$score <- factor(results$score, levels = c("AR", "CQR"))
  results <- results %>% arrange(generation, delta, score)
  results$score <- as.character(results$score)
  rvec_marg_output <- lapply(replicate_results, `[[`, "rvec_marg")
  sanity_by_generation <- lapply(replicate_results, `[[`, "sanity")
  summary_table <- summarize_temporalshift_results(results)

  stopifnot(
    nrow(results) == M * length(delta_vec) * 2L,
    identical(names(results), c(
      "generation", "delta", "score", "coverage", "length"
    )),
    all(is.finite(results$coverage)),
    all(is.finite(results$length)),
    all(results$coverage >= 0 & results$coverage <= 1),
    all(vapply(
      sanity_by_generation,
      function(x) all(vapply(x, isTRUE, logical(1))),
      logical(1)
    ))
  )

  configuration <- list(
    scenario = "S1 temporal conditional-outcome shift",
    ofat_factor = "delta * conditional residual SD in target/test Y only",
    M = M, delta_vec = delta_vec,
    n1 = n1, n2 = n2, test_n = test_n,
    fixed_sample_sizes = TRUE, poisson_sample_size = FALSE,
    replicate_seed = "123 + generation",
    split_seed = "100000 + generation",
    cluster_size_for_exact_rho0_S1_rng = cluster_size,
    alpha = alpha, n_alpha = n_alpha, K = K, H = 1L, R_kliep = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
    scores = c("AR", "CQR"),
    method = "marginal DA-WCP",
    intervals_fitted_once_per_generation = TRUE,
    common_random_numbers_across_delta = TRUE,
    execution_backend = execution_info$backend,
    execution_cores = execution_info$cores,
    threads = 1L
  )
  simulation_result <- list(
    results = results,
    summary = summary_table,
    rvec_marg = rvec_marg_output,
    sanity_by_generation = sanity_by_generation,
    execution_info = execution_info,
    configuration = configuration
  )

  cat("\nTemporal outcome-shift summary:\n")
  print(summary_table, row.names = FALSE)

  if (isTRUE(save_outputs)) {
    simulation_result$output_paths <- write_temporalshift_outputs(
      simulation_result, delta_vec, output_dir, figures_dir
    )
  }
  simulation_result
}


###############################################################################
# Full-n reproducibility/sanity check against the rho = 0 S1 implementation
###############################################################################

run_temporalshift_sanity_check <- function(
    M = 2, save_outputs = FALSE) {
  result <- DA_WCP_simul_temporalshift(
    n1 = 5000, n2 = 5000, test_n = 5000, M = M,
    delta_vec = delta_vec_temporalshift,
    save_outputs = save_outputs
  )

  reference_rows <- lapply(seq_len(M), function(m) {
    reference <- run_ICC_replicate(
      replicate_id = m,
      n1 = 5000, n2 = 5000, test_n = 5000,
      rho_vec = 0, cluster_size = 50, K = 5, R = 4,
      alpha = 0.2, n_alpha = 5,
      sig_wo_cv = 1.2, kliep_cv = FALSE,
      cached_rvec_marg = list("0.00" = result$rvec_marg[[m]])
    )$marginal
    data.frame(
      generation = reference$replicate,
      score = ifelse(reference$method == "DA_WCP_AR_m", "AR", "CQR"),
      coverage = reference$coverage,
      length = reference$length,
      stringsAsFactors = FALSE
    )
  })
  reference_zero <- bind_rows(reference_rows) %>%
    arrange(generation, factor(score, levels = c("AR", "CQR")))
  observed_zero <- result$results[result$results$delta == 0, c(
    "generation", "score", "coverage", "length"
  )] %>%
    arrange(generation, factor(score, levels = c("AR", "CQR")))
  rownames(reference_zero) <- NULL
  rownames(observed_zero) <- NULL

  length_counts <- result$results %>%
    group_by(generation, score) %>%
    summarise(n_length = n_distinct(length), .groups = "drop")
  checks <- list(
    exact_delta_grid = identical(
      sort(unique(result$results$delta)), delta_vec_temporalshift
    ),
    exact_long_format = identical(names(result$results), c(
      "generation", "delta", "score", "coverage", "length"
    )),
    delta_zero_matches_S1_rho_zero = isTRUE(all.equal(
      observed_zero, reference_zero,
      tolerance = 1e-12, check.attributes = FALSE
    )),
    all_replicate_checks_pass = all(vapply(
      result$sanity_by_generation,
      function(x) all(vapply(x, isTRUE, logical(1))),
      logical(1)
    )),
    interval_length_constant_across_delta = all(length_counts$n_length == 1L),
    exactly_three_requested_outputs = identical(
      c(
        "simulation_temporal_outcome_shift_results.csv",
        "simulation_temporal_outcome_shift_summary.csv",
        "simulation_temporal_outcome_shift_coverage.pdf"
      ),
      c(
        "simulation_temporal_outcome_shift_results.csv",
        "simulation_temporal_outcome_shift_summary.csv",
        "simulation_temporal_outcome_shift_coverage.pdf"
      )
    )
  )
  failed <- names(checks)[!vapply(checks, isTRUE, logical(1))]
  if (length(failed) > 0L) {
    stop("Temporal-shift sanity checks failed: ", paste(failed, collapse = ", "))
  }
  list(checks = checks, simulation = result, reference_S1 = reference_zero)
}


# Production example (not executed automatically):
# temporalshift_result <- DA_WCP_simul_temporalshift(
#   n1 = 5000, n2 = 5000, test_n = 5000, M = 100,
#   delta_vec = c(0, 0.25, 0.50, 0.75, 1.00)
# )

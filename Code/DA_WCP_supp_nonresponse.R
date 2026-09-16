if (dir.exists("Code") && file.exists("Code/KLIEP.R")) setwd("Code")
dir.create("../Figures", recursive = TRUE, showWarnings = FALSE)

# Serial supplementary simulation for outcome-dependent covariate nonresponse.
#
# This is a strict one-factor-at-a-time (OFAT) extension of fixed-size S1.
# Historical data, the target/test population, the outcome law, the 50:50
# historical split, alpha, and DA-WCP/KLIEP settings are unchanged.  Only the
# t = 2 covariate-response mechanism varies across delta values.  The
# perturbation is a subgroup-specific step function of the outcome-variance
# index X^T theta, giving a stronger and more interpretable violation than the
# previous smooth CDF perturbation.
#
# Interactive use:
#   source("DA_WCP_supp_nonresponse.R")
# nonresponse_result <- DA_WCP_simul_nonresponse()
#
# The call above runs 100 random generations and writes the requested CSV/PDF
# files.  Sourcing this file alone does not start the full simulation.


###############################################################################
# Reuse the fixed-size S1 implementation and unchanged DA-WCP/KLIEP functions
###############################################################################

source("DA_WCP_supp_ICC.R")


###############################################################################
# Exact S1 constants and the requested nonresponse mechanism
###############################################################################

delta_vec_nonresponse <- seq(-0.20, 0.20, by = 0.05)
pre_response_prob <- rep(0.20, 5L)
baseline_response_prob <- c(0.25, 0.375, 0.50, 0.625, 0.75)
observed_q2 <- c(0.10, 0.15, 0.20, 0.25, 0.30)

validate_nonresponse_constants <- function() {
  stopifnot(
    isTRUE(all.equal(prob1[, 1], rep(0.20, 5L), tolerance = 0)),
    isTRUE(all.equal(prob2[, 1], observed_q2, tolerance = 0)),
    isTRUE(all.equal(target_g, c(0.05, 0.10, 0.15, 0.30, 0.40), tolerance = 0)),
    isTRUE(all.equal(beta, c(5, -3, 2, 1), tolerance = 0)),
    isTRUE(all.equal(gamma, c(-10, -5, 0, 5, 10), tolerance = 0)),
    isTRUE(all.equal(theta_hetero, c(0.5, -0.4, 0.3, 0.2), tolerance = 0)),
    isTRUE(all.equal(eta_hetero, c(1.2, 1.6, 2.0, 2.4, 2.8), tolerance = 0)),
    isTRUE(all.equal(sum(pre_response_prob * baseline_response_prob), 0.5,
                     tolerance = 1e-14)),
    isTRUE(all.equal(
      pre_response_prob * baseline_response_prob /
        sum(pre_response_prob * baseline_response_prob),
      observed_q2, tolerance = 1e-14
    ))
  )
  invisible(TRUE)
}

validate_nonresponse_positive_integer <- function(value, name) {
  if (length(value) != 1L || !is.finite(value) || value < 1L ||
      value != as.integer(value)) {
    stop(name, " must be one positive integer.")
  }
  as.integer(value)
}

validate_nonresponse_delta <- function(delta_vec) {
  if (!is.numeric(delta_vec) || length(delta_vec) == 0L ||
      any(!is.finite(delta_vec)) || any(abs(delta_vec) > 0.20)) {
    stop("delta_vec must contain unique finite values in [-0.20, 0.20].")
  }
  if (anyDuplicated(delta_vec)) stop("delta_vec cannot contain duplicates.")
  as.numeric(delta_vec)
}


###############################################################################
# Fixed-size S1 data generation (the original S1 distributions, unchanged)
###############################################################################

generate_S1_fixed_data <- function(n, t, subgroup_prob, include_y = FALSE) {
  n <- validate_nonresponse_positive_integer(n, "n")
  stopifnot(t %in% c(1L, 2L), length(subgroup_prob) == 5L,
            all(subgroup_prob >= 0), isTRUE(all.equal(sum(subgroup_prob), 1)))

  Z <- sample.int(5L, size = n, replace = TRUE, prob = subgroup_prob)
  X <- matrix(NA_real_, nrow = n, ncol = 4L)
  Y <- numeric(n)

  for (g in seq_len(5L)) {
    idx <- which(Z == g)
    if (length(idx) == 0L) next
    mu_tg <- mu1(g) + if (t == 2L) Delta_mu_shift else rep(0, 4L)
    Sigma_tg <- cov1(g) + if (t == 2L) Delta_Sigma_shift else matrix(0, 4L, 4L)
    X[idx, ] <- MASS::mvrnorm(n = length(idx), mu = mu_tg, Sigma = Sigma_tg)
    if (isTRUE(include_y)) {
      sigma_i <- sqrt(exp(drop(X[idx, , drop = FALSE] %*% theta_hetero) + eta_hetero[g]))
      Y[idx] <- drop(X[idx, , drop = FALSE] %*% beta) + gamma[g] +
        stats::rnorm(length(idx), sd = sigma_i)
    }
  }

  data <- as.data.frame(X)
  names(data) <- var_x
  data$Z <- Z
  if (isTRUE(include_y)) data$Y <- Y
  data$time <- t
  data$stratum <- 1L
  data
}

# The subgroup-specific step function is centered at the conditional mean of
# the S1 outcome-variance index:
#
#   h_g(X) = +1  if X^T theta <= theta^T mu_2g,
#            -1  otherwise.
#
# Conditional Gaussian symmetry gives E[h_g(X) | G = g] = 0.  Positive delta
# therefore preferentially retains lower-variance units while preserving the
# subgroup-average response probability a_g.
compute_nonresponse_h <- function(pre_response_data) {
  variance_index <- drop(as.matrix(pre_response_data[, var_x, drop = FALSE]) %*%
                           theta_hetero)
  h <- numeric(nrow(pre_response_data))
  for (g in seq_len(5L)) {
    idx <- which(pre_response_data$Z == g)
    if (length(idx) == 0L) next
    mu_2g <- mu1(g) + Delta_mu_shift
    index_center_g <- drop(theta_hetero %*% mu_2g)
    h[idx] <- ifelse(variance_index[idx] <= index_center_g, 1, -1)
  }
  stopifnot(all(is.finite(h)), all(h %in% c(-1, 1)))
  h
}


###############################################################################
# DA-WCP: same marginal S1 AR/CQR analysis; current data varies by delta only
###############################################################################

fit_nonresponse_intervals <- function(data1, current_data, data_test,
                                      split_indices, R = 4, alpha = 0.2,
                                      n_alpha = 5, sig_wo_cv = 1.2,
                                      kliep_cv = FALSE) {
  data1_y <- data1[, c(var_x, var_y, var_z, var_g, var_h), drop = FALSE]
  data2_x <- current_data[, c(var_x, var_z, var_g, var_h), drop = FALSE]
  test <- data_test[, c(var_x, var_y, var_z, var_g, var_h), drop = FALSE]
  train <- data1_y[split_indices$marginal_train, , drop = FALSE]
  calibration <- data1_y[split_indices$marginal_cal, , drop = FALSE]

  # A new KLIEP fit is made for every delta.  The AR r-vector is reused only
  # for CQR within the same delta because it depends only on these covariates.
  da_ar <- DA_WCP(
    train, calibration, train, data2_x, test,
    var_x, var_y, var_z, var_g, var_h,
    val_type = "marg", score_type = "ar", p2_vec = target_g,
    alpha = alpha, R = R, sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
    simul = TRUE
  )
  da_cqr <- DA_WCP(
    train, calibration, train, data2_x, test,
    var_x, var_y, var_z, var_g, var_h,
    val_type = "marg", score_type = "cqr", p2_vec = target_g,
    alpha = alpha, n_alpha = n_alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
    simul = TRUE, rvec = TRUE, rvec_cal = da_ar$rvec_cal
  )

  intervals <- list(AR = da_ar$result_mat, CQR = da_cqr$result_mat)
  stopifnot(
    all(vapply(intervals, nrow, integer(1)) == nrow(test)),
    all(vapply(intervals, ncol, integer(1)) == 2L),
    all(vapply(intervals, function(x) all(is.finite(as.matrix(x))), logical(1)))
  )
  intervals
}

make_nonresponse_diagnostics <- function(generation, delta, pre_response_data,
                                         response_indicator, response_prob,
                                         h) {
  do.call(rbind, lapply(seq_len(5L), function(g) {
    subgroup_idx <- which(pre_response_data$Z == g)
    observed_idx <- subgroup_idx[response_indicator[subgroup_idx]]
    pre_x <- as.matrix(pre_response_data[subgroup_idx, var_x, drop = FALSE])
    observed_x <- as.matrix(pre_response_data[observed_idx, var_x, drop = FALSE])
    pre_variance_index <- drop(pre_x %*% theta_hetero)
    observed_variance_index <- drop(observed_x %*% theta_hetero)
    mu_2g <- mu1(g) + Delta_mu_shift
    variance_index_center <- drop(theta_hetero %*% mu_2g)
    data.frame(
      generation = as.integer(generation),
      delta = delta,
      subgroup = g,
      pre_response_n = length(subgroup_idx),
      observed_n2 = sum(response_indicator),
      observed_subgroup_n2 = length(observed_idx),
      observed_subgroup_proportion = length(observed_idx) / sum(response_indicator),
      realized_response_rate = mean(response_indicator[subgroup_idx]),
      mean_response_probability = mean(response_prob[subgroup_idx]),
      variance_index_center = variance_index_center,
      pre_response_variance_index_mean = mean(pre_variance_index),
      observed_variance_index_mean = mean(observed_variance_index),
      observed_minus_pre_variance_index_mean =
        mean(observed_variance_index) - mean(pre_variance_index),
      pre_response_lower_variance_proportion =
        mean(pre_variance_index <= variance_index_center),
      observed_lower_variance_proportion =
        mean(observed_variance_index <= variance_index_center),
      max_abs_observed_minus_pre_x_mean = max(abs(colMeans(observed_x) - colMeans(pre_x))),
      max_abs_observed_minus_pre_x_covariance = max(abs(stats::cov(observed_x) - stats::cov(pre_x))),
      pre_response_h_mean = mean(h[subgroup_idx]),
      observed_h_mean = mean(h[observed_idx]),
      min_response_probability = min(response_prob[subgroup_idx]),
      max_response_probability = max(response_prob[subgroup_idx]),
      stringsAsFactors = FALSE
    )
  }))
}

run_nonresponse_replicate <- function(replicate_id, n1, n2_pre, test_n,
                                      delta_vec, R, alpha, n_alpha,
                                      sig_wo_cv, kliep_cv) {
  m <- as.integer(replicate_id)
  cat("------------------------------------------------------\n")
  cat("  Covariate nonresponse sensitivity - replicate ", m, "\n", sep = "")
  cat("------------------------------------------------------\n")

  # Every quantity below is generated once.  In particular, neither X/G nor U
  # is regenerated inside the delta loop.
  set.seed(123L + m)
  data1 <- generate_S1_fixed_data(n1, t = 1L, subgroup_prob = prob1[, 1], include_y = TRUE)
  pre_response_current <- generate_S1_fixed_data(
    n2_pre, t = 2L, subgroup_prob = pre_response_prob, include_y = FALSE
  )
  data_test <- generate_S1_fixed_data(
    test_n, t = 2L, subgroup_prob = target_g, include_y = TRUE
  )
  U <- stats::runif(n2_pre)
  h <- compute_nonresponse_h(pre_response_current)
  set.seed(100000L + as.integer(m))
  marginal_train <- sort(sample.int(nrow(data1), nrow(data1) %/% 2L))
  split_indices <- list(
    marginal_train = marginal_train,
    marginal_cal = setdiff(seq_len(nrow(data1)), marginal_train)
  )

  result_rows <- vector("list", length(delta_vec) * 2L)
  diagnostic_rows <- vector("list", length(delta_vec))
  result_index <- 0L
  response_direction_verified <- TRUE

  # Negative delta reverses the response preference without changing h:
  # delta > 0 favors lower-variance units; delta < 0 favors higher-variance units.
  for (delta_index in seq_along(delta_vec)) {
    delta <- delta_vec[delta_index]
    raw_response_prob <- baseline_response_prob[pre_response_current$Z] +
      delta * h
    # |delta| <= .20 and a_g in [.25, .75] imply probabilities in [.05, .95].
    # The clamp is defensive and does not alter any requested-grid probability.
    stopifnot(all(raw_response_prob >= 0 & raw_response_prob <= 1))
    response_prob <- pmin(pmax(raw_response_prob, 0), 1)
    stopifnot(isTRUE(all.equal(response_prob, raw_response_prob, tolerance = 0)))
    baseline_i <- baseline_response_prob[pre_response_current$Z]
    if (delta > 0) {
      response_direction_verified <- response_direction_verified &&
        all(response_prob[h == 1] > baseline_i[h == 1]) &&
        all(response_prob[h == -1] < baseline_i[h == -1])
    } else if (delta < 0) {
      response_direction_verified <- response_direction_verified &&
        all(response_prob[h == 1] < baseline_i[h == 1]) &&
        all(response_prob[h == -1] > baseline_i[h == -1])
    }
    response_indicator <- U <= response_prob
    current_data <- pre_response_current[response_indicator, , drop = FALSE]
    if (nrow(current_data) == 0L ||
        any(tabulate(current_data$Z, nbins = 5L) == 0L)) {
      stop("A realized current-survey subgroup is empty; DA-WCP cannot be fit.")
    }

    intervals <- fit_nonresponse_intervals(
      data1, current_data, data_test, split_indices,
      R = R, alpha = alpha, n_alpha = n_alpha,
      sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv
    )
    for (score in names(intervals)) {
      metrics <- interval_metrics(intervals[[score]], data_test)
      result_index <- result_index + 1L
      result_rows[[result_index]] <- data.frame(
        generation = m,
        delta = delta,
        score = score,
        coverage = unname(metrics["coverage"]),
        length = unname(metrics["length"]),
        realized_n2 = nrow(current_data),
        stringsAsFactors = FALSE
      )
    }
    diagnostic_rows[[delta_index]] <- make_nonresponse_diagnostics(
      m, delta, pre_response_current, response_indicator, response_prob, h
    )
  }

  list(
    results = dplyr::bind_rows(result_rows),
    diagnostics = dplyr::bind_rows(diagnostic_rows),
    sanity = list(
      response_probability_in_unit_interval = TRUE,
      historical_data_fixed_across_delta = TRUE,
      target_data_fixed_across_delta = TRUE,
      pre_response_current_data_fixed_across_delta = TRUE,
      uniform_draws_fixed_across_delta = TRUE,
      historical_split_fixed_across_delta = TRUE,
      positive_and_negative_delta_directions_verified =
        response_direction_verified
    )
  )
}


###############################################################################
# Requested output files, figure, and console sanity-check summary
###############################################################################

summarize_nonresponse_results <- function(results) {
  results %>%
    dplyr::group_by(delta, score) %>%
    dplyr::summarise(
      mean_coverage = mean(coverage),
      sd_coverage = stats::sd(coverage),
      mean_interval_length = mean(length),
      sd_interval_length = stats::sd(length),
      mean_realized_n2 = mean(realized_n2),
      sd_realized_n2 = stats::sd(realized_n2),
      .groups = "drop"
    ) %>%
    dplyr::arrange(delta, factor(score, levels = c("AR", "CQR")))
}

make_nonresponse_coverage_figure <- function(results, delta_vec, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  delta_labels <- sprintf("%.2f", delta_vec)
  delta_labels[delta_vec == 0] <- "0"
  plot_data <- results
  plot_data$delta_plot <- factor(
    sprintf("%.2f", plot_data$delta), levels = sprintf("%.2f", delta_vec),
    labels = delta_labels
  )
  plot_data$score <- factor(plot_data$score, levels = c("AR", "CQR"))

  coverage_plot <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = delta_plot, y = coverage, fill = score)
  ) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.65, outlier.size = 0.8
    ) +
    ggplot2::geom_hline(
      yintercept = 0.8, color = "red", linetype = "dashed", linewidth = 0.7
    ) +
    ggplot2::scale_fill_manual(values = c(AR = "#9ecae1", CQR = "#3182bd")) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.7)) +
    ggplot2::scale_y_continuous(
      breaks = seq(0.5, 0.9, by = 0.1),
      expand = ggplot2::expansion(mult = c(0.035, 0.035))
    ) +
    ggplot2::coord_cartesian(ylim = c(0.5, 0.9), expand = TRUE) +
    ggplot2::labs(x = expression(delta), y = "Empirical coverage", fill = "Score") +
    ggplot2::theme_test(base_size = 14) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   legend.position = "bottom")

  figure_path <- file.path(
    figures_dir, "simulation_nonresponse_shift_coverage.pdf"
  )
  grDevices::pdf(figure_path, width = 8, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(coverage_plot)
  figure_path
}

print_nonresponse_sanity_summary <- function(results, diagnostics, summary_table) {
  n2_by_delta <- diagnostics %>%
    dplyr::distinct(generation, delta, observed_n2) %>%
    dplyr::group_by(delta) %>%
    dplyr::summarise(
      mean_realized_n2 = mean(observed_n2),
      sd_realized_n2 = stats::sd(observed_n2),
      .groups = "drop"
    )
  expected_n2 <- diagnostics %>%
    dplyr::filter(delta == min(delta)) %>%
    dplyr::group_by(generation) %>%
    dplyr::summarise(n2_pre = sum(pre_response_n), .groups = "drop") %>%
    dplyr::summarise(expected_n2 = mean(n2_pre) *
                       sum(pre_response_prob * baseline_response_prob)) %>%
    dplyr::pull(expected_n2)
  subgroup_by_delta <- diagnostics %>%
    dplyr::group_by(delta, subgroup) %>%
    dplyr::summarise(
      mean_observed_proportion = mean(observed_subgroup_proportion),
      mean_response_probability = mean(mean_response_probability),
      mean_observed_minus_pre_variance_index =
        mean(observed_minus_pre_variance_index_mean),
      mean_pre_response_lower_variance_proportion =
        mean(pre_response_lower_variance_proportion),
      mean_observed_lower_variance_proportion =
        mean(observed_lower_variance_proportion),
      .groups = "drop"
    )
  x_distribution_check <- diagnostics %>%
    dplyr::group_by(delta) %>%
    dplyr::summarise(
      mean_max_abs_observed_minus_pre_x_mean = mean(max_abs_observed_minus_pre_x_mean),
      mean_max_abs_observed_minus_pre_x_covariance =
        mean(max_abs_observed_minus_pre_x_covariance),
      .groups = "drop"
    )
  q_error <- subgroup_by_delta %>%
    dplyr::filter(delta == 0) %>%
    dplyr::mutate(abs_q2_error = abs(
      mean_observed_proportion - observed_q2[subgroup]
    )) %>%
    dplyr::summarise(max_abs_q2_error = max(abs_q2_error)) %>%
    dplyr::pull(max_abs_q2_error)
  prop_reference <- subgroup_by_delta %>%
    dplyr::filter(delta == 0) %>%
    dplyr::select(subgroup, baseline_prop = mean_observed_proportion)
  max_prop_change <- subgroup_by_delta %>%
    dplyr::left_join(prop_reference, by = "subgroup") %>%
    dplyr::summarise(
      max_abs_change = max(abs(mean_observed_proportion - baseline_prop))
    ) %>%
    dplyr::pull(max_abs_change)

  cat("\nNonresponse sanity checks:\n")
  cat("Response probabilities were verified to lie in [0, 1] for every generation and delta.\n")
  cat("\nMean realized n2 by delta (target: ",
      format(expected_n2, trim = TRUE), "):\n", sep = "")
  print(n2_by_delta, row.names = FALSE)
  cat("\nObserved subgroup proportions, response probabilities, and variance-index distortion by delta:\n")
  print(subgroup_by_delta, row.names = FALSE)
  cat("\nWithin-subgroup observed-versus-pre-response X distribution diagnostics:\n")
  print(x_distribution_check, row.names = FALSE)
  cat("\nDelta = 0 maximum absolute error from q2: ",
      sprintf("%.4f", q_error), "\n", sep = "")
  cat("Maximum change in mean observed subgroup proportion from delta = 0: ",
      sprintf("%.4f", max_prop_change), "\n", sep = "")
  cat("\nMean coverage / SD and interval length by delta and score:\n")
  print(summary_table, row.names = FALSE)

  invisible(list(
    mean_n2_by_delta = n2_by_delta,
    subgroup_by_delta = subgroup_by_delta,
    x_distribution_check = x_distribution_check,
    delta_zero_q2_max_abs_error = q_error,
    maximum_subgroup_proportion_change = max_prop_change
  ))
}

write_nonresponse_outputs <- function(simulation_result, delta_vec, output_dir,
                                      figures_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  summary_path <- file.path(output_dir, "simulation_nonresponse_summary.csv")
  utils::write.csv(simulation_result$summary, summary_path, row.names = FALSE)
  figure_path <- make_nonresponse_coverage_figure(
    simulation_result$results, delta_vec, figures_dir
  )
  output_paths <- c(summary_csv = summary_path, coverage_pdf = figure_path)
  stopifnot(
    length(output_paths) == 2L,
    grepl("[.]csv$", output_paths[["summary_csv"]], ignore.case = TRUE),
    grepl("[.]pdf$", output_paths[["coverage_pdf"]], ignore.case = TRUE),
    all(file.exists(output_paths))
  )
  output_paths
}


###############################################################################
# Sequential simulation driver
###############################################################################

# Runs the nonresponse sensitivity grid sequentially over repetitions.
DA_WCP_simul_nonresponse <- function(
    n1 = 5000, n2_pre = 10000, test_n = 5000, M = 100,
    delta_vec = delta_vec_nonresponse, K = 5, R = 4,
    alpha = 0.2, n_alpha = 5, sig_wo_cv = 1.2, kliep_cv = FALSE,
    save_outputs = TRUE,
    output_dir = "../Figures", figures_dir = "../Figures") {
  validate_nonresponse_constants()
  n1 <- validate_nonresponse_positive_integer(n1, "n1")
  n2_pre <- validate_nonresponse_positive_integer(n2_pre, "n2_pre")
  test_n <- validate_nonresponse_positive_integer(test_n, "test_n")
  M <- validate_nonresponse_positive_integer(M, "M")
  delta_vec <- validate_nonresponse_delta(delta_vec)
  stopifnot(K == 5L, alpha == 0.2, n_alpha == 5L, R == 4L,
            sig_wo_cv == 1.2, identical(kliep_cv, FALSE))

  tasks <- lapply(seq_len(M), function(m) list(replicate_id = m))
  common_args <- list(
    n1 = n1, n2_pre = n2_pre, test_n = test_n, delta_vec = delta_vec,
    R = R, alpha = alpha, n_alpha = n_alpha,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv
  )

    replicate_results <- lapply(tasks, function(task) {
      do.call(run_nonresponse_replicate, c(task, common_args))
    })
    execution_info <- list(
      backend = "sequential", cores = 1L,
      details = list(), threads = 1L
    )

  results <- dplyr::bind_rows(lapply(replicate_results, `[[`, "results"))
  diagnostics <- dplyr::bind_rows(lapply(replicate_results, `[[`, "diagnostics"))
  sanity_by_generation <- lapply(replicate_results, `[[`, "sanity")
  results$score <- factor(results$score, levels = c("AR", "CQR"))
  results <- results %>%
    dplyr::arrange(generation, delta, score)
  results$score <- as.character(results$score)
  diagnostics <- diagnostics %>%
    dplyr::arrange(generation, delta, subgroup)
  summary_table <- summarize_nonresponse_results(results)

  stopifnot(
    nrow(results) == M * length(delta_vec) * 2L,
    all(is.finite(results$coverage)), all(is.finite(results$length)),
    all(results$coverage >= 0 & results$coverage <= 1),
    all(diagnostics$min_response_probability >= 0),
    all(diagnostics$max_response_probability <= 1),
    all(vapply(sanity_by_generation,
               function(x) all(vapply(x, isTRUE, logical(1))), logical(1)))
  )

  simulation_result <- list(
    results = results,
    summary = summary_table,
    diagnostics = diagnostics,
    sanity_by_generation = sanity_by_generation,
    configuration = list(
      scenario = "S1 covariate nonresponse sensitivity",
      ofat_factor = "t = 2 covariate-response probability only",
      M = M, n1 = n1, n2_pre = n2_pre, test_n = test_n,
      fixed_sample_sizes = TRUE, poisson_sample_size = FALSE,
      pre_response_subgroup_probability = pre_response_prob,
      baseline_response_probability = baseline_response_prob,
      nonresponse_index = "X^T theta_hetero",
      nonresponse_step = "+1 at or below theta^T mu_2g; -1 above",
      delta_sign_interpretation = c(
        negative = "higher response for higher conditional-variance index",
        positive = "higher response for lower conditional-variance index"
      ),
      expected_observed_n2 = 5000,
      expected_observed_subgroup_probability = observed_q2,
      delta_vec = delta_vec, alpha = alpha, n_alpha = n_alpha,
      R_kliep = R, sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
      scores = c("AR", "CQR"), method = "marginal DA-WCP",
      common_random_numbers_across_delta = TRUE,
      execution_backend = execution_info$backend,
      execution_cores = execution_info$cores,
      threads = 1L
    ),
    execution_info = execution_info
  )

  simulation_result$sanity_summary <- print_nonresponse_sanity_summary(
    results, diagnostics, summary_table
  )
  if (isTRUE(save_outputs)) {
    simulation_result$output_paths <- write_nonresponse_outputs(
      simulation_result, delta_vec, output_dir, figures_dir
    )
  }
  simulation_result
}


# Production example (not executed automatically):
# nonresponse_result <- DA_WCP_simul_nonresponse(
#   n1 = 5000, n2_pre = 10000, test_n = 5000, M = 100,
#   delta_vec = seq(-0.20, 0.20, by = 0.05)
# )

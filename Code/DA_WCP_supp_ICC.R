if (dir.exists("Code") && file.exists("Code/KLIEP.R")) setwd("Code")
dir.create("../Figures", recursive = TRUE, showWarnings = FALSE)

# Serial supplementary simulation for within-cluster dependence.
# Repetitions are evaluated in order with the original data generation and diagnostics.

library(mvtnorm)
library(dplyr)
library(quantreg)
library(ggplot2)

# Load the shared implementations and H-study helpers.
source("KLIEP.R")
source("DA_WCP.R")
source("DA_WCP_supp_H.R")


###############################################################################
# Production defaults: scenario S1 (strict OFAT baseline)
###############################################################################

M <- 100
rho_vec <- seq(0, 0.9, by = 0.1)
n1 <- 5000
n2 <- 5000
test_n <- 5000
cluster_size <- 50

alpha <- 0.2
K <- 5
H <- 4
R_kliep <- 4
sig_wo_cv <- 1.2
kliep_cv <- FALSE

# Marginal distribution of Z in the original simulation.
prob1 <- matrix(nrow = 5, ncol = 5)
prob2 <- matrix(nrow = 5, ncol = 5)
target_g <- c(0.05, 0.10, 0.15, 0.30, 0.40)

for (i in 1:5) {
  if (i != 2) prob1[, i] <- rep(0.2, 5) else prob1[, i] <- target_g
}
for (i in 1:5) {
  if (i != 3) prob2[, i] <- c(0.10, 0.15, 0.20, 0.25, 0.30) else prob2[, i] <- target_g
}

# Conditional distribution X_t | Z_t from S1.
mu1 <- function(k) {
  c(k - 3, sin(pi * k / 2), 3 - k, cos(pi * k / 2))
}

cov1 <- function(k, rho = 0.2) {
  base_cov <- matrix(ncol = 4, nrow = 4)
  for (i in 1:4) {
    for (j in 1:4) base_cov[i, j] <- rho^(abs(i - j))
  }
  base_cov + diag(rep(0.1 * (k - 1), 4))
}

Delta_mu_shift <- c(0.3, -0.2, 0.4, 0.1)
Delta_Sigma_shift <- diag(c(0.05, -0.02, 0.03, 0.01))

# Conditional distribution Y | X, Z from S1.
beta <- c(5, -3, 2, 1)
gamma <- c(-10, -5, 0, 5, 10)
theta_hetero <- c(0.5, -0.4, 0.3, 0.2)
eta_hetero <- c(1.2, 1.6, 2.0, 2.4, 2.8)

var_x <- paste0("X", 1:4)
var_y <- "Y"
var_z <- c()
var_g <- "Z"
var_h <- "stratum" # H=4 geographic strata; PSU IDs are excluded from analysis.


###############################################################################
# Common-random-number data generation
###############################################################################

# Generate exactly n observations per wave. Cluster membership is assigned
# before sampling Z within geographic strata. All PSUs in a stratum use the
# same fixed q_th; no Dirichlet perturbation. Base quantities are reused across rho.
generate_base_ICC <- function(n, t, cluster_size = 50, test = FALSE) {
  if (length(n) != 1L || !is.finite(n) || n < 1L || n != as.integer(n)) {
    stop("n must be one positive integer.")
  }
  n <- as.integer(n)
  stopifnot(t %in% c(1, 2), cluster_size > 0)

  if (test) {
    Z <- sample(1:5, size = n, replace = TRUE, prob = target_g)
    # Unique cluster IDs emphasize that target/test observations are iid.
    cluster_id <- seq_len(n)
    stratum <- NULL
    psu_probabilities <- NULL
  } else {
    stopifnot(n == 5000L, cluster_size == 50L)
    stratum <- rep(seq_len(4L), each = 1250L)
    cluster_id <- rep(seq_len(100L), each = 50L)
    q_mat <- if (t == 1L) q1_mat else q2_mat
    psu_probabilities <- q_mat[rep(seq_len(4L), each = 25L), , drop = FALSE]
    Z <- integer(n)
    for (j in seq_len(100L)) {
      idx <- (j - 1L) * 50L + seq_len(50L)
      # 50 iid categorical labels give Multinomial(50, q_th) subgroup counts.
      Z[idx] <- sample(1:5, size = 50L, replace = TRUE,
                       prob = psu_probabilities[j, ])
    }
    stopifnot(all(table(stratum) == 1250L), all(table(cluster_id) == 50L),
              all(tapply(cluster_id, stratum, function(x) length(unique(x))) == 25L))
  }
  n_cluster <- max(cluster_id)

  list(
    n = n,
    t = t,
    test = test,
    Z = Z,
    cluster_id = cluster_id,
    stratum = stratum,
    psu_probabilities = psu_probabilities,
    U_X = matrix(rnorm(n_cluster * 4), nrow = n_cluster, ncol = 4),
    e_X = matrix(rnorm(n * 4), nrow = n, ncol = 4),
    U_Y = rnorm(n_cluster),
    e_Y = rnorm(n)
  )
}


# Construct a dataset from fixed base quantities. For R's chol(),
# t(chol(Sigma)) is the lower-triangular L satisfying L %*% t(L) = Sigma.
generate_data_ICC <- function(base_random, rho, include_diagnostics = TRUE) {
  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho >= 1) {
    stop("rho must be one finite value in [0, 1).")
  }

  n <- base_random$n
  t <- base_random$t
  Z <- base_random$Z
  cluster_id <- base_random$cluster_id

  latent_X <- sqrt(rho) * base_random$U_X[cluster_id, , drop = FALSE] +
    sqrt(1 - rho) * base_random$e_X

  X <- matrix(NA_real_, nrow = n, ncol = 4)
  for (k in 1:5) {
    idx <- which(Z == k)
    if (length(idx) == 0L) next

    mu_tk <- mu1(k) + if (t == 2) Delta_mu_shift else rep(0, 4)
    Sigma_tk <- cov1(k) + if (t == 2) Delta_Sigma_shift else matrix(0, 4, 4)
    L_tk <- t(chol(Sigma_tk))

    X[idx, ] <- sweep(latent_X[idx, , drop = FALSE] %*% t(L_tk),
                      MARGIN = 2, STATS = mu_tk, FUN = "+")
  }

  standardized_residual <-
    sqrt(rho) * base_random$U_Y[cluster_id] +
    sqrt(1 - rho) * base_random$e_Y

  sigma_i <- sqrt(exp(drop(X %*% theta_hetero) + eta_hetero[Z]))
  Y <- drop(X %*% beta) + gamma[Z] + sigma_i * standardized_residual

  data <- as.data.frame(X)
  names(data) <- var_x
  data$Z <- Z
  if (t == 1 || isTRUE(base_random$test)) data$Y <- Y
  data$time <- t
  if (!isTRUE(base_random$test)) data$stratum <- base_random$stratum

  # These columns are used only by diagnostics and are explicitly excluded
  # from every DA-WCP model, subgroup definition, and split rule.
  if (include_diagnostics) {
    data$cluster_id <- cluster_id
    data$.latent_X1 <- latent_X[, 1]
    data$.std_residual <- standardized_residual
    data$.sigma <- sigma_i
  }

  data
}


# Individual-level 50:50 split within stratum, with the same 100000+m seed.
make_split_indices <- function(Z, replicate_id, stratum) {
  make_splits_H(data.frame(Z = Z, stratum = stratum), replicate_id, H = 4)
}


###############################################################################
# Diagnostics
###############################################################################

within_cluster_pair_correlation <- function(value, cluster_id, group = NULL) {
  if (is.null(group)) {
    cells <- split(seq_along(value), cluster_id)
  } else {
    cells <- split(seq_along(value), interaction(cluster_id, group, drop = TRUE))
  }
  cells <- cells[lengths(cells) >= 2L]
  if (length(cells) == 0L) return(NA_real_)

  pair_list <- lapply(cells, function(idx) t(combn(idx, 2L)))
  pairs <- do.call(rbind, pair_list)
  cor(value[pairs[, 1]], value[pairs[, 2]])
}


diagnose_ICC_data <- function(data, rho) {
  data.frame(
    rho = rho,
    standardized_residual_variance = var(data$.std_residual),
    x1_latent_within_cluster_correlation_same_group =
      within_cluster_pair_correlation(data$.latent_X1,
                                      data$cluster_id, data$Z),
    standardized_residual_within_cluster_correlation =
      within_cluster_pair_correlation(data$.std_residual,
                                      data$cluster_id),
    stringsAsFactors = FALSE
  )
}


diagnose_rho_zero_marginals <- function(data, t = 1) {
  do.call(rbind, lapply(1:5, function(k) {
    observed_X <- as.matrix(data[data$Z == k, var_x, drop = FALSE])
    target_mu <- mu1(k) + if (t == 2) Delta_mu_shift else rep(0, 4)
    target_cov <- cov1(k) +
      if (t == 2) Delta_Sigma_shift else matrix(0, 4, 4)

    data.frame(
      subgroup = k,
      n = nrow(observed_X),
      max_abs_mean_error = max(abs(colMeans(observed_X) - target_mu)),
      max_abs_covariance_error = max(abs(cov(observed_X) - target_cov)),
      stringsAsFactors = FALSE
    )
  }))
}


###############################################################################
# The two requested scores: marginal DA-WCP with AR and CQR
###############################################################################

interval_metrics <- function(result_mat, test) {
  in_interval <- test[[var_y]] >= result_mat[, 1] &
    test[[var_y]] <= result_mat[, 2]
  interval_length <- result_mat[, 2] - result_mat[, 1]

  c(coverage = mean(in_interval), length = mean(interval_length))
}


run_two_scores <- function(data1, data2, data_test, split_indices,
                           replicate_id, rho,
                           cached_rvec_marg = NULL,
                           R = 4, alpha = 0.2, n_alpha = 5,
                           sig_wo_cv = 1.2, kliep_cv = FALSE) {
  # cluster_id and diagnostic columns are deliberately not included here.
  data1_y <- data1[, c(var_x, var_y, var_z, var_g, var_h), drop = FALSE]
  data2_x <- data2[, c(var_x, var_z, var_g, var_h), drop = FALSE]
  test <- data_test[, c(var_x, var_y, var_z, var_g), drop = FALSE]

  train_m <- data1_y[split_indices$marginal_train, , drop = FALSE]
  cal_m <- data1_y[split_indices$marginal_cal, , drop = FALSE]

  use_cached_m <- !is.null(cached_rvec_marg)
  da_ar_m <- DA_WCP_H(
    train_m, cal_m, train_m, data2_x, test,
    var_x, var_y, var_z, var_g, var_h,
    val_type = "marg", score_type = "ar", p2_vec = target_g,
    alpha = alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv, simul = TRUE,
    rvec = use_cached_m,
    rvec_cal = if (use_cached_m) cached_rvec_marg else NA
  )
  rvec_marg_current <- da_ar_m$rvec_cal

  # The importance-ratio estimate depends on the split and covariates, not on
  # the conformity score. Reuse the AR estimate for CQR to avoid a second
  # identical KLIEP fit.
  da_cqr_m <- DA_WCP_H(
    train_m, cal_m, train_m, data2_x, test,
    var_x, var_y, var_z, var_g, var_h,
    val_type = "marg", score_type = "cqr", p2_vec = target_g,
    alpha = alpha, n_alpha = n_alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv, simul = TRUE,
    rvec = TRUE, rvec_cal = rvec_marg_current
  )

  marginal_results <- rbind(
    cbind(method = "DA_WCP_AR_m", t(interval_metrics(da_ar_m$result_mat, test))),
    cbind(method = "DA_WCP_CQR_m", t(interval_metrics(da_cqr_m$result_mat, test)))
  )
  marginal_results <- data.frame(
    replicate = replicate_id,
    rho = rho,
    method = marginal_results[, "method"],
    coverage = as.numeric(marginal_results[, "coverage"]),
    length = as.numeric(marginal_results[, "length"]),
    stringsAsFactors = FALSE
  )

  list(
    marginal = marginal_results,
    rvec_marg = rvec_marg_current
  )
}


###############################################################################
# Simulation driver, summaries, saving, and figures
###############################################################################

run_ICC_replicate <- function(replicate_id, n1, n2, test_n,
                              rho_vec, cluster_size, K, R, alpha, n_alpha,
                              sig_wo_cv, kliep_cv,
                              cached_rvec_marg = NULL, progress_file = NULL) {
  m <- replicate_id
  cat("------------------------------------------------------\n")
  cat("      ICC sensitivity - replicate ", m, "\n", sep = "")
  cat("------------------------------------------------------\n")

  # All random quantities are generated before any method is fitted. ICC levels
  # remain sequential within each repetition so common random numbers are preserved.
  set.seed(123 + m)
  base1 <- generate_base_ICC(n1, t = 1, cluster_size = cluster_size)
  base2 <- generate_base_ICC(n2, t = 2, cluster_size = cluster_size)
  base_test <- generate_base_ICC(test_n, t = 2,
                                 cluster_size = 1, test = TRUE)
  data_test <- generate_data_ICC(base_test, rho = 0)
  splits <- make_split_indices(base1$Z, replicate_id = m, stratum = base1$stratum)

  stopifnot(
    base1$n == n1,
    base2$n == n2,
    length(base1$Z) == base1$n,
    length(base2$Z) == base2$n,
    identical(splits$marginal_cal,
              setdiff(seq_len(base1$n), splits$marginal_train))
  )

  marginal_list <- vector("list", length(rho_vec))
  diagnostic_list <- vector("list", length(rho_vec))
  rho_zero_marginal_diagnostics <- NULL
  rvec_marg_output <- list()

  for (rho_index in seq_along(rho_vec)) {
    rho <- rho_vec[rho_index]
    rho_key <- sprintf("%.2f", rho)
    cat("replicate ", m, ", rho = ", rho_key, "\n", sep = "")
    flush.console()
    if (!is.null(progress_file)) writeLines(
      sprintf("replicate %d: rho = %s calculating", m, rho_key), progress_file
    )

    data1 <- generate_data_ICC(base1, rho = rho)
    data2 <- generate_data_ICC(base2, rho = rho)
    if (m == 1L) {
      diagnostic_list[[rho_index]] <- diagnose_ICC_data(data1, rho)
      if (rho == 0) {
        rho_zero_marginal_diagnostics <-
          diagnose_rho_zero_marginals(data1, t = 1)
      }
    }

    cached_m <- if (is.null(cached_rvec_marg)) NULL else
      cached_rvec_marg[[rho_key]]
    one_result <- run_two_scores(
      data1, data2, data_test, splits,
      replicate_id = m, rho = rho,
      cached_rvec_marg = cached_m,
      R = R, alpha = alpha, n_alpha = n_alpha,
      sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv
    )
    marginal_list[[rho_index]] <- one_result$marginal
    rvec_marg_output[[rho_key]] <- one_result$rvec_marg
  }

  list(
    marginal = bind_rows(marginal_list),
    rvec_marg = rvec_marg_output,
    diagnostics = bind_rows(diagnostic_list),
    rho_zero_marginal_diagnostics = rho_zero_marginal_diagnostics
  )
}


# Runs the ICC sensitivity study sequentially over repetitions.
DA_WCP_simul_ICC <- function(n1 = 5000, n2 = 5000, test_n = 5000,
                             M = 100,
                             rho_vec = seq(0, 0.9, by = 0.1),
                             cluster_size = 50, K = 5, R = 4,
                             alpha = 0.2, n_alpha = 5,
                             rvec = FALSE,
                             rvec_marg = NULL,
                             sig_wo_cv = 1.2, kliep_cv = FALSE,
                             save_outputs = TRUE,
                             output_dir = "../Figures",
                             figures_dir = "../Figures") {
  stopifnot(K == 5, n1 == 5000L, n2 == 5000L, cluster_size == 50L,
            all(rho_vec >= 0), all(rho_vec < 1))
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
    rho_vec = rho_vec, cluster_size = cluster_size,
    K = K, R = R, alpha = alpha, n_alpha = n_alpha,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv
  )

  cat(sprintf("ICC simulation starting: M = %d, %d rho values\n", M, length(rho_vec)))
  flush.console()
    replicate_results <- lapply(
      tasks, function(task) {
        result <- do.call(run_ICC_replicate, c(task, common_args))
        cat(sprintf("%d/%d completed\n", task$replicate_id, M))
        flush.console()
        result
      }
    )
    execution_info <- list(
      backend = "sequential", cores = 1L,
      details = list(), threads = 1L
    )

  results_marginal <- bind_rows(lapply(replicate_results, `[[`, "marginal"))
  rvec_marg_output <- lapply(replicate_results, `[[`, "rvec_marg")
  diagnostics <- bind_rows(lapply(replicate_results, `[[`, "diagnostics"))
  rho_zero_candidates <- lapply(
    replicate_results, `[[`, "rho_zero_marginal_diagnostics"
  )
  rho_zero_candidates <- rho_zero_candidates[!vapply(
    rho_zero_candidates, is.null, logical(1)
  )]
  rho_zero_marginal_diagnostics <- if (length(rho_zero_candidates) > 0L)
    rho_zero_candidates[[1L]] else NULL

  marginal_levels <- c("DA_WCP_AR_m", "DA_WCP_CQR_m")
  results_marginal$method <- factor(results_marginal$method,
                                    levels = marginal_levels)

  stopifnot(
    nrow(results_marginal) == M * length(rho_vec) * 2L,
    all(is.finite(results_marginal$coverage)),
    all(is.finite(results_marginal$length))
  )

  simulation_result <- list(
    results_marginal = results_marginal,
    df_m = results_marginal,
    rvec_marg = rvec_marg_output,
    diagnostics = diagnostics,
    rho_zero_marginal_diagnostics = rho_zero_marginal_diagnostics,
    execution_info = execution_info
  )

  if (isTRUE(save_outputs)) {
    configuration <- list(
      scenario = "S1",
      M = M, rho_vec = rho_vec,
      n1 = n1, n2 = n2, test_n = test_n,
      cluster_size = cluster_size,
      alpha = alpha, K = K, H = 4L, R_kliep = R,
      psus_per_stratum = 25L, observations_per_stratum = 1250L,
      q1_mat = q1_mat, q2_mat = q2_mat,
      n_alpha = n_alpha,
      sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
      methods = c("DA_WCP_AR_m", "DA_WCP_CQR_m"),
      execution_backend = execution_info$backend,
      execution_cores = execution_info$cores,
      threads = 1L,
      test_sample_clustered = FALSE,
      split_type = "within-geographic-stratum individual-level 50:50"
    )
    simulation_result$output_paths <- write_ICC_outputs(
      simulation_result, rho_vec, configuration,
      output_dir, figures_dir
    )
  }

  simulation_result
}


summarize_ICC_results <- function(results_marginal) {
  results_marginal %>%
    group_by(rho, method) %>%
    summarise(
      mean_coverage = mean(coverage),
      sd_coverage = sd(coverage),
      mean_length = mean(length),
      sd_length = sd(length),
      .groups = "drop"
    )
}


save_pdf_plot <- function(plot_object, path, width, height) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_object)
}


make_ICC_figures <- function(results_marginal, rho_vec, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  rho_levels <- format(rho_vec, trim = TRUE, scientific = FALSE)
  plot_m <- results_marginal
  plot_m$rho_plot <- ordered(format(plot_m$rho, trim = TRUE,
                                    scientific = FALSE), levels = rho_levels)

  palette_m <- c(
    DA_WCP_AR_m = "#9ecae1",
    DA_WCP_CQR_m = "#3182bd"
  )

  common_theme <- theme_test(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank()
    )

  coverage_m_plot <- ggplot(plot_m,
                            aes(x = rho_plot, y = coverage, fill = method)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
    geom_hline(yintercept = 0.8, color = "red", linewidth = 0.7) +
    scale_fill_manual(values = palette_m) +
    labs(x = expression(rho), y = "Empirical coverage") +
    common_theme

  length_m_plot <- ggplot(plot_m,
                          aes(x = rho_plot, y = length, fill = method)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
    scale_fill_manual(values = palette_m) +
    labs(x = expression(rho), y = "Average interval length") +
    common_theme

  save_pdf_plot(coverage_m_plot,
                file.path(figures_dir, "simulation_ICC_coverage_m.pdf"),
                width = 10, height = 5)
  save_pdf_plot(length_m_plot,
                file.path(figures_dir, "simulation_ICC_length_m.pdf"),
                width = 10, height = 5)

  invisible(c(
    file.path(figures_dir, "simulation_ICC_coverage_m.pdf"),
    file.path(figures_dir, "simulation_ICC_length_m.pdf")
  ))
}


write_ICC_outputs <- function(simulation_result, rho_vec,
                              configuration, output_dir, figures_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results_marginal <- simulation_result$results_marginal
  df_m <- simulation_result$df_m
  rvec_marg <- simulation_result$rvec_marg
  diagnostics <- simulation_result$diagnostics
  execution_info <- simulation_result$execution_info
  rho_zero_marginal_diagnostics <-
    simulation_result$rho_zero_marginal_diagnostics

  summary_marginal <- summarize_ICC_results(results_marginal)

  marginal_csv <- file.path(
    output_dir, "DA_WCP_supp_ICC_summary_marginal.csv")
  rdata_path <- file.path(output_dir, "DA_WCP_supp_ICC.RData")

  write.csv(summary_marginal, marginal_csv, row.names = FALSE)

  figure_paths <- make_ICC_figures(results_marginal, rho_vec, figures_dir)

  save(
    results_marginal, df_m, summary_marginal, rvec_marg,
    diagnostics, rho_zero_marginal_diagnostics,
    execution_info, configuration,
    file = rdata_path
  )

  output_paths <- c(rdata_path, marginal_csv, figure_paths)
  stopifnot(all(file.exists(output_paths)))
  output_paths
}


# Nothing runs automatically when this file is sourced. Start the production
# simulation explicitly from R/RStudio with:
# ICC_result <- DA_WCP_simul_ICC(M = 100)

if (dir.exists("Code") && file.exists("Code/KLIEP.R")) setwd("Code")
dir.create("../Figures", recursive = TRUE, showWarnings = FALSE)

# Serial supplementary simulations for survey strata and PSU composition.
# This file combines the H and H-PSU studies and runs every repetition sequentially.
# The data-generating mechanisms, estimators, seeds, diagnostics, and outputs match the source studies.

library(mvtnorm)
library(dplyr)
library(quantreg)
library(ggplot2)

# Load the shared implementations.
source("KLIEP.R")
source("DA_WCP.R")


###############################################################################
# Production defaults and scenario S1 constants
###############################################################################

H <- 4
M <- 100
n1 <- 5000
n2 <- 5000
test_n <- 5000

K <- 5
alpha <- 0.2
R_kliep <- 4
n_alpha <- 5
sig_wo_cv <- 1.2
kliep_cv <- FALSE

target_g <- c(0.05, 0.10, 0.15, 0.30, 0.40)

q1_mat <- rbind(
  c(0.30, 0.25, 0.20, 0.15, 0.10),
  c(0.10, 0.15, 0.20, 0.25, 0.30),
  c(0.25, 0.15, 0.10, 0.20, 0.30),
  c(0.15, 0.25, 0.30, 0.20, 0.10)
)

q2_mat <- rbind(
  c(0.12, 0.17, 0.20, 0.24, 0.27),
  c(0.08, 0.13, 0.20, 0.26, 0.33),
  c(0.11, 0.14, 0.22, 0.23, 0.30),
  c(0.09, 0.16, 0.18, 0.27, 0.30)
)

historical_S1_g <- rep(0.20, 5)
current_S1_g <- c(0.10, 0.15, 0.20, 0.25, 0.30)

stopifnot(
  H == 4,
  isTRUE(all.equal(colMeans(q1_mat), historical_S1_g, tolerance = 1e-15)),
  isTRUE(all.equal(colMeans(q2_mat), current_S1_g, tolerance = 1e-15)),
  isTRUE(all.equal(rowSums(q1_mat), rep(1, H), tolerance = 1e-15)),
  isTRUE(all.equal(rowSums(q2_mat), rep(1, H), tolerance = 1e-15)),
  identical(target_g, c(0.05, 0.10, 0.15, 0.30, 0.40)),
  n1 == 5000,
  n2 == 5000,
  n1 %% H == 0,
  n2 %% H == 0
)

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

beta <- c(5, -3, 2, 1)
gamma <- c(-10, -5, 0, 5, 10)
theta_hetero <- c(0.5, -0.4, 0.3, 0.2)
eta_hetero <- c(1.2, 1.6, 2.0, 2.4, 2.8)

var_x <- c("X1", "X2", "X3", "X4")
var_y <- "Y"
var_z <- c()
var_g <- "Z"
var_h <- "stratum"


###############################################################################
# H = 4 data generation
###############################################################################

generate_XY_S1 <- function(Z, t) {
  n <- length(Z)
  X <- matrix(0, nrow = n, ncol = 4)
  Y <- numeric(n)

  for (k in 1:5) {
    idx <- which(Z == k)
    if (length(idx) == 0L) next
    mu_k <- mu1(k) + if (t == 2) Delta_mu_shift else 0
    Sigma_k <- cov1(k) + if (t == 2) Delta_Sigma_shift else 0
    X[idx, ] <- MASS::mvrnorm(n = length(idx), mu = mu_k, Sigma = Sigma_k)
    for (i in idx) {
      Y[i] <- X[i, ] %*% beta + gamma[k] +
        rnorm(1, sd = sqrt(exp(X[i, ] %*% theta_hetero + eta_hetero[k])))
    }
  }

  list(X = X, Y = Y)
}


generate_data_H <- function(n, t, test = FALSE, H = 4) {
  if (length(n) != 1L || !is.finite(n) || n < 1L || n != as.integer(n)) {
    stop("n must be one positive integer.")
  }
  n <- as.integer(n)
  stopifnot(t %in% c(1, 2))

  if (test) {
    Z <- sample(1:5, size = n, replace = TRUE, prob = target_g)
    generated <- generate_XY_S1(Z, t = 2)
    data <- as.data.frame(generated$X)
    names(data) <- var_x
    data$Z <- Z
    data$Y <- generated$Y
    data$time <- 2
    attr(data, "stratum_sizes") <- NULL
    return(data)
  }

  if (!H %in% c(1L, 4L)) {
    stop("Survey generation in this supplementary design requires H = 1 or 4.")
  }
  if (n %% H != 0L) stop("n must be divisible by H for equal fixed stratum sizes.")
  fixed_n_h <- as.integer(n / H)
  q_mat <- if (H == 1L) {
    matrix(if (t == 1) historical_S1_g else current_S1_g, nrow = 1L)
  } else if (t == 1) {
    q1_mat
  } else {
    q2_mat
  }
  stratum_data <- vector("list", H)
  n_h <- rep(fixed_n_h, H)

  for (h in seq_len(H)) {
    Z <- sample(1:5, size = n_h[h], replace = TRUE, prob = q_mat[h, ])
    generated <- generate_XY_S1(Z, t = t)

    data_h <- as.data.frame(generated$X)
    names(data_h) <- var_x
    data_h$Z <- Z
    if (t == 1) data_h$Y <- generated$Y
    data_h$time <- t
    data_h$stratum <- h
    stratum_data[[h]] <- data_h
  }

  data <- bind_rows(stratum_data)
  stopifnot(nrow(data) == n)
  attr(data, "stratum_sizes") <- n_h
  data
}


###############################################################################
# DA_WCP_H: compatibility wrapper for marginal DA-WCP
###############################################################################

# Compatibility wrapper: all interval construction is delegated to DA_WCP.R.
# Keep the existing supplementary return fields, including weight diagnostics.
DA_WCP_H <- function(dtra, dcal, dttra, d2, dtest,
                     var_x, var_y, var_z, var_g, var_h,
                     val_type, score_type,
                     var_wt_x = NA, var_wt_y = NA,
                     alpha = 0.2, p2_vec = NA, n_alpha = 5,
                     rvec = FALSE, rvec_cal = NA, R = 4,
                     sig_wo_cv = 1.2, kliep_cv = FALSE,
                     simul = FALSE, eval = TRUE) {
  if (val_type != "marg" || !score_type %in% c("ar", "cqr")) {
    stop("DA_WCP_H supports marginal AR and CQR only.")
  }
  DA_WCP(
    dtra, dcal, dttra, d2, dtest, var_x, var_y, var_z, var_g, var_h,
    val_type, score_type, var_wt_x = var_wt_x, var_wt_y = var_wt_y,
    alpha = alpha, p2_vec = p2_vec, n_alpha = n_alpha,
    rvec = rvec, rvec_cal = rvec_cal, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
    simul = simul, eval = FALSE, return_diagnostics = TRUE
  )
}


###############################################################################
# Reproducible split construction and method evaluation
###############################################################################

make_splits_H <- function(data1, replicate_id, H = 4) {
  set.seed(100000L + as.integer(replicate_id))
  marginal_train <- integer(0)
  marginal_cal <- integer(0)
  for (h in seq_len(H)) {
    idx <- which(data1$stratum == h)
    local_train <- sort(sample.int(length(idx), length(idx) %/% 2L))
    local_cal <- setdiff(seq_along(idx), local_train)
    marginal_train <- c(marginal_train, idx[local_train])
    marginal_cal <- c(marginal_cal, idx[local_cal])
  }

  stopifnot(
    length(intersect(marginal_train, marginal_cal)) == 0L,
    setequal(c(marginal_train, marginal_cal), seq_len(nrow(data1))),
    setequal(unique(data1$stratum[marginal_train]), seq_len(H)),
    setequal(unique(data1$stratum[marginal_cal]), seq_len(H))
  )

  list(
    marginal_train = marginal_train,
    marginal_cal = marginal_cal
  )
}


interval_metrics <- function(result_mat, test) {
  covered <- test[[var_y]] >= result_mat[, 1] &
    test[[var_y]] <= result_mat[, 2]
  interval_length <- result_mat[, 2] - result_mat[, 1]
  c(coverage = mean(covered), length = mean(interval_length))
}


run_H1_equivalence_check <- function(train, cal, data2_x, test,
                                     rvec_cal, alpha, n_alpha, R,
                                     sig_wo_cv, kliep_cv,
                                     tolerance = 1e-10) {
  train_h1 <- train
  cal_h1 <- cal
  data2_h1 <- data2_x
  train_h1$stratum <- 1L
  cal_h1$stratum <- 1L
  data2_h1$stratum <- 1L

  custom <- DA_WCP_H(
    train_h1, cal_h1, train_h1, data2_h1, test,
    var_x, var_y, var_z, var_g, var_h,
    val_type = "marg", score_type = "ar",
    p2_vec = target_g, alpha = alpha, n_alpha = n_alpha,
    rvec = TRUE, rvec_cal = rvec_cal,
    R = R, sig_wo_cv = sig_wo_cv,
    kliep_cv = kliep_cv, simul = TRUE
  )
  original <- DA_WCP(
    train_h1, cal_h1, train_h1, data2_h1, test,
    var_x, var_y, var_z, var_g, var_h,
    val_type = "marg", score_type = "ar",
    p2_vec = target_g, alpha = alpha, n_alpha = n_alpha,
    rvec = TRUE, rvec_cal = rvec_cal,
    R = R, sig_wo_cv = sig_wo_cv,
    kliep_cv = kliep_cv, simul = TRUE
  )

  max_abs_difference <- max(abs(custom$result_mat - original$result_mat))
  stopifnot(is.finite(max_abs_difference), max_abs_difference <= tolerance)
  max_abs_difference
}


run_four_marginal_methods_H <- function(data1, data2, data_test, splits,
                                        replicate_id, H = 4,
                                        cached_rvec_marg = NULL,
                                        alpha = 0.2, n_alpha = 5, R = 4,
                                        sig_wo_cv = 1.2,
                                        kliep_cv = FALSE,
                                        run_h1_check = FALSE) {
  data1_y <- data1[, c(var_x, var_y, var_g, var_h), drop = FALSE]
  data2_x <- data2[, c(var_x, var_g, var_h), drop = FALSE]
  test <- data_test[, c(var_x, var_y, var_g), drop = FALSE]

  train_m <- data1_y[splits$marginal_train, , drop = FALSE]
  cal_m <- data1_y[splits$marginal_cal, , drop = FALSE]

  use_cached_m <- !is.null(cached_rvec_marg)
  da_ar_m <- DA_WCP(
    train_m, cal_m, train_m, data2_x, test,
    var_x, var_y, var_z, var_g, var_h,
    val_type = "marg", score_type = "ar",
    p2_vec = target_g, alpha = alpha, n_alpha = n_alpha,
    rvec = use_cached_m,
    rvec_cal = if (use_cached_m) cached_rvec_marg else NA,
    R = R, sig_wo_cv = sig_wo_cv,
    kliep_cv = kliep_cv, simul = TRUE, return_diagnostics = TRUE
  )
  rvec_marg_current <- da_ar_m$rvec_cal

  da_cqr_m <- DA_WCP(
    train_m, cal_m, train_m, data2_x, test,
    var_x, var_y, var_z, var_g, var_h,
    val_type = "marg", score_type = "cqr",
    p2_vec = target_g, alpha = alpha, n_alpha = n_alpha,
    rvec = TRUE, rvec_cal = rvec_marg_current,
    R = R, sig_wo_cv = sig_wo_cv,
    kliep_cv = kliep_cv, simul = TRUE, return_diagnostics = TRUE
  )

  # Use the unchanged GWCP implementation from DA_WCP.R. As in the main
  # simulation, marginal GWCP uses demographic composition weights only and
  # therefore neither estimates nor reuses a KLIEP density ratio.
  gwcp_ar_m <- DA_WCP(
    train_m, cal_m, train_m, data2_x, test,
    var_x, var_y, var_z, var_g, var_h,
    gwcp = TRUE, val_type = "marg", score_type = "ar",
    p2_vec = target_g, alpha = alpha, n_alpha = n_alpha,
    R = R, sig_wo_cv = sig_wo_cv,
    kliep_cv = kliep_cv, simul = TRUE
  )
  gwcp_cqr_m <- DA_WCP(
    train_m, cal_m, train_m, data2_x, test,
    var_x, var_y, var_z, var_g, var_h,
    gwcp = TRUE, val_type = "marg", score_type = "cqr",
    p2_vec = target_g, alpha = alpha, n_alpha = n_alpha,
    R = R, sig_wo_cv = sig_wo_cv,
    kliep_cv = kliep_cv, simul = TRUE
  )

  marginal_methods <- list(
    DA_WCP_AR_m = da_ar_m$result_mat,
    GWCP_AR_m = gwcp_ar_m$result_mat,
    DA_WCP_CQR_m = da_cqr_m$result_mat,
    GWCP_CQR_m = gwcp_cqr_m$result_mat
  )
  results_marginal <- do.call(rbind, lapply(names(marginal_methods), function(method) {
    metrics <- interval_metrics(marginal_methods[[method]], test)
    data.frame(
      replicate = replicate_id,
      H = H,
      method = method,
      coverage = unname(metrics["coverage"]),
      length = unname(metrics["length"]),
      stringsAsFactors = FALSE
    )
  }))

  h1_max_abs_difference <- NA_real_
  if (run_h1_check) {
    h1_max_abs_difference <- run_H1_equivalence_check(
      train_m, cal_m, data2_x, test, rvec_marg_current,
      alpha, n_alpha, R, sig_wo_cv, kliep_cv
    )
  }

  list(
    marginal = results_marginal,
    rvec_marg = rvec_marg_current,
    weight_mass_by_stratum = da_ar_m$weight_mass_by_stratum,
    q1h_hat = da_ar_m$q1h_mat,
    h1_max_abs_difference = h1_max_abs_difference
  )
}


###############################################################################
# Simulation driver
###############################################################################

run_H_replicate <- function(replicate_id, n1, n2, test_n, H_values,
                            cached_rvec_marg = NULL,
                            alpha = 0.2, n_alpha = 5, R = 4,
                            sig_wo_cv = 1.2, kliep_cv = FALSE,
                            run_h1_check = FALSE) {
  m <- replicate_id
  cat("------------------------------------------------------\n")
  cat("     H = 1 versus H = 4 - replicate ", m, "\n", sep = "")
  cat("------------------------------------------------------\n")

  H_results <- vector("list", length(H_values))
  H_sample_sizes <- vector("list", length(H_values))
  H_rvec <- setNames(vector("list", length(H_values)), paste0("H", H_values))
  H_diagnostics <- setNames(
    vector("list", length(H_values)), paste0("H", H_values)
  )

  for (H_index in seq_along(H_values)) {
    H_current <- H_values[H_index]
    cat("H = ", H_current, "\n", sep = "")

    # Reset to the same original S1 replicate seed for each H arm. This keeps
    # the random-seed convention identical while changing only the stratum
    # design and its associated demographic compositions.
    set.seed(123 + m)
    data1 <- generate_data_H(n1, t = 1, H = H_current)
    n1h <- attr(data1, "stratum_sizes")
    data2 <- generate_data_H(n2, t = 2, H = H_current)
    n2h <- attr(data2, "stratum_sizes")
    data_test <- generate_data_H(test_n, t = 2, test = TRUE)

    stopifnot(
      setequal(unique(data1$stratum), seq_len(H_current)),
      setequal(unique(data2$stratum), seq_len(H_current)),
      nrow(data1) == n1,
      nrow(data2) == n2,
      all(n1h == n1 / H_current),
      all(n2h == n2 / H_current),
      length(n1h) == H_current,
      length(n2h) == H_current,
      !"stratum" %in% names(data_test)
    )

    cached_H_rvec <- NULL
    if (!is.null(cached_rvec_marg)) {
      H_name <- paste0("H", H_current)
      if (is.list(cached_rvec_marg) && H_name %in% names(cached_rvec_marg)) {
        cached_H_rvec <- cached_rvec_marg[[H_name]]
      } else if (H_current == 4L && is.numeric(cached_rvec_marg)) {
        # Backward compatibility with rvec objects saved by the earlier
        # H = 4-only version of this supplementary simulation.
        cached_H_rvec <- cached_rvec_marg
      }
    }

    splits <- make_splits_H(data1, replicate_id = m, H = H_current)
    one_result <- run_four_marginal_methods_H(
      data1, data2, data_test, splits,
      replicate_id = m, H = H_current,
      cached_rvec_marg = cached_H_rvec,
      alpha = alpha, n_alpha = n_alpha, R = R,
      sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
      run_h1_check = run_h1_check && H_current == 1L
    )

    H_results[[H_index]] <- one_result$marginal
    H_sample_sizes[[H_index]] <- data.frame(
      replicate = m, H = H_current, stratum = seq_len(H_current),
      n1h = n1h, n2h = n2h
    )
    H_rvec[[paste0("H", H_current)]] <- one_result$rvec_marg
    H_diagnostics[[paste0("H", H_current)]] <- list(
      q1h_hat = one_result$q1h_hat,
      weight_mass_by_stratum = one_result$weight_mass_by_stratum,
      h1_max_abs_difference = one_result$h1_max_abs_difference
    )
  }

  list(
    marginal = bind_rows(H_results),
    rvec_marg = H_rvec,
    sample_sizes = bind_rows(H_sample_sizes),
    method_diagnostics = list(
      replicate = m, by_H = H_diagnostics
    )
  )
}


# Runs the H = 1 and H = 4 designs sequentially over repetitions.
DA_WCP_simul_H <- function(n1 = 5000, n2 = 5000, test_n = 5000,
                           H = 4, M = 100, K = 5,
                           alpha = 0.2, n_alpha = 5, R = 4,
                           sig_wo_cv = 1.2, kliep_cv = FALSE,
                           rvec = FALSE,
                           rvec_marg = NULL,
                           run_h1_check = TRUE,
                           save_outputs = TRUE,
                           output_dir = "../Figures",
                           figures_dir = "../Figures") {
  stopifnot(H == 4, K == 5)
  H_values <- c(1L, 4L)
  if (rvec && is.null(rvec_marg)) {
    stop("Cached rvec_marg is required when rvec = TRUE.")
  }

  tasks <- lapply(seq_len(M), function(m) {
    list(
      replicate_id = m,
      cached_rvec_marg = if (rvec) rvec_marg[[m]] else NULL,
      run_h1_check = run_h1_check && m == 1L
    )
  })
  common_args <- list(
    n1 = n1, n2 = n2, test_n = test_n, H_values = H_values,
    alpha = alpha, n_alpha = n_alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv
  )

    replicate_results <- lapply(
      tasks, function(task) do.call(run_H_replicate, c(task, common_args))
    )
    execution_info <- list(
      backend = "sequential", cores = 1L,
      details = list(), threads = 1L
    )

  results_marginal <- bind_rows(lapply(replicate_results, `[[`, "marginal"))
  rvec_marg_output <- lapply(replicate_results, `[[`, "rvec_marg")
  sample_sizes <- bind_rows(lapply(replicate_results, `[[`, "sample_sizes"))
  method_diagnostics <- lapply(
    replicate_results, `[[`, "method_diagnostics"
  )

  marginal_levels <- c(
    "DA_WCP_AR_m", "GWCP_AR_m", "DA_WCP_CQR_m", "GWCP_CQR_m"
  )
  results_marginal$method <- factor(results_marginal$method, levels = marginal_levels)

  stopifnot(
    nrow(results_marginal) == M * length(H_values) * 4L,
    setequal(unique(results_marginal$H), H_values),
    all(is.finite(results_marginal$coverage)),
    all(is.finite(results_marginal$length))
  )

  simulation_result <- list(
    results_marginal = results_marginal,
    df_m = results_marginal,
    rvec_marg = rvec_marg_output,
    sample_sizes = sample_sizes,
    method_diagnostics = method_diagnostics,
    execution_info = execution_info
  )

  if (isTRUE(save_outputs)) {
    configuration <- list(
      scenario = "S1 comparison of H = 1 and H = 4",
      H = H_values, M = M, n1 = n1, n2 = n2, test_n = test_n,
      fixed_stratum_size_n1 = setNames(n1 / H_values, paste0("H", H_values)),
      fixed_stratum_size_n2 = setNames(n2 / H_values, paste0("H", H_values)),
      H1_historical_g = historical_S1_g,
      H1_current_g = current_S1_g,
      H4_historical_g_by_stratum = q1_mat,
      H4_current_g_by_stratum = q2_mat,
      replicate_seed = "123 + replicate, reset separately for H = 1 and H = 4",
      alpha = alpha, n_alpha = n_alpha, K = K, R_kliep = R,
      sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
      methods = marginal_levels,
      marginal_split = "within stratum 50:50",
      target_has_survey_stratum = FALSE,
      execution_backend = execution_info$backend,
      execution_cores = execution_info$cores,
      threads = 1L
    )
    simulation_result$output_paths <- write_H_outputs(
      simulation_result, configuration, output_dir, figures_dir
    )
  }

  simulation_result
}


###############################################################################
# Summaries, figures, and saving
###############################################################################

summarize_H_results <- function(results_marginal) {
  results_marginal %>%
    group_by(H, method) %>%
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


make_H_figures <- function(results_marginal, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  plot_data <- results_marginal
  plot_data$H_label <- factor(
    paste0("H = ", plot_data$H), levels = c("H = 1", "H = 4")
  )
  palette_H <- c("H = 1" = "#9ecae1", "H = 4" = "#e6550d")

  coverage_m <- ggplot(
    plot_data, aes(x = coverage, y = method, fill = H_label)
  ) +
    geom_boxplot(
      position = position_dodge2(width = 0.75, preserve = "single"),
      width = 0.65, outlier.size = 0.8
    ) +
    geom_vline(
      xintercept = 0.8, color = "red", linetype = "dashed", linewidth = 0.7
    ) +
    scale_fill_manual(values = palette_H) +
    labs(x = "Empirical coverage", y = "", fill = "Survey strata") +
    theme_test(base_size = 14) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())

  length_m <- ggplot(
    plot_data, aes(x = length, y = method, fill = H_label)
  ) +
    geom_boxplot(
      position = position_dodge2(width = 0.75, preserve = "single"),
      width = 0.65, outlier.size = 0.8
    ) +
    scale_fill_manual(values = palette_H) +
    labs(x = "Average interval length", y = "", fill = "Survey strata") +
    theme_test(base_size = 14) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())

  paths <- c(
    file.path(figures_dir, "simulation_H1_H4_coverage_m.pdf"),
    file.path(figures_dir, "simulation_H1_H4_length_m.pdf")
  )
  save_pdf_plot(coverage_m, paths[1], 8, 4.5)
  save_pdf_plot(length_m, paths[2], 8, 4.5)
  invisible(paths)
}


write_H_outputs <- function(simulation_result, configuration,
                            output_dir, figures_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results_marginal <- simulation_result$results_marginal
  df_m <- simulation_result$df_m
  rvec_marg <- simulation_result$rvec_marg
  sample_sizes <- simulation_result$sample_sizes
  method_diagnostics <- simulation_result$method_diagnostics
  execution_info <- simulation_result$execution_info
  summary_marginal <- summarize_H_results(results_marginal)

  marginal_csv <- file.path(output_dir, "DA_WCP_supp_H_summary_marginal.csv")
  rdata_path <- file.path(output_dir, "DA_WCP_supp_H.RData")
  write.csv(summary_marginal, marginal_csv, row.names = FALSE)
  figure_paths <- make_H_figures(results_marginal, figures_dir)

  save(
    results_marginal, df_m, summary_marginal,
    rvec_marg, sample_sizes, method_diagnostics,
    execution_info, q1_mat, q2_mat, target_g, configuration,
    file = rdata_path
  )

  output_paths <- c(rdata_path, marginal_csv, figure_paths)
  stopifnot(all(file.exists(output_paths)))
  output_paths
}


# Nothing runs automatically when this file is sourced. Start the production
# simulation explicitly from R/RStudio with:
# H_result <- DA_WCP_simul_H()


###############################################################################
# PSU-composition sensitivity study
###############################################################################

# Edit this object to change the sensitivity grid.
kappa_grid <- c(1, 2, 5, 10, 25, 50, 100)

generate_data_H_PSU <- function(t, kappa) {
  stopifnot(t %in% c(1L, 2L), length(kappa) == 1L,
            !is.na(kappa), kappa > 0)
  q_mat <- if (t == 1L) q1_mat else q2_mat
  strata <- vector("list", 4L)
  diagnostics <- vector("list", 100L)
  for (h in seq_len(4L)) {
    labels <- integer(1250L)
    for (j in seq_len(25L)) {
      q <- if (is.infinite(kappa)) q_mat[h, ] else {
        draw <- rgamma(5L, shape = kappa * q_mat[h, ], rate = 1)
        if (!is.finite(sum(draw)) || sum(draw) <= 0)
          stop("Invalid Dirichlet draw; check kappa.")
        draw / sum(draw)
      }
      idx <- (j - 1L) * 50L + seq_len(50L)
      labels[idx] <- sample.int(5L, 50L, replace = TRUE, prob = q)
      observed_q <- tabulate(labels[idx], nbins = 5L) / 50
      one <- data.frame(wave = t, stratum = h, psu = j, n = 50L,
                        l1_composition = sum(abs(q - q_mat[h, ])),
                        l1_observed = sum(abs(observed_q - q_mat[h, ])))
      for (k in seq_len(5L)) {
        one[[paste0("q", k)]] <- q[k]
        one[[paste0("baseline_q", k)]] <- q_mat[h, k]
        one[[paste0("observed_q", k)]] <- observed_q[k]
      }
      diagnostics[[(h - 1L) * 25L + j]] <- one
    }
    # Original independent S1 X,Y draws conditional on each individual's Z.
    # Generate across the pooled stratum; PSU identity never enters this helper.
    xy <- generate_XY_S1(labels, t)
    d <- as.data.frame(xy$X)
    names(d) <- var_x
    d$Z <- labels
    if (t == 1L) d$Y <- xy$Y
    d$time <- t
    d$stratum <- h
    d$psu <- rep(seq_len(25L), each = 50L)
    strata[[h]] <- d
  }
  data <- dplyr::bind_rows(strata)
  diagnostic <- dplyr::bind_rows(diagnostics)
  stopifnot(nrow(data) == 5000L,
            all(table(data$stratum) == 1250L),
            all(table(data$stratum, data$psu) == 50L),
            nrow(diagnostic) == 100L,
            max(abs(rowSums(diagnostic[paste0("q", 1:5)]) - 1)) < 1e-12)
  if (is.infinite(kappa)) stopifnot(all(diagnostic$l1_composition == 0))
  list(data = data, diagnostic = diagnostic)
}

analyze_pooled_H_PSU <- function(data1, data2, data_test, replicate_id) {
  # Explicit whitelist removes PSU labels before BOTH splitting and estimation.
  data1 <- data1[, c(var_x, var_y, var_g, var_h), drop = FALSE]
  data2 <- data2[, c(var_x, var_g, var_h), drop = FALSE]
  test <- data_test[, c(var_x, var_y, var_g), drop = FALSE]
  stopifnot(!"psu" %in% names(data1), !"psu" %in% names(data2))
  splits <- make_splits_H(data1, replicate_id = replicate_id, H = 4)
  train <- data1[splits$marginal_train, , drop = FALSE]
  cal <- data1[splits$marginal_cal, , drop = FALSE]
  stopifnot(all(table(train$stratum) == 625L),
            all(table(cal$stratum) == 625L))
  # A zero subgroup count is not silently repaired or resampled.
  if (any(table(factor(train$stratum, levels = 1:4),
                factor(train$Z, levels = 1:5)) == 0L))
    stop("A training stratum/subgroup cell is empty; original weights undefined.")

  common <- list(dtra = train, dcal = cal, dttra = train, d2 = data2,
                 dtest = test, var_x = var_x, var_y = var_y, var_z = var_z,
                 var_g = var_g, var_h = var_h, val_type = "marg",
                 p2_vec = target_g, alpha = 0.2, n_alpha = 5,
                 R = 4, sig_wo_cv = 1.2, kliep_cv = FALSE, simul = TRUE)
  ar <- do.call(DA_WCP, c(common, list(return_diagnostics = TRUE, score_type = "ar", rvec = FALSE)))
  cqr <- do.call(DA_WCP, c(common, list(return_diagnostics = TRUE, score_type = "cqr", rvec = TRUE,
                                        rvec_cal = ar$rvec_cal)))
  metrics <- rbind(AR = interval_metrics(ar$result_mat, test),
                   CQR = interval_metrics(cqr$result_mat, test))
  stopifnot(all(is.finite(metrics)), all(metrics[, "coverage"] >= 0),
            all(metrics[, "coverage"] <= 1),
            all(abs(ar$weight_mass_by_stratum - 0.25) < 1e-10),
            all(abs(cqr$weight_mass_by_stratum - 0.25) < 1e-10))
  list(marginal = data.frame(score = rownames(metrics), coverage = metrics[, "coverage"],
                             length = metrics[, "length"], row.names = NULL),
       rvec_marg = ar$rvec_cal,
       method_diagnostics = list(q1h_hat = ar$q1h_mat,
                                  weight_mass_by_stratum = ar$weight_mass_by_stratum))
}

summarize_H_PSU_results <- function(results_marginal) {
  results_marginal %>%
    dplyr::group_by(kappa, method, score) %>%
    dplyr::summarise(
      mean_coverage = mean(coverage), sd_coverage = sd(coverage),
      mean_length = mean(length), sd_length = sd(length), .groups = "drop")
}

make_H_PSU_figures <- function(results_marginal, kappa_values, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  plot_m <- results_marginal
  plot_m$kappa_plot <- ordered(as.character(plot_m$kappa),
                               levels = as.character(kappa_values))
  # Same palette, boxplots, theme and dimensions as make_ICC_figures().
  palette_m <- c(DA_WCP_AR_m = "#9ecae1", DA_WCP_CQR_m = "#3182bd")
  common_theme <- ggplot2::theme_test(base_size = 14) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   legend.position = "bottom", legend.title = ggplot2::element_blank())
  coverage_m_plot <- ggplot2::ggplot(
    plot_m, ggplot2::aes(x = kappa_plot, y = coverage, fill = method)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8),
                          outlier.size = 0.8) +
    ggplot2::geom_hline(yintercept = 0.8, color = "red", linewidth = 0.7) +
    ggplot2::scale_y_continuous(breaks = seq(0.60, 0.90, by = 0.05),
                                labels = function(x) sprintf("%.2f", x),
                                expand = ggplot2::expansion(mult = 0.03)) +
    ggplot2::coord_cartesian(ylim = c(0.60, 0.90)) +
    ggplot2::scale_fill_manual(values = palette_m) +
    ggplot2::labs(x = expression(kappa), y = "Empirical coverage") + common_theme
  length_m_plot <- ggplot2::ggplot(
    plot_m, ggplot2::aes(x = kappa_plot, y = length, fill = method)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8),
                          outlier.size = 0.8) +
    ggplot2::scale_fill_manual(values = palette_m) +
    ggplot2::labs(x = expression(kappa), y = "Average interval length") + common_theme
  paths <- file.path(figures_dir, c("simulation_H_PSU_coverage_m.pdf",
                                    "simulation_H_PSU_length_m.pdf"))
  save_pdf_plot(coverage_m_plot, paths[1], width = 10, height = 5)
  save_pdf_plot(length_m_plot, paths[2], width = 10, height = 5)
  invisible(paths)
}

write_H_PSU_outputs <- function(simulation_result, configuration, output_dir,
                                figures_dir, make_figure = TRUE) {
  # Match write_H_outputs(): named R objects in .RData plus a summary CSV.
  # Repetition-level results and diagnostics are inside .RData, not separate CSVs.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  results_marginal <- simulation_result$results_marginal
  df_m <- simulation_result$df_m
  rvec_marg <- simulation_result$rvec_marg
  sample_sizes <- simulation_result$sample_sizes
  method_diagnostics <- simulation_result$method_diagnostics
  diagnostics <- simulation_result$diagnostics
  execution_info <- simulation_result$execution_info
  summary_marginal <- summarize_H_PSU_results(results_marginal)
  diagnostics_by_repetition <- diagnostics %>%
    dplyr::group_by(kappa_index, kappa, replicate) %>%
    dplyr::summarise(mean_l1_composition = mean(l1_composition),
                     mean_l1_observed = mean(l1_observed), .groups = "drop")
  kappa_grid <- configuration$kappa_grid
  q1_mat <- configuration$q1_mat
  q2_mat <- configuration$q2_mat
  target_g <- configuration$target_g
  marginal_csv <- file.path(output_dir, "DA_WCP_supp_H_PSU_summary_marginal.csv")
  rdata_path <- file.path(output_dir, "DA_WCP_supp_H_PSU.RData")
  write.csv(summary_marginal, marginal_csv, row.names = FALSE)
  save(results_marginal, df_m, summary_marginal, rvec_marg, sample_sizes,
       method_diagnostics, diagnostics, diagnostics_by_repetition,
       execution_info, kappa_grid, q1_mat, q2_mat, target_g, configuration,
       file = rdata_path)
  figure_paths <- if (make_figure)
    make_H_PSU_figures(results_marginal, kappa_grid, figures_dir) else character()
  output_paths <- c(rdata_path, marginal_csv, figure_paths)
  stopifnot(all(file.exists(output_paths)))
  invisible(output_paths)
}

run_H_PSU_task <- function(task, seed) {
  ki <- task$kappa_index
  m <- task$replicate
  kappa <- task$kappa
  data_seed <- as.integer(seed + ki * 1000000 + m * 3)
  # Seeds depend only on the task, only on the task index.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  start <- proc.time()[["elapsed"]]
  task_warnings <- character()
  tryCatch(withCallingHandlers({
    set.seed(data_seed)
    wave1 <- generate_data_H_PSU(1L, kappa)
    set.seed(data_seed + 1L)
    wave2 <- generate_data_H_PSU(2L, kappa)
    set.seed(data_seed + 2L)
    test <- generate_data_H(5000L, t = 2L, test = TRUE)
    analysis <- analyze_pooled_H_PSU(wave1$data, wave2$data, test, m)
    metric <- analysis$marginal
    diag <- dplyr::bind_rows(wave1$diagnostic, wave2$diagnostic)
    metric$kappa_index <- ki
    metric$kappa <- kappa
    metric$replicate <- m
    metric$data_seed <- data_seed
    metric$method <- factor(paste0("DA_WCP_", metric$score, "_m"),
                            levels = c("DA_WCP_AR_m", "DA_WCP_CQR_m"))
    metric$H <- 4L
    diag$kappa_index <- ki
    diag$kappa <- kappa
    diag$replicate <- m
    list(metric = metric, diagnostic = diag, rvec_marg = analysis$rvec_marg,
         method_diagnostics = analysis$method_diagnostics,
         sample_sizes = data.frame(kappa = kappa, replicate = m, H = 4L,
                                    stratum = 1:4, n1h = 1250L, n2h = 1250L),
         warnings = task_warnings,
         elapsed_seconds = proc.time()[["elapsed"]] - start)
  }, warning = function(w) {
    task_warnings <<- c(task_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e), warnings = task_warnings))
}

# Runs the PSU-composition sensitivity grid sequentially.
DA_WCP_simul_H_PSU <- function(M = 100L, kappa_values = kappa_grid,
                               seed = 20260912L,
                               output_dir = "../Figures",
                               figures_dir = "../Figures",
                               save_outputs = TRUE) {
  stopifnot(length(M) == 1L, is.finite(M), M >= 1, M == floor(M), M <= 100000L,
            is.numeric(kappa_values), length(kappa_values) > 0L,
            all(!is.na(kappa_values)), all(kappa_values > 0),
            !anyDuplicated(kappa_values), length(seed) == 1L,
            is.finite(seed), seed >= 0, seed == floor(seed),
            seed + length(kappa_values) * 1000000 + M * 3 < .Machine$integer.max)
  tasks <- expand.grid(replicate = seq_len(M), kappa_index = seq_along(kappa_values))
  tasks$kappa <- kappa_values[tasks$kappa_index]
  if (is.null(output_dir)) output_dir <- "../Figures"
  if (isTRUE(save_outputs)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
  }
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  configuration <- list(M = M, kappa_grid = kappa_values, seed = seed,
    H = 4L, psus_per_stratum = 25L, individuals_per_psu = 50L,
    n_per_stratum = 1250L, n1 = 5000L, n2 = 5000L, test_n = 5000L,
    alpha = 0.2, n_alpha = 5L, R = 4L, sigma = 1.2, kliep_cv = FALSE,
    q1_mat = q1_mat, q2_mat = q2_mat, target_g = target_g,
    execution = "sequential", execution_unit = "kappa x repetition", conditional_PSU_random_effects = FALSE,
    seed_rule = "seed + kappa_index*1000000 + replicate*3 + (wave-1); test uses +2; split 100000+replicate",
    diagnostic = "L1 of latent PSU probability vs stratum probability; observed label proportions stored separately")
  configuration$session_info <- sessionInfo()
  results <- vector("list", nrow(tasks))
  diagnostic <- vector("list", nrow(tasks))
  rvec_output <- vector("list", nrow(tasks))
  sample_size_output <- vector("list", nrow(tasks))
  method_diagnostic_output <- vector("list", nrow(tasks))
  execution_info <- list(backend = "sequential", mode = "sequential",
                        threads = 1L,
                        execution_unit = "kappa x repetition")
  collect_result <- function() {
    results_marginal <- dplyr::bind_rows(results)
    # Same named objects as H/ICC; ratios are indexed by kappa, then repetition.
    rvec_marg <- setNames(lapply(seq_along(kappa_values), function(ki)
      rvec_output[tasks$kappa_index == ki]), paste0("kappa_", kappa_values))
    list(results_marginal = results_marginal, df_m = results_marginal,
         summary_marginal = summarize_H_PSU_results(results_marginal),
         rvec_marg = rvec_marg, sample_sizes = dplyr::bind_rows(sample_size_output),
         method_diagnostics = method_diagnostic_output,
         diagnostics = dplyr::bind_rows(diagnostic), execution_info = execution_info,
         configuration = configuration)
  }
  completed <- 0L
  log <- function(message) {
    line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message)
    cat(line, "\n"); flush.console()
  }
  log(paste("Output:", output_dir))
  log(sprintf("%d task(s); sequential execution", nrow(tasks)))
  record_result <- function(one, task_id) {
    task <- tasks[task_id, ]
    for (w in one$warnings)
      log(sprintf("WARNING kappa=%s repetition=%d: %s", task$kappa, task$replicate, w))
    if (!is.null(one$error)) {
      message <- sprintf("ERROR kappa=%s repetition=%d: %s", task$kappa, task$replicate, one$error)
      log(message)
      stop(message, call. = FALSE)
    }
    completed <<- completed + 1L
    # Store by task ID, preserving kappa/repetition order in all checkpoints.
    results[[task_id]] <<- one$metric
    diagnostic[[task_id]] <<- one$diagnostic
    rvec_output[[task_id]] <<- one$rvec_marg
    sample_size_output[[task_id]] <<- one$sample_sizes
    method_diagnostic_output[[task_id]] <<- one$method_diagnostics
    if (isTRUE(save_outputs))
      write_H_PSU_outputs(collect_result(), configuration, output_dir,
                           figures_dir, make_figure = FALSE)
    log(sprintf("%d/%d completed; kappa=%s, repetition=%d; AR=%.4f CQR=%.4f; %.1f s",
                completed, nrow(tasks), task$kappa, task$replicate,
                one$metric$coverage[1], one$metric$coverage[2],
                one$elapsed_seconds))
  }
    for (task_id in seq_len(nrow(tasks))) {
      task <- tasks[task_id, ]
      log(sprintf("START kappa=%s, repetition=%d/%d", task$kappa, task$replicate, M))
      record_result(run_H_PSU_task(task, seed), task_id)
    }
  simulation_result <- collect_result()
  if (isTRUE(save_outputs))
    simulation_result$output_paths <- write_H_PSU_outputs(
      simulation_result, configuration, output_dir, figures_dir)
  log("COMPLETE (SD is NA when M=1, since between-repetition SD is undefined).")
  print(simulation_result$summary_marginal, n = Inf)
  invisible(simulation_result)
}


# Production example for the PSU-based design (not executed automatically):
# H_PSU_result <- DA_WCP_simul_H_PSU(M = 100)

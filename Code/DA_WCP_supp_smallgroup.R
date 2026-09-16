if (dir.exists("Code") && file.exists("Code/KLIEP.R")) setwd("Code")
dir.create("../Figures", recursive = TRUE, showWarnings = FALSE)

# Serial supplementary simulations for small demographic subgroups.
# This file combines the baseline and stabilization studies in one sequential script.
# All repetitions, settings, and tuning values are evaluated in their original order.

library(mvtnorm)
library(dplyr)
library(quantreg)
library(ggplot2)

source("KLIEP.R")
source("DA_WCP.R")


###############################################################################
# S1 data-generating mechanism with two mapped demographic subgroups
###############################################################################

smallgroup_settings <- list(
  Balanced = list(q1 = c(0.55, 0.45), q2 = c(0.45, 0.55)),
  Small = list(q1 = c(0.995, 0.005), q2 = c(0.99, 0.01))
)

target_g <- c(0.5, 0.5)
original_subgroup_map <- c(1L, 5L)

var_x <- paste0("X", 1:4)
var_y <- "Y"
var_z <- character(0)
var_g <- "G"
var_h <- character(0) # H = 1

mu_S1 <- function(k) {
  c(k - 3, sin(pi * k / 2), 3 - k, cos(pi * k / 2))
}

cov_S1 <- function(k, rho = 0.2) {
  base_cov <- matrix(nrow = 4, ncol = 4)
  for (i in 1:4) {
    for (j in 1:4) base_cov[i, j] <- rho^(abs(i - j))
  }
  base_cov + diag(rep(0.1 * (k - 1), 4))
}

Delta_mu_S1 <- c(0.3, -0.2, 0.4, 0.1)
Delta_Sigma_S1 <- diag(c(0.05, -0.02, 0.03, 0.01))
beta_S1 <- c(5, -3, 2, 1)
gamma_S1 <- c(-10, -5, 0, 5, 10)
theta_S1 <- c(0.5, -0.4, 0.3, 0.2)
eta_S1 <- c(1.2, 1.6, 2.0, 2.4, 2.8)

stopifnot(
  identical(original_subgroup_map, c(1L, 5L)),
  identical(target_g, c(0.5, 0.5)),
  all(vapply(smallgroup_settings, function(x) sum(x$q1), numeric(1)) == 1),
  all(vapply(smallgroup_settings, function(x) sum(x$q2), numeric(1)) == 1),
  length(var_h) == 0L
)

validate_fixed_n <- function(n, name) {
  if (length(n) != 1L || !is.finite(n) || n < 1L || n != as.integer(n)) {
    stop(name, " must be one positive integer.")
  }
  as.integer(n)
}

generate_smallgroup_data <- function(n, t, subgroup_probabilities = NULL,
                                     test = FALSE) {
  n <- validate_fixed_n(n, "n")
  if (!t %in% c(1L, 2L)) stop("t must be 1 or 2.")
  probabilities <- if (test) target_g else subgroup_probabilities
  if (is.null(probabilities) || length(probabilities) != 2L ||
      any(probabilities < 0) || !isTRUE(all.equal(sum(probabilities), 1))) {
    stop("Subgroup probabilities must be a length-two probability vector.")
  }

  G <- sample(1:2, size = n, replace = TRUE, prob = probabilities)
  original_k <- original_subgroup_map[G]
  X <- matrix(0, nrow = n, ncol = 4)
  Y <- numeric(n)

  for (g in 1:2) {
    idx <- which(G == g)
    if (length(idx) == 0L) next
    k <- original_subgroup_map[g]
    mu_k <- mu_S1(k) + if (t == 2L) Delta_mu_S1 else 0
    Sigma_k <- cov_S1(k) + if (t == 2L) Delta_Sigma_S1 else 0
    X[idx, ] <- MASS::mvrnorm(length(idx), mu = mu_k, Sigma = Sigma_k)
    outcome_mean <- drop(X[idx, , drop = FALSE] %*% beta_S1) + gamma_S1[k]
    outcome_sd <- sqrt(exp(
      drop(X[idx, , drop = FALSE] %*% theta_S1) + eta_S1[k]
    ))
    Y[idx] <- outcome_mean + rnorm(length(idx), sd = outcome_sd)
  }

  data <- as.data.frame(X)
  names(data) <- var_x
  data$G <- G
  if (t == 1L || isTRUE(test)) data$Y <- Y
  data$time <- t
  attr(data, "original_subgroup") <- original_k
  stopifnot(nrow(data) == n)
  data
}

group_counts <- function(G) {
  as.integer(table(factor(G, levels = 1:2)))
}

primary_weight_diagnostic_metrics <- c(
  "kliep_mae_g1",
  "kliep_mae_g2",
  "kliep_mae_all",
  "kliep_medae_g1",
  "kliep_medae_g2",
  "kliep_medae_all",
  "subgroup_ratio_abs_error_g1",
  "subgroup_ratio_abs_error_g2",
  "wcp_weight_mae_g1",
  "wcp_weight_mae_g2",
  "wcp_weight_mae_all",
  "wcp_weight_medae_g1",
  "wcp_weight_medae_g2",
  "wcp_weight_medae_all"
)

true_KLIEP_ratio <- function(calibration) {
  true_ratio <- rep(NA_real_, nrow(calibration))
  X_calibration <- as.matrix(calibration[, var_x, drop = FALSE])

  for (g in 1:2) {
    index <- which(calibration$G == g)
    if (length(index) == 0L) next
    k <- original_subgroup_map[g]
    log_density_historical <- mvtnorm::dmvnorm(
      X_calibration[index, , drop = FALSE],
      mean = mu_S1(k), sigma = cov_S1(k), log = TRUE
    )
    log_density_current <- mvtnorm::dmvnorm(
      X_calibration[index, , drop = FALSE],
      mean = mu_S1(k) + Delta_mu_S1,
      sigma = cov_S1(k) + Delta_Sigma_S1, log = TRUE
    )
    true_ratio[index] <- exp(log_density_current - log_density_historical)
  }

  if (length(true_ratio) != nrow(calibration) ||
      any(!is.finite(true_ratio)) || any(true_ratio <= 0)) {
    stop("Analytic true KLIEP ratios must be finite and strictly positive.")
  }
  true_ratio
}

make_smallgroup_splits <- function(data1, repetition) {
  n <- nrow(data1)

  set.seed(100000L + as.integer(repetition))
  marginal_train <- sort(sample.int(n, n %/% 2L))
  marginal_cal <- setdiff(seq_len(n), marginal_train)

  stopifnot(
    length(intersect(marginal_train, marginal_cal)) == 0L,
    setequal(c(marginal_train, marginal_cal), seq_len(n))
  )

  list(
    marginal_train = marginal_train,
    marginal_cal = marginal_cal
  )
}


###############################################################################
# Diagnostics and failure-tolerant DA-WCP evaluation
###############################################################################

make_sample_diagnostics <- function(repetition, setting, specification,
                                    data1, data2, splits) {
  n1_counts <- group_counts(data1$G)
  n2_counts <- group_counts(data2$G)
  marginal_train_counts <- group_counts(data1$G[splits$marginal_train])
  marginal_cal_counts <- group_counts(data1$G[splits$marginal_cal])

  data.frame(
    repetition = repetition,
    setting = setting,
    q1_subgroup_1 = specification$q1[1L],
    q1_subgroup_2 = specification$q1[2L],
    q2_subgroup_1 = specification$q2[1L],
    q2_subgroup_2 = specification$q2[2L],
    historical_total = nrow(data1),
    current_total = nrow(data2),
    historical_subgroup_1_count = n1_counts[1L],
    historical_subgroup_2_count = n1_counts[2L],
    current_subgroup_1_count = n2_counts[1L],
    current_subgroup_2_count = n2_counts[2L],
    historical_subgroup_1_proportion = n1_counts[1L] / nrow(data1),
    historical_subgroup_2_proportion = n1_counts[2L] / nrow(data1),
    current_subgroup_1_proportion = n2_counts[1L] / nrow(data2),
    current_subgroup_2_proportion = n2_counts[2L] / nrow(data2),
    historical_subgroup_1_training_count = marginal_train_counts[1L],
    historical_subgroup_1_calibration_count = marginal_cal_counts[1L],
    historical_subgroup_2_training_count = marginal_train_counts[2L],
    historical_subgroup_2_calibration_count = marginal_cal_counts[2L],
    historical_subgroup_2_zero = n1_counts[2L] == 0L,
    current_subgroup_2_zero = n2_counts[2L] == 0L,
    stringsAsFactors = FALSE
  )
}

empty_sample_diagnostics <- function(repetition, setting, specification) {
  template <- make_sample_diagnostics(
    repetition = repetition,
    setting = setting,
    specification = specification,
    data1 = data.frame(G = c(1L, 2L)),
    data2 = data.frame(G = c(1L, 2L)),
    splits = list(
      marginal_train = 1L, marginal_cal = 2L
    )
  )
  count_columns <- setdiff(
    names(template),
    c("repetition", "setting", "q1_subgroup_1", "q1_subgroup_2",
      "q2_subgroup_1", "q2_subgroup_2")
  )
  template[count_columns] <- NA
  template
}

safe_compute <- function(expression) {
  warning_messages <- character(0)
  tryCatch({
    value <- withCallingHandlers(
      force(expression),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(
      ok = TRUE, value = value, error = NA_character_,
      warnings = paste(unique(warning_messages), collapse = " | ")
    )
  }, error = function(e) {
    list(
      ok = FALSE, value = NULL, error = conditionMessage(e),
      warnings = paste(unique(warning_messages), collapse = " | ")
    )
  })
}

empty_weight_diagnostics <- function(repetition, setting, specification,
                                     training = NULL, calibration = NULL,
                                     reason = NA_character_) {
  q1_true <- specification$q1
  subgroup_ratio_true <- target_g / q1_true
  training_counts <- if (is.null(training)) c(NA_integer_, NA_integer_) else
    group_counts(training$G)
  calibration_counts <- if (is.null(calibration)) {
    c(NA_integer_, NA_integer_)
  } else {
    group_counts(calibration$G)
  }
  q1_hat <- if (is.null(training) || nrow(training) == 0L) {
    c(NA_real_, NA_real_)
  } else {
    training_counts / nrow(training)
  }
  subgroup_ratio_hat <- ifelse(
    is.finite(q1_hat) & q1_hat > 0, target_g / q1_hat, NA_real_
  )

  output <- data.frame(
    repetition = repetition,
    setting = setting,
    kliep_mae_g1 = NA_real_,
    kliep_mae_g2 = NA_real_,
    kliep_mae_all = NA_real_,
    kliep_medae_g1 = NA_real_,
    kliep_medae_g2 = NA_real_,
    kliep_medae_all = NA_real_,
    subgroup_ratio_abs_error_g1 = if (is.finite(subgroup_ratio_hat[1L])) {
      abs(subgroup_ratio_hat[1L] - subgroup_ratio_true[1L])
    } else NA_real_,
    subgroup_ratio_abs_error_g2 = if (is.finite(subgroup_ratio_hat[2L])) {
      abs(subgroup_ratio_hat[2L] - subgroup_ratio_true[2L])
    } else NA_real_,
    wcp_weight_mae_g1 = NA_real_,
    wcp_weight_mae_g2 = NA_real_,
    wcp_weight_mae_all = NA_real_,
    wcp_weight_medae_g1 = NA_real_,
    wcp_weight_medae_g2 = NA_real_,
    wcp_weight_medae_all = NA_real_,
    q1_true_g1 = q1_true[1L],
    q1_true_g2 = q1_true[2L],
    q1_hat_g1 = q1_hat[1L],
    q1_hat_g2 = q1_hat[2L],
    subgroup_ratio_true_g1 = subgroup_ratio_true[1L],
    subgroup_ratio_true_g2 = subgroup_ratio_true[2L],
    subgroup_ratio_hat_g1 = subgroup_ratio_hat[1L],
    subgroup_ratio_hat_g2 = subgroup_ratio_hat[2L],
    training_g1_count = training_counts[1L],
    training_g2_count = training_counts[2L],
    calibration_g1_count = calibration_counts[1L],
    calibration_g2_count = calibration_counts[2L],
    true_ratio_finite_positive = NA,
    calibration_row_alignment_verified = NA,
    same_rvec_reused_by_ar_cqr = NA,
    raw_weight_formula_verified = NA,
    raw_weight_scale = "before proportions() normalization",
    diagnostic_ok = FALSE,
    diagnostic_failure_reason = reason,
    stringsAsFactors = FALSE
  )
  output
}

compute_weight_diagnostics <- function(repetition, setting, specification,
                                       training, calibration, r_hat,
                                       same_rvec_reused_by_ar_cqr) {
  output <- empty_weight_diagnostics(
    repetition, setting, specification, training, calibration
  )

  calculation <- safe_compute({
    q1_hat <- c(output$q1_hat_g1, output$q1_hat_g2)
    if (any(!is.finite(q1_hat)) || any(q1_hat <= 0)) {
      stop("A training-sample q1_hat is zero or non-finite.")
    }
    calibration_counts <- c(
      output$calibration_g1_count, output$calibration_g2_count
    )
    if (any(calibration_counts == 0L)) {
      stop("A calibration subgroup count is zero.")
    }
    if (is.null(r_hat) || length(r_hat) != nrow(calibration)) {
      stop("Estimated KLIEP ratio length does not match calibration rows.")
    }
    if (any(!is.finite(r_hat))) {
      stop("Estimated KLIEP ratios contain non-finite values.")
    }

    r_true <- true_KLIEP_ratio(calibration)
    if (length(r_true) != length(r_hat)) {
      stop("True and estimated KLIEP ratios are not row-aligned.")
    }

    subgroup_ratio_true <- c(
      output$subgroup_ratio_true_g1, output$subgroup_ratio_true_g2
    )
    subgroup_ratio_hat <- target_g / q1_hat
    estimated_raw_weight <- r_hat * subgroup_ratio_hat[calibration$G]
    true_raw_weight <- r_true * subgroup_ratio_true[calibration$G]

    if (any(!is.finite(estimated_raw_weight)) ||
        any(!is.finite(true_raw_weight))) {
      stop("Raw estimated or true WCP weights contain non-finite values.")
    }

    absolute_KLIEP_error <- abs(r_hat - r_true)
    absolute_weight_error <- abs(estimated_raw_weight - true_raw_weight)
    group_mae <- function(absolute_error, g) {
      mean(absolute_error[calibration$G == g])
    }
    group_medae <- function(absolute_error, g) {
      median(absolute_error[calibration$G == g])
    }
    metrics <- c(
      kliep_mae_g1 = group_mae(absolute_KLIEP_error, 1L),
      kliep_mae_g2 = group_mae(absolute_KLIEP_error, 2L),
      kliep_mae_all = mean(absolute_KLIEP_error),
      kliep_medae_g1 = group_medae(absolute_KLIEP_error, 1L),
      kliep_medae_g2 = group_medae(absolute_KLIEP_error, 2L),
      kliep_medae_all = median(absolute_KLIEP_error),
      subgroup_ratio_abs_error_g1 = abs(
        subgroup_ratio_hat[1L] - subgroup_ratio_true[1L]
      ),
      subgroup_ratio_abs_error_g2 = abs(
        subgroup_ratio_hat[2L] - subgroup_ratio_true[2L]
      ),
      wcp_weight_mae_g1 = group_mae(absolute_weight_error, 1L),
      wcp_weight_mae_g2 = group_mae(absolute_weight_error, 2L),
      wcp_weight_mae_all = mean(absolute_weight_error),
      wcp_weight_medae_g1 = group_medae(absolute_weight_error, 1L),
      wcp_weight_medae_g2 = group_medae(absolute_weight_error, 2L),
      wcp_weight_medae_all = median(absolute_weight_error)
    )

    if (any(!is.finite(metrics)) || any(metrics < 0)) {
      stop("Primary weight diagnostics must be finite and nonnegative.")
    }
    pooled_KLIEP_mae <- sum(
      calibration_counts * metrics[c("kliep_mae_g1", "kliep_mae_g2")]
    ) / sum(calibration_counts)
    pooled_weight_mae <- sum(
      calibration_counts * metrics[c(
        "wcp_weight_mae_g1", "wcp_weight_mae_g2"
      )]
    ) / sum(calibration_counts)
    if (!isTRUE(all.equal(
      unname(metrics["kliep_mae_all"]), unname(pooled_KLIEP_mae),
      tolerance = 1e-12
    )) || !isTRUE(all.equal(
      unname(metrics["wcp_weight_mae_all"]), unname(pooled_weight_mae),
      tolerance = 1e-12
    ))) {
      stop("Overall MAE is not the observation-level pooled MAE.")
    }

    list(
      metrics = metrics,
      subgroup_ratio_hat = subgroup_ratio_hat,
      true_ratio_finite_positive = all(is.finite(r_true) & r_true > 0),
      row_alignment_verified = length(r_hat) == nrow(calibration),
      raw_weight_formula_verified = isTRUE(all.equal(
        estimated_raw_weight,
        r_hat * subgroup_ratio_hat[calibration$G], tolerance = 0
      )) && isTRUE(all.equal(
        true_raw_weight,
        r_true * subgroup_ratio_true[calibration$G], tolerance = 0
      ))
    )
  })

  if (!calculation$ok) {
    output$diagnostic_failure_reason <- calculation$error
    return(output)
  }

  output[primary_weight_diagnostic_metrics] <- as.list(
    calculation$value$metrics[primary_weight_diagnostic_metrics]
  )
  output$subgroup_ratio_hat_g1 <- calculation$value$subgroup_ratio_hat[1L]
  output$subgroup_ratio_hat_g2 <- calculation$value$subgroup_ratio_hat[2L]
  output$true_ratio_finite_positive <-
    calculation$value$true_ratio_finite_positive
  output$calibration_row_alignment_verified <-
    calculation$value$row_alignment_verified
  output$same_rvec_reused_by_ar_cqr <- same_rvec_reused_by_ar_cqr
  output$raw_weight_formula_verified <-
    calculation$value$raw_weight_formula_verified
  output$diagnostic_ok <- TRUE
  output$diagnostic_failure_reason <- NA_character_
  output
}

validate_DA_inputs <- function(train, calibration, current, target,
                               context) {
  counts <- rbind(
    training = group_counts(train$G),
    calibration = group_counts(calibration$G),
    current = group_counts(current$G),
    target = group_counts(target$G)
  )
  empty_cells <- which(counts == 0L, arr.ind = TRUE)
  if (nrow(empty_cells) > 0L) {
    labels <- apply(empty_cells, 1L, function(index) {
      paste0(rownames(counts)[index[1L]], " subgroup ", index[2L])
    })
    stop(context, " cannot run because these samples are empty: ",
         paste(labels, collapse = ", "))
  }
  invisible(TRUE)
}

fit_DA_WCP_safely <- function(train, calibration, current, target,
                              val_type, score_type,
                              alpha, n_alpha, R,
                              sig_wo_cv, kliep_cv,
                              cached_rvec = NULL) {
  safe_compute({
    validate_DA_inputs(
      train, calibration, current, target,
      paste("DA-WCP", val_type, score_type)
    )
    use_cached <- !is.null(cached_rvec)
    if (use_cached && length(cached_rvec) != nrow(calibration)) {
      stop("Cached KLIEP ratio length does not match calibration sample size.")
    }
    DA_WCP(
      train, calibration, train, current, target,
      var_x, var_y, var_z, var_g, var_h,
      gwcp = FALSE, val_type = val_type, score_type = score_type,
      p2_vec = target_g,
      alpha = alpha, n_alpha = n_alpha,
      rvec = use_cached,
      rvec_cal = if (use_cached) cached_rvec else NA,
      R = R, sig_wo_cv = sig_wo_cv,
      kliep_cv = kliep_cv, simul = TRUE
    )
  })
}

fit_score_pair <- function(train, calibration, current, target,
                           val_type, alpha, n_alpha, R,
                           sig_wo_cv, kliep_cv,
                           cached_rvec = NULL) {
  ar_fit <- fit_DA_WCP_safely(
    train, calibration, current, target,
    val_type = val_type, score_type = "ar",
    alpha = alpha, n_alpha = n_alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
    cached_rvec = cached_rvec
  )

  reusable_rvec <- if (ar_fit$ok) {
    ar_fit$value$rvec_cal
  } else {
    cached_rvec
  }
  cqr_fit <- fit_DA_WCP_safely(
    train, calibration, current, target,
    val_type = val_type, score_type = "cqr",
    alpha = alpha, n_alpha = n_alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
    cached_rvec = reusable_rvec
  )

  if (is.null(reusable_rvec) && cqr_fit$ok) {
    reusable_rvec <- cqr_fit$value$rvec_cal
  }
  same_rvec_reused_by_ar_cqr <- if (ar_fit$ok && cqr_fit$ok) {
    isTRUE(all.equal(
      ar_fit$value$rvec_cal, cqr_fit$value$rvec_cal, tolerance = 0
    ))
  } else {
    NA
  }
  list(
    AR = ar_fit, CQR = cqr_fit, rvec = reusable_rvec,
    same_rvec_reused_by_ar_cqr = same_rvec_reused_by_ar_cqr
  )
}

interval_metrics <- function(result_mat, target, subgroup = NULL) {
  index <- if (is.null(subgroup)) {
    seq_len(nrow(target))
  } else {
    which(target$G == subgroup)
  }
  if (length(index) == 0L) stop("The requested target subgroup is empty.")
  endpoints <- as.matrix(result_mat[index, 1:2, drop = FALSE])
  if (nrow(endpoints) != length(index) || any(!is.finite(endpoints))) {
    stop("Prediction interval endpoints are missing or non-finite.")
  }
  covered <- target$Y[index] >= endpoints[, 1L] &
    target$Y[index] <= endpoints[, 2L]
  interval_length <- endpoints[, 2L] - endpoints[, 1L]
  c(
    empirical_coverage = mean(covered),
    average_prediction_length = mean(interval_length)
  )
}

result_row_from_fit <- function(fit, repetition, setting, score,
                                target_type, target, diagnostics,
                                subgroup = NULL) {
  metric_result <- if (fit$ok) {
    safe_compute(interval_metrics(fit$value$result_mat, target, subgroup))
  } else {
    list(ok = FALSE, value = NULL, error = fit$error, warnings = fit$warnings)
  }
  computation_ok <- isTRUE(metric_result$ok)
  metrics <- if (computation_ok) metric_result$value else
    c(empirical_coverage = NA_real_, average_prediction_length = NA_real_)
  warning_message <- paste(
    unique(c(fit$warnings, metric_result$warnings)), collapse = " | "
  )
  warning_message <- gsub("^( \\| )+|( \\| )+$", "", warning_message)

  row <- data.frame(
    repetition = repetition,
    setting = setting,
    method = "DA_WCP",
    score = score,
    target_type = target_type,
    empirical_coverage = unname(metrics["empirical_coverage"]),
    average_prediction_length = unname(metrics["average_prediction_length"]),
    computation_ok = computation_ok,
    failure_indicator = !computation_ok,
    failure_reason = if (computation_ok) NA_character_ else metric_result$error,
    warning_message = if (nzchar(warning_message)) warning_message else NA_character_,
    stringsAsFactors = FALSE
  )

  diagnostic_columns <- setdiff(names(diagnostics), c("repetition", "setting"))
  cbind(row, diagnostics[rep(1L, nrow(row)), diagnostic_columns, drop = FALSE])
}

failure_rows_for_setting <- function(repetition, setting, reason, diagnostics) {
  design <- data.frame(
    score = c("AR", "CQR"),
    target_type = c("marginal", "marginal"),
    stringsAsFactors = FALSE
  )
  rows <- data.frame(
    repetition = repetition,
    setting = setting,
    method = "DA_WCP",
    score = design$score,
    target_type = design$target_type,
    empirical_coverage = NA_real_,
    average_prediction_length = NA_real_,
    computation_ok = FALSE,
    failure_indicator = TRUE,
    failure_reason = reason,
    warning_message = NA_character_,
    stringsAsFactors = FALSE
  )
  diagnostic_columns <- setdiff(names(diagnostics), c("repetition", "setting"))
  cbind(rows, diagnostics[rep(1L, nrow(rows)), diagnostic_columns, drop = FALSE])
}

run_smallgroup_setting <- function(repetition, setting, specification,
                                   target, n1, n2,
                                   alpha, n_alpha, R,
                                   sig_wo_cv, kliep_cv,
                                   cached_rvec_marg = NULL) {
  cat("Setting: ", setting, "\n", sep = "")

  # Resetting these seeds for each setting supplies common random numbers;
  # only q1 and q2 differ between Balanced and Small.
  set.seed(123L + as.integer(repetition))
  data1 <- generate_smallgroup_data(n1, t = 1L, specification$q1)
  set.seed(100123L + as.integer(repetition))
  data2 <- generate_smallgroup_data(n2, t = 2L, specification$q2)
  splits <- make_smallgroup_splits(data1, repetition)
  diagnostics <- make_sample_diagnostics(
    repetition, setting, specification, data1, data2, splits
  )

  data1_y <- data1[, c(var_x, var_y, var_g), drop = FALSE]
  data2_x <- data2[, c(var_x, var_g), drop = FALSE]
  target_y <- target[, c(var_x, var_y, var_g), drop = FALSE]

  marginal_train <- data1_y[splits$marginal_train, , drop = FALSE]
  marginal_cal <- data1_y[splits$marginal_cal, , drop = FALSE]

  marginal_fits <- fit_score_pair(
    marginal_train, marginal_cal, data2_x, target_y,
    val_type = "marg", alpha = alpha, n_alpha = n_alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
    cached_rvec = cached_rvec_marg
  )
  weight_diagnostics <- compute_weight_diagnostics(
    repetition = repetition,
    setting = setting,
    specification = specification,
    training = marginal_train,
    calibration = marginal_cal,
    r_hat = marginal_fits$rvec,
    same_rvec_reused_by_ar_cqr =
      marginal_fits$same_rvec_reused_by_ar_cqr
  )
  results <- bind_rows(
    result_row_from_fit(
      marginal_fits$AR, repetition, setting, "AR", "marginal",
      target_y, diagnostics
    ),
    result_row_from_fit(
      marginal_fits$CQR, repetition, setting, "CQR", "marginal",
      target_y, diagnostics
    )
  )

  if (any(results$failure_indicator)) {
    failed <- results[results$failure_indicator,
                      c("score", "target_type", "failure_reason"), drop = FALSE]
    for (i in seq_len(nrow(failed))) {
      cat(
        "Failure [", setting, ", ", failed$score[i], ", ",
        failed$target_type[i], "]: ", failed$failure_reason[i], "\n",
        sep = ""
      )
    }
  }

  list(
    results = results,
    diagnostics = diagnostics,
    weight_diagnostics = weight_diagnostics,
    rvec_marg = marginal_fits$rvec
  )
}

run_smallgroup_replicate <- function(repetition, n1, n2, n_test,
                                     alpha, n_alpha, R,
                                     sig_wo_cv, kliep_cv,
                                     cached_rvec_marg = NULL) {
  cat("------------------------------------------------------\n")
  cat("     Small-subgroup sensitivity - replicate ", repetition, "\n",
      sep = "")
  cat("------------------------------------------------------\n")

  set.seed(200123L + as.integer(repetition))
  target <- generate_smallgroup_data(n_test, t = 2L, test = TRUE)

  setting_results <- vector("list", length(smallgroup_settings))
  setting_diagnostics <- vector("list", length(smallgroup_settings))
  setting_weight_diagnostics <- vector("list", length(smallgroup_settings))
  rvec_marg_output <- list()

  setting_names <- names(smallgroup_settings)
  for (setting_index in seq_along(setting_names)) {
    setting <- setting_names[setting_index]
    specification <- smallgroup_settings[[setting]]
    cached_m <- if (is.null(cached_rvec_marg)) NULL else
      cached_rvec_marg[[setting]]

    one_setting <- tryCatch(
      run_smallgroup_setting(
        repetition, setting, specification, target,
        n1, n2, alpha, n_alpha, R, sig_wo_cv, kliep_cv,
        cached_rvec_marg = cached_m
      ),
      error = function(e) {
        diagnostics <- empty_sample_diagnostics(
          repetition, setting, specification
        )
        list(
          results = failure_rows_for_setting(
            repetition, setting,
            paste("Setting-level failure:", conditionMessage(e)),
            diagnostics
          ),
          diagnostics = diagnostics,
          weight_diagnostics = empty_weight_diagnostics(
            repetition, setting, specification,
            reason = paste("Setting-level failure:", conditionMessage(e))
          ),
          rvec_marg = NULL
        )
      }
    )

    setting_results[[setting_index]] <- one_setting$results
    setting_diagnostics[[setting_index]] <- one_setting$diagnostics
    setting_weight_diagnostics[[setting_index]] <-
      one_setting$weight_diagnostics
    rvec_marg_output[[setting]] <- one_setting$rvec_marg
  }

  list(
    results = bind_rows(setting_results),
    sample_diagnostics = bind_rows(setting_diagnostics),
    weight_diagnostics = bind_rows(setting_weight_diagnostics),
    rvec_marg = rvec_marg_output
  )
}


###############################################################################
# Sequential simulation driver
###############################################################################

# Runs the baseline small-subgroup experiment sequentially.
DA_WCP_simul_smallgroup_1 <- function(
    n1 = 5000, n2 = 5000, n_test = 5000,
    M = 100, alpha = 0.2, n_alpha = 5, R = 4,
    sig_wo_cv = 1.2, kliep_cv = FALSE,
    rvec = FALSE, rvec_marg = NULL,
    save_outputs = TRUE,
    output_dir = "../Figures",
    figures_dir = "../Figures") {
  n1 <- validate_fixed_n(n1, "n1")
  n2 <- validate_fixed_n(n2, "n2")
  n_test <- validate_fixed_n(n_test, "n_test")
  M <- validate_fixed_n(M, "M")
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1).")
  }
  if (rvec && is.null(rvec_marg)) {
    stop("rvec_marg is required when rvec = TRUE.")
  }

  tasks <- lapply(seq_len(M), function(repetition) {
    list(
      repetition = repetition,
      cached_rvec_marg = if (rvec) rvec_marg[[repetition]] else NULL
    )
  })
  common_args <- list(
    n1 = n1, n2 = n2, n_test = n_test,
    alpha = alpha, n_alpha = n_alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv
  )

    replicate_results <- lapply(
      tasks,
      function(task) do.call(run_smallgroup_replicate, c(task, common_args))
    )
    execution_info <- list(
      backend = "sequential",
      cores = 1L,
      details = list(),
      threads = 1L
    )

  results <- bind_rows(lapply(replicate_results, `[[`, "results"))
  sample_size_diagnostics <- bind_rows(lapply(
    replicate_results, `[[`, "sample_diagnostics"
  ))
  weight_diagnostics <- bind_rows(lapply(
    replicate_results, `[[`, "weight_diagnostics"
  ))
  rvec_marg_output <- lapply(replicate_results, `[[`, "rvec_marg")

  results$setting <- factor(
    results$setting, levels = c("Balanced", "Small")
  )
  results$score <- factor(results$score, levels = c("AR", "CQR"))
  results$target_type <- factor(
    results$target_type, levels = "marginal"
  )
  sample_size_diagnostics$setting <- factor(
    sample_size_diagnostics$setting, levels = c("Balanced", "Small")
  )
  weight_diagnostics$setting <- factor(
    weight_diagnostics$setting, levels = c("Balanced", "Small")
  )

  stopifnot(
    nrow(results) == M * length(smallgroup_settings) * 2L,
    nrow(sample_size_diagnostics) == M * length(smallgroup_settings),
    nrow(weight_diagnostics) == M * length(smallgroup_settings),
    all(sample_size_diagnostics$historical_total == n1),
    all(sample_size_diagnostics$current_total == n2)
  )

  configuration <- list(
    scenario = "S1 small-demographic-subgroup sensitivity",
    settings = smallgroup_settings,
    only_setting_difference = c("q1", "q2"),
    original_subgroup_map = original_subgroup_map,
    target_g = target_g,
    n1 = n1, n2 = n2, n_test = n_test,
    fixed_sample_sizes = TRUE,
    poisson_sample_size = FALSE,
    M = M, H = 1L,
    within_cluster_dependence = FALSE,
    alpha = alpha, n_alpha = n_alpha,
    R_kliep = R, sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
    methods = "DA_WCP",
    scores = c("AR", "CQR"),
    target_types = "marginal",
    marginal_split = "individual-level random 50:50",
    primary_weight_diagnostic_metrics = primary_weight_diagnostic_metrics,
    true_KLIEP_ratio =
      "f_X2_given_G(x) / f_X1_given_G(x), evaluated by Gaussian log densities",
    subgroup_ratio_hat_source = "same marginal DA-WCP training sample",
    raw_weight_diagnostic_scale = "before proportions() normalization",
    execution_backend = execution_info$backend,
    execution_cores = execution_info$cores,
    threads = 1L
  )

  simulation_result <- list(
    results = results,
    sample_size_diagnostics = sample_size_diagnostics,
    weight_diagnostics = weight_diagnostics,
    rvec_marg = rvec_marg_output,
    execution_info = execution_info,
    configuration = configuration
  )

  if (isTRUE(save_outputs)) {
    simulation_result$output_paths <- write_smallgroup_1_outputs(
      simulation_result, output_dir, figures_dir
    )
  }
  simulation_result
}


###############################################################################
# Summaries, output files, and coverage figure
###############################################################################

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) NA_real_ else stats::sd(x)
}

summarize_smallgroup_results <- function(results) {
  results %>%
    group_by(setting, method, score, target_type) %>%
    summarise(
      mean_empirical_coverage = mean_or_na(empirical_coverage),
      sd_empirical_coverage = sd_or_na(empirical_coverage),
      monte_carlo_se = ifelse(
        sum(computation_ok & is.finite(empirical_coverage)) > 1L,
        sd_or_na(empirical_coverage) /
          sqrt(sum(computation_ok & is.finite(empirical_coverage))),
        NA_real_
      ),
      mean_prediction_length = mean_or_na(average_prediction_length),
      sd_prediction_length = sd_or_na(average_prediction_length),
      n_valid_repetitions = sum(
        computation_ok & is.finite(empirical_coverage) &
          is.finite(average_prediction_length)
      ),
      n_failed_repetitions = sum(failure_indicator),
      .groups = "drop"
    )
}

summarize_sample_sizes <- function(sample_size_diagnostics, results) {
  failure_by_repetition <- results %>%
    group_by(repetition, setting) %>%
    summarise(
      any_computation_failure = any(failure_indicator),
      .groups = "drop"
    )

  diagnostics <- sample_size_diagnostics %>%
    left_join(failure_by_repetition, by = c("repetition", "setting"))

  diagnostics %>%
    group_by(setting) %>%
    summarise(
      mean_historical_subgroup_2_count = mean_or_na(
        historical_subgroup_2_count
      ),
      min_historical_subgroup_2_count = ifelse(
        all(is.na(historical_subgroup_2_count)), NA_real_,
        min(historical_subgroup_2_count, na.rm = TRUE)
      ),
      max_historical_subgroup_2_count = ifelse(
        all(is.na(historical_subgroup_2_count)), NA_real_,
        max(historical_subgroup_2_count, na.rm = TRUE)
      ),
      mean_historical_subgroup_2_training_count = mean_or_na(
        historical_subgroup_2_training_count
      ),
      mean_historical_subgroup_2_calibration_count = mean_or_na(
        historical_subgroup_2_calibration_count
      ),
      mean_current_subgroup_2_count = mean_or_na(
        current_subgroup_2_count
      ),
      historical_zero_count_frequency = mean(
        historical_subgroup_2_zero, na.rm = TRUE
      ),
      current_zero_count_frequency = mean(
        current_subgroup_2_zero, na.rm = TRUE
      ),
      any_zero_count_frequency = mean(
        historical_subgroup_2_zero | current_subgroup_2_zero,
        na.rm = TRUE
      ),
      computation_failure_frequency = mean(
        any_computation_failure, na.rm = TRUE
      ),
      .groups = "drop"
    )
}

summarize_weight_diagnostics <- function(weight_diagnostics) {
  setting_names <- c("Balanced", "Small")
  summary_rows <- lapply(setting_names, function(setting_name) {
    setting_data <- weight_diagnostics[
      as.character(weight_diagnostics$setting) == setting_name, , drop = FALSE
    ]
    output <- data.frame(setting = setting_name, stringsAsFactors = FALSE)
    for (metric in primary_weight_diagnostic_metrics) {
      values <- setting_data[[metric]]
      output[[paste0("mean_", metric)]] <- mean_or_na(values)
      output[[paste0("sd_", metric)]] <- sd_or_na(values)
      output[[paste0("median_", metric)]] <- stats::median(
        values[is.finite(values)]
      )
      output[[paste0("n_valid_", metric)]] <- sum(is.finite(values))
    }
    output
  })
  bind_rows(summary_rows)
}

save_pdf_plot <- function(plot_object, path, width, height) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_object)
}

make_smallgroup_1_coverage_figure <- function(results, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  plot_data <- results[results$computation_ok, , drop = FALSE]
  plot_data$setting <- factor(
    plot_data$setting, levels = c("Balanced", "Small")
  )
  plot_data$score <- factor(plot_data$score, levels = c("AR", "CQR"))

  coverage_plot <- ggplot(
    plot_data,
    aes(x = setting, y = empirical_coverage, fill = score)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.65, outlier.size = 0.7
    ) +
    geom_hline(
      yintercept = 0.8, color = "red", linetype = "dashed", linewidth = 0.7
    ) +
    scale_fill_manual(values = c(AR = "#9ecae1", CQR = "#3182bd")) +
    labs(
      x = "Survey subgroup composition",
      y = "Empirical coverage",
      fill = "Score"
    ) +
    theme_test(base_size = 13) +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )

  figure_path <- file.path(
    figures_dir, "simulation_smallgroup_1_coverage.pdf"
  )
  save_pdf_plot(coverage_plot, figure_path, width = 7, height = 4.5)
  figure_path
}

build_smallgroup_1_all_information <- function(simulation_result) {
  results <- simulation_result$results
  performance_summary <- summarize_smallgroup_results(results)
  sample_size_summary <- summarize_sample_sizes(
    simulation_result$sample_size_diagnostics, results
  )
  weight_diagnostics_summary <- summarize_weight_diagnostics(
    simulation_result$weight_diagnostics
  )

  results %>%
    left_join(
      simulation_result$weight_diagnostics,
      by = c("repetition", "setting")
    ) %>%
    left_join(
      performance_summary,
      by = c("setting", "method", "score", "target_type")
    ) %>%
    left_join(sample_size_summary, by = "setting") %>%
    left_join(weight_diagnostics_summary, by = "setting") %>%
    arrange(repetition, setting, score)
}

write_smallgroup_1_outputs <- function(simulation_result,
                                       output_dir, figures_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  all_information <- build_smallgroup_1_all_information(simulation_result)

  csv_path <- file.path(
    output_dir, "simulation_smallgroup_1_results.csv"
  )
  write.csv(all_information, csv_path, row.names = FALSE)
  figure_path <- make_smallgroup_1_coverage_figure(
    simulation_result$results, figures_dir
  )

  output_paths <- c(csv_path, figure_path)
  stopifnot(length(output_paths) == 2L, all(file.exists(output_paths)))
  output_paths
}


###############################################################################
# Fast full-n sanity option (M = 2 by default)
###############################################################################

run_smallgroup_1_sanity_check <- function(M = 2, save_outputs = FALSE) {
  sanity_result <- DA_WCP_simul_smallgroup_1(
    n1 = 5000, n2 = 5000, n_test = 5000,
    M = M,
    save_outputs = save_outputs
  )
  diagnostics <- sanity_result$sample_size_diagnostics
  weight_diagnostics <- sanity_result$weight_diagnostics

  mean_diagnostics <- diagnostics %>%
    group_by(setting) %>%
    summarise(
      historical_subgroup_1_proportion = mean(
        historical_subgroup_1_proportion
      ),
      historical_subgroup_2_proportion = mean(
        historical_subgroup_2_proportion
      ),
      current_subgroup_1_proportion = mean(current_subgroup_1_proportion),
      current_subgroup_2_proportion = mean(current_subgroup_2_proportion),
      mean_historical_subgroup_2_count = mean(
        historical_subgroup_2_count
      ),
      mean_current_subgroup_2_count = mean(current_subgroup_2_count),
      .groups = "drop"
    )

  balanced <- mean_diagnostics[
    as.character(mean_diagnostics$setting) == "Balanced", , drop = FALSE
  ]
  small <- mean_diagnostics[
    as.character(mean_diagnostics$setting) == "Small", , drop = FALSE
  ]

  valid_weight_diagnostics <- weight_diagnostics[
    weight_diagnostics$diagnostic_ok, , drop = FALSE
  ]
  metric_matrix <- as.matrix(
    valid_weight_diagnostics[, primary_weight_diagnostic_metrics, drop = FALSE]
  )
  q1_hat_from_training <- with(
    valid_weight_diagnostics,
    abs(q1_hat_g1 - training_g1_count /
          (training_g1_count + training_g2_count)) < 1e-15 &
      abs(q1_hat_g2 - training_g2_count /
          (training_g1_count + training_g2_count)) < 1e-15
  )
  calibration_total <- with(
    valid_weight_diagnostics,
    calibration_g1_count + calibration_g2_count
  )
  KLIEP_pooled <- with(
    valid_weight_diagnostics,
    (calibration_g1_count * kliep_mae_g1 +
       calibration_g2_count * kliep_mae_g2) / calibration_total
  )
  weight_pooled <- with(
    valid_weight_diagnostics,
    (calibration_g1_count * wcp_weight_mae_g1 +
       calibration_g2_count * wcp_weight_mae_g2) / calibration_total
  )
  pooled_MAE_verified <- isTRUE(all.equal(
    valid_weight_diagnostics$kliep_mae_all, KLIEP_pooled,
    tolerance = 1e-12
  )) && isTRUE(all.equal(
    valid_weight_diagnostics$wcp_weight_mae_all, weight_pooled,
    tolerance = 1e-12
  ))

  checks <- list(
    fixed_n1_n2_n_test = all(diagnostics$historical_total == 5000) &&
      all(diagnostics$current_total == 5000) &&
      sanity_result$configuration$n_test == 5000,
    balanced_realized_proportions_close =
      abs(balanced$historical_subgroup_1_proportion - 0.55) < 0.03 &&
      abs(balanced$current_subgroup_1_proportion - 0.45) < 0.03,
    small_realized_proportions_close =
      abs(small$historical_subgroup_2_proportion - 0.005) < 0.005 &&
      abs(small$current_subgroup_2_proportion - 0.01) < 0.007,
    small_expected_counts_close =
      abs(small$mean_historical_subgroup_2_count - 25) < 25 &&
      abs(small$mean_current_subgroup_2_count - 50) < 35,
    target_distribution_fixed = identical(target_g, c(0.5, 0.5)),
    conditional_distributions_shared = identical(
      sanity_result$configuration$only_setting_difference,
      c("q1", "q2")
    ),
    subgroup_mapping_correct = identical(
      original_subgroup_map, c(1L, 5L)
    ),
    H_equals_one = sanity_result$configuration$H == 1L,
    no_within_cluster_dependence =
      !sanity_result$configuration$within_cluster_dependence,
    all_weight_diagnostics_available =
      nrow(valid_weight_diagnostics) == M * length(smallgroup_settings),
    true_ratios_finite_positive =
      all(valid_weight_diagnostics$true_ratio_finite_positive),
    calibration_rows_aligned =
      all(valid_weight_diagnostics$calibration_row_alignment_verified),
    q1_hat_uses_same_training_sample = all(q1_hat_from_training),
    small_true_subgroup_ratio_g2_is_100 = all(
      weight_diagnostics$subgroup_ratio_true_g2[
        as.character(weight_diagnostics$setting) == "Small"
      ] == 100
    ),
    raw_weight_formulas_verified =
      all(valid_weight_diagnostics$raw_weight_formula_verified),
    raw_weights_not_normalized = all(
      valid_weight_diagnostics$raw_weight_scale ==
        "before proportions() normalization"
    ),
    MAE_values_nonnegative =
      length(metric_matrix) > 0L && all(is.finite(metric_matrix)) &&
        all(metric_matrix >= 0),
    overall_MAE_is_observation_level_pooled = pooled_MAE_verified,
    AR_CQR_reuse_same_rvec =
      all(valid_weight_diagnostics$same_rvec_reused_by_ar_cqr)
  )
  failed_checks <- names(checks)[!vapply(checks, isTRUE, logical(1))]
  if (length(failed_checks) > 0L) {
    stop("Sanity checks failed: ", paste(failed_checks, collapse = ", "))
  }

  list(
    checks = checks,
    mean_diagnostics = mean_diagnostics,
    weight_diagnostics_summary = summarize_weight_diagnostics(
      weight_diagnostics
    ),
    simulation = sanity_result
  )
}


# Nothing runs automatically when this file is sourced. Production example:
# smallgroup_1_result <- DA_WCP_simul_smallgroup_1(
#   n1 = 5000, n2 = 5000, n_test = 5000, M = 100\n# )


###############################################################################
# Stabilization study for the same two small-group settings
###############################################################################

# These are the actual objects read from the current _1 file.
smallgroup_2_settings <- smallgroup_settings
target_g_2 <- target_g
truncation_c_vec_2 <- c(2, 4, 6, 8, 10)
regularization_lambda_vec_2 <- seq(0, 0.9, by = 0.1)

stopifnot(
  identical(names(smallgroup_2_settings), c("Balanced", "Small")),
  identical(
    smallgroup_2_settings$Balanced$q1, c(0.55, 0.45)
  ),
  identical(
    smallgroup_2_settings$Balanced$q2, c(0.45, 0.55)
  ),
  identical(
    smallgroup_2_settings$Small$q1, c(0.995, 0.005)
  ),
  identical(
    smallgroup_2_settings$Small$q2, c(0.99, 0.01)
  ),
  identical(target_g_2, c(0.5, 0.5)),
  identical(original_subgroup_map, c(1L, 5L))
)


###############################################################################
# Base KLIEP estimate and score models (same formulas as _1 / DA_WCP)
###############################################################################

estimate_group_KLIEP_2 <- function(training, calibration, current,
                                   R, sig_wo_cv, kliep_cv) {
  validate_DA_inputs(
    training, calibration, current,
    data.frame(G = c(1L, 2L)),
    "Two-group KLIEP"
  )
  r_hat <- rep(NA_real_, nrow(calibration))
  for (g in 1:2) {
    training_g <- training[training$G == g, var_x, drop = FALSE]
    calibration_g <- calibration[calibration$G == g, var_x, drop = FALSE]
    current_g <- current[current$G == g, var_x, drop = FALSE]

    cat(
      g, "th group KLIEP progressing (n1 = ", nrow(training_g),
      ", n2 = ", nrow(current_g), ") \n", sep = ""
    )
    current_mean <- colMeans(current_g)
    current_sd <- apply(current_g, 2L, sd)
    if (any(!is.finite(current_sd)) || any(current_sd <= 0)) {
      stop("KLIEP normalization SD is invalid for G=", g, ".")
    }
    normalize <- function(data) {
      as.matrix(sweep(
        sweep(data, 2L, current_mean, "-"),
        2L, current_sd, "/"
      ))
    }
    estimate <- rvec_fun(
      normalize(calibration_g), normalize(training_g), normalize(current_g),
      R = R, sigma = sig_wo_cv, cv = kliep_cv
    )
    if (length(estimate) != nrow(calibration_g) ||
        any(!is.finite(estimate)) || any(estimate <= 0)) {
      stop("KLIEP returned invalid density ratios for G=", g, ".")
    }
    r_hat[calibration$G == g] <- estimate
  }
  if (any(!is.finite(r_hat)) || any(r_hat <= 0)) {
    stop("KLIEP ratios are not aligned with calibration rows.")
  }
  r_hat
}

estimate_base_weight_2 <- function(training, calibration, current,
                                   specification,
                                   R, sig_wo_cv, kliep_cv) {
  validate_DA_inputs(
    training, calibration, current,
    data.frame(G = c(1L, 2L)),
    "Two-group DA-WCP weight"
  )
  r_hat <- estimate_group_KLIEP_2(
    training, calibration, current,
    R = R, sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv
  )
  training_counts <- group_counts(training$G)
  q1_hat <- training_counts / nrow(training)
  if (any(!is.finite(q1_hat)) || any(q1_hat <= 0)) {
    stop("A training-sample q1_hat is zero or non-finite.")
  }
  subgroup_factor_hat <- target_g_2 / q1_hat
  r_true <- true_KLIEP_ratio(calibration)
  subgroup_factor_true <- target_g_2 / specification$q1
  estimated_raw_weight <- r_hat * subgroup_factor_hat[calibration$G]
  true_raw_weight <- r_true * subgroup_factor_true[calibration$G]
  if (any(!is.finite(estimated_raw_weight)) ||
      any(!is.finite(true_raw_weight)) ||
      any(estimated_raw_weight <= 0) || any(true_raw_weight <= 0)) {
    stop("Estimated or true raw WCP weights are invalid.")
  }
  list(
    r_hat = r_hat,
    r_true = r_true,
    q1_hat = q1_hat,
    subgroup_factor_hat = subgroup_factor_hat,
    subgroup_factor_true = subgroup_factor_true,
    calibration_G = calibration$G,
    estimated_raw_weight = estimated_raw_weight,
    true_raw_weight = true_raw_weight
  )
}

make_alpha_grid_2 <- function(n_alpha) {
  alpha_l <- seq(0, 0.5, length.out = n_alpha + 2L)[
    -c(1L, n_alpha + 2L)
  ]
  alpha_h <- seq(0.5, 1, length.out = n_alpha + 2L)[
    -c(1L, n_alpha + 2L)
  ]
  as.matrix(expand.grid(alpha_l, alpha_h))
}

fit_score_components_2 <- function(training, calibration, target,
                                   score, n_alpha) {
  formula_m <- as.formula(
    paste(var_y, "~", paste(c(var_x, var_g), collapse = " + "))
  )
  if (score == "AR") {
    model <- lm(formula_m, data = training)
    calibration_prediction <- as.vector(predict(model, calibration))
    target_prediction <- as.vector(predict(model, target))
    calibration_score <- abs(
      calibration_prediction - calibration[[var_y]]
    )
    if (any(!is.finite(calibration_score)) ||
        any(!is.finite(target_prediction))) {
      stop("AR score fitting produced non-finite values.")
    }
    return(list(
      score = "AR",
      calibration_score = calibration_score,
      target_prediction = target_prediction
    ))
  }

  if (score != "CQR") stop("score must be AR or CQR.")
  alpha_grid <- make_alpha_grid_2(n_alpha)
  grid_size <- nrow(alpha_grid)
  score_matrix <- matrix(NA_real_, nrow(calibration), grid_size)
  target_lower <- matrix(NA_real_, nrow(target), grid_size)
  target_upper <- matrix(NA_real_, nrow(target), grid_size)
  for (grid_index in seq_len(grid_size)) {
    lower_model <- rq(
      formula_m, tau = alpha_grid[grid_index, 1L], data = training
    )
    upper_model <- rq(
      formula_m, tau = alpha_grid[grid_index, 2L], data = training
    )
    calibration_lower <- as.vector(predict(lower_model, calibration))
    calibration_upper <- as.vector(predict(upper_model, calibration))
    target_lower[, grid_index] <- as.vector(predict(lower_model, target))
    target_upper[, grid_index] <- as.vector(predict(upper_model, target))
    score_matrix[, grid_index] <- pmax(
      calibration_lower - calibration[[var_y]],
      calibration[[var_y]] - calibration_upper
    )
  }
  if (any(!is.finite(score_matrix)) || any(!is.finite(target_lower)) ||
      any(!is.finite(target_upper))) {
    stop("CQR score fitting produced non-finite values.")
  }
  list(
    score = "CQR", alpha_grid = alpha_grid,
    calibration_score = score_matrix,
    target_lower = target_lower, target_upper = target_upper
  )
}

construct_interval_2 <- function(score_components, raw_weight,
                                 target, alpha) {
  n_calibration <- if (is.null(dim(score_components$calibration_score))) {
    length(score_components$calibration_score)
  } else {
    nrow(score_components$calibration_score)
  }
  if (length(raw_weight) != n_calibration ||
      any(!is.finite(raw_weight)) || any(raw_weight <= 0)) {
    stop("Raw weights are invalid or do not match calibration rows.")
  }
  normalized_weight <- proportions(raw_weight)

  if (score_components$score == "AR") {
    cutoff <- return_cutoff(
      normalized_weight, score_components$calibration_score, alpha
    )
    result_mat <- cbind(
      score_components$target_prediction - cutoff,
      score_components$target_prediction + cutoff
    )
  } else {
    grid_size <- ncol(score_components$calibration_score)
    cutoffs <- numeric(grid_size)
    average_lengths <- numeric(grid_size)
    for (grid_index in seq_len(grid_size)) {
      cutoffs[grid_index] <- return_cutoff(
        normalized_weight,
        score_components$calibration_score[, grid_index], alpha
      )
      average_lengths[grid_index] <- mean(
        2 * cutoffs[grid_index] +
          score_components$target_upper[, grid_index] -
          score_components$target_lower[, grid_index]
      )
    }
    selected <- which.min(average_lengths)
    result_mat <- cbind(
      score_components$target_lower[, selected] - cutoffs[selected],
      score_components$target_upper[, selected] + cutoffs[selected]
    )
  }
  if (any(!is.finite(result_mat))) {
    stop("Prediction interval endpoints are non-finite.")
  }
  covered <- target[[var_y]] >= result_mat[, 1L] &
    target[[var_y]] <= result_mat[, 2L]
  list(
    empirical_coverage = mean(covered),
    average_interval_length = mean(result_mat[, 2L] - result_mat[, 1L])
  )
}


###############################################################################
# Configurations and density-ratio-only transformations
###############################################################################

format_parameter_2 <- function(value) {
  format(value, trim = TRUE, scientific = FALSE, digits = 12)
}

make_configuration_table_2 <- function(c_choices, lambda_choices) {
  output <- bind_rows(
    data.frame(
      procedure = "Original", parameter = NA_real_,
      configuration = "Original", stringsAsFactors = FALSE
    ),
    data.frame(
      procedure = rep("Truncation", length(c_choices)),
      parameter = c_choices,
      configuration = paste0(
        "Truncation_c_",
        vapply(c_choices, format_parameter_2, character(1))
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      procedure = rep("Regularization", length(lambda_choices)),
      parameter = lambda_choices,
      configuration = paste0(
        "Regularization_lambda_",
        vapply(lambda_choices, format_parameter_2, character(1))
      ),
      stringsAsFactors = FALSE
    )
  )
  output$configuration_order <- seq_len(nrow(output))
  output
}

failed_safe_result_2 <- function(reason) {
  list(ok = FALSE, value = NULL, error = reason, warnings = "")
}

derive_weight_fit_2 <- function(base_fit, procedure, parameter = NA_real_) {
  if (!base_fit$ok) return(failed_safe_result_2(base_fit$error))
  safe_compute({
    base <- base_fit$value
    r_used <- switch(
      procedure,
      Original = base$r_hat,
      Truncation = pmin(pmax(base$r_hat, 1 / parameter), parameter),
      Regularization = (1 - parameter) * base$r_hat + parameter,
      stop("Unknown stabilization procedure.")
    )
    estimated_raw_weight <-
      r_used * base$subgroup_factor_hat[base$calibration_G]
    if (any(!is.finite(r_used)) || any(r_used <= 0) ||
        any(!is.finite(estimated_raw_weight)) ||
        any(estimated_raw_weight <= 0)) {
      stop("Transformed density ratios or raw weights are invalid.")
    }
    transformation_verified <- switch(
      procedure,
      Original = isTRUE(all.equal(r_used, base$r_hat, tolerance = 0)),
      Truncation = all(r_used >= 1 / parameter) &&
        all(r_used <= parameter) && isTRUE(all.equal(
          r_used, pmin(pmax(base$r_hat, 1 / parameter), parameter),
          tolerance = 0
        )),
      Regularization = isTRUE(all.equal(
        r_used, (1 - parameter) * base$r_hat + parameter,
        tolerance = 0
      ))
    )
    weight_formula_verified <- isTRUE(all.equal(
      estimated_raw_weight,
      r_used * base$subgroup_factor_hat[base$calibration_G],
      tolerance = 0
    ))
    list(
      estimated_raw_weight = estimated_raw_weight,
      true_raw_weight = base$true_raw_weight,
      r_hat = base$r_hat,
      r_used = r_used,
      q1_hat = base$q1_hat,
      subgroup_factor_hat = base$subgroup_factor_hat,
      subgroup_factor_true = base$subgroup_factor_true,
      calibration_G = base$calibration_G,
      transformation_verified = transformation_verified,
      weight_formula_verified = weight_formula_verified,
      subgroup_factor_unchanged = TRUE,
      transformation_scope = "KLIEP density ratio only",
      raw_weight_scale = "before proportions() normalization"
    )
  })
}


###############################################################################
# Per-repetition diagnostics and evaluation
###############################################################################

sample_count_row_2 <- function(repetition, setting, data1, training,
                               calibration, current) {
  count_one <- function(data, g) sum(data$G == g)
  data.frame(
    repetition = repetition, setting = setting,
    historical_g1_count = count_one(data1, 1L),
    historical_g2_count = count_one(data1, 2L),
    training_g1_count = count_one(training, 1L),
    training_g2_count = count_one(training, 2L),
    calibration_g1_count = count_one(calibration, 1L),
    calibration_g2_count = count_one(calibration, 2L),
    current_g1_count = count_one(current, 1L),
    current_g2_count = count_one(current, 2L),
    historical_g1_zero = count_one(data1, 1L) == 0L,
    historical_g2_zero = count_one(data1, 2L) == 0L,
    training_g1_zero = count_one(training, 1L) == 0L,
    training_g2_zero = count_one(training, 2L) == 0L,
    calibration_g1_zero = count_one(calibration, 1L) == 0L,
    calibration_g2_zero = count_one(calibration, 2L) == 0L,
    current_g1_zero = count_one(current, 1L) == 0L,
    current_g2_zero = count_one(current, 2L) == 0L,
    stringsAsFactors = FALSE
  )
}

empty_sample_count_row_2 <- function(repetition, setting) {
  output <- data.frame(repetition = repetition, setting = setting)
  for (sample in c("historical", "training", "calibration", "current")) {
    for (g in 1:2) {
      output[[paste0(sample, "_g", g, "_count")]] <- NA_integer_
      output[[paste0(sample, "_g", g, "_zero")]] <- NA
    }
  }
  output
}

make_weight_diagnostic_row_2 <- function(
    repetition, setting, configuration_row, weight_fit,
    calibration_G, sample_counts) {
  count_columns <- sample_counts[
    setdiff(names(sample_counts), c("repetition", "setting"))
  ]
  base_columns <- list(
    repetition = repetition,
    setting = setting,
    procedure = configuration_row$procedure,
    parameter = configuration_row$parameter,
    configuration = configuration_row$configuration,
    configuration_order = configuration_row$configuration_order
  )
  if (!weight_fit$ok) {
    return(data.frame(
      base_columns,
      weight_mae_g1 = NA_real_, weight_medae_g1 = NA_real_,
      weight_mae_g2 = NA_real_, weight_medae_g2 = NA_real_,
      weight_mae_all = NA_real_, weight_medae_all = NA_real_,
      computation_ok = FALSE,
      failure_reason = weight_fit$error,
      oracle_error_metric = "pointwise absolute error: MAE and MedAE",
      transformation_verified = FALSE,
      weight_formula_verified = FALSE,
      subgroup_factor_unchanged = FALSE,
      transformation_scope = "KLIEP density ratio only",
      base_r_hat_reused = TRUE,
      overall_mae_pools_all_rows = TRUE,
      overall_medae_pools_all_rows = TRUE,
      raw_weight_scale = "before proportions() normalization",
      q1_hat_g1 = NA_real_, q1_hat_g2 = NA_real_,
      subgroup_factor_hat_g1 = NA_real_, subgroup_factor_hat_g2 = NA_real_,
      count_columns,
      stringsAsFactors = FALSE
    ))
  }

  estimated <- weight_fit$value$estimated_raw_weight
  truth <- weight_fit$value$true_raw_weight
  if (length(estimated) != length(calibration_G) ||
      length(truth) != length(calibration_G)) {
    return(make_weight_diagnostic_row_2(
      repetition, setting, configuration_row,
      failed_safe_result_2("Raw weights do not align with calibration rows."),
      calibration_G, sample_counts
    ))
  }
  absolute_error <- abs(estimated - truth)
  group_metric <- function(g, fun) {
    index <- which(calibration_G == g)
    if (length(index) == 0L) NA_real_ else fun(absolute_error[index])
  }
  metrics <- c(
    weight_mae_g1 = group_metric(1L, mean),
    weight_medae_g1 = group_metric(1L, stats::median),
    weight_mae_g2 = group_metric(2L, mean),
    weight_medae_g2 = group_metric(2L, stats::median),
    weight_mae_all = mean(absolute_error),
    weight_medae_all = stats::median(absolute_error)
  )
  if (any(!is.finite(metrics))) {
    return(make_weight_diagnostic_row_2(
      repetition, setting, configuration_row,
      failed_safe_result_2("A raw-weight MAE or MedAE is non-finite."),
      calibration_G, sample_counts
    ))
  }
  calibration_counts <- group_counts(calibration_G)
  pooled_mae <- sum(
    calibration_counts * metrics[c("weight_mae_g1", "weight_mae_g2")]
  ) / sum(calibration_counts)
  if (!isTRUE(all.equal(
    unname(metrics["weight_mae_all"]), unname(pooled_mae),
    tolerance = 1e-12
  ))) {
    stop("Overall weight MAE is not pooled over calibration rows.")
  }

  data.frame(
    base_columns,
    weight_mae_g1 = unname(metrics["weight_mae_g1"]),
    weight_medae_g1 = unname(metrics["weight_medae_g1"]),
    weight_mae_g2 = unname(metrics["weight_mae_g2"]),
    weight_medae_g2 = unname(metrics["weight_medae_g2"]),
    weight_mae_all = unname(metrics["weight_mae_all"]),
    weight_medae_all = unname(metrics["weight_medae_all"]),
    computation_ok = TRUE,
    failure_reason = NA_character_,
    oracle_error_metric = "pointwise absolute error: MAE and MedAE",
    transformation_verified = weight_fit$value$transformation_verified,
    weight_formula_verified = weight_fit$value$weight_formula_verified,
    subgroup_factor_unchanged = weight_fit$value$subgroup_factor_unchanged,
    transformation_scope = weight_fit$value$transformation_scope,
    base_r_hat_reused = TRUE,
    overall_mae_pools_all_rows = TRUE,
    overall_medae_pools_all_rows = TRUE,
    raw_weight_scale = weight_fit$value$raw_weight_scale,
    q1_hat_g1 = weight_fit$value$q1_hat[1L],
    q1_hat_g2 = weight_fit$value$q1_hat[2L],
    subgroup_factor_hat_g1 = weight_fit$value$subgroup_factor_hat[1L],
    subgroup_factor_hat_g2 = weight_fit$value$subgroup_factor_hat[2L],
    count_columns,
    stringsAsFactors = FALSE
  )
}

make_performance_row_2 <- function(
    repetition, setting, configuration_row, score,
    weight_fit, score_fit, target, weight_diagnostic,
    sample_counts, alpha) {
  failure_reasons <- character(0)
  if (!weight_fit$ok) {
    failure_reasons <- c(
      failure_reasons, paste0("Weight estimation: ", weight_fit$error)
    )
  }
  if (!score_fit$ok) {
    failure_reasons <- c(
      failure_reasons, paste0("Score fitting: ", score_fit$error)
    )
  }
  interval_fit <- if (length(failure_reasons) == 0L) {
    safe_compute(construct_interval_2(
      score_fit$value, weight_fit$value$estimated_raw_weight,
      target, alpha
    ))
  } else {
    failed_safe_result_2(paste(failure_reasons, collapse = " | "))
  }
  interval_ok <- interval_fit$ok
  count_columns <- sample_counts[
    setdiff(names(sample_counts), c("repetition", "setting"))
  ]
  data.frame(
    repetition = repetition,
    setting = setting,
    procedure = configuration_row$procedure,
    parameter = configuration_row$parameter,
    configuration = configuration_row$configuration,
    configuration_order = configuration_row$configuration_order,
    score = score,
    empirical_coverage = if (interval_ok) {
      interval_fit$value$empirical_coverage
    } else NA_real_,
    average_interval_length = if (interval_ok) {
      interval_fit$value$average_interval_length
    } else NA_real_,
    computation_ok = interval_ok,
    failure_reason = if (interval_ok) NA_character_ else interval_fit$error,
    weight_mae_g1 = weight_diagnostic$weight_mae_g1,
    weight_medae_g1 = weight_diagnostic$weight_medae_g1,
    weight_mae_g2 = weight_diagnostic$weight_mae_g2,
    weight_medae_g2 = weight_diagnostic$weight_medae_g2,
    weight_mae_all = weight_diagnostic$weight_mae_all,
    weight_medae_all = weight_diagnostic$weight_medae_all,
    weight_computation_ok = weight_diagnostic$computation_ok,
    weight_failure_reason = weight_diagnostic$failure_reason,
    count_columns,
    stringsAsFactors = FALSE
  )
}

run_smallgroup_2_setting_replicate <- function(
    repetition, setting, specification, target,
    n1, n2, alpha, n_alpha, R, sig_wo_cv, kliep_cv,
    c_choices, lambda_choices) {
  cat("Setting: ", setting, "\n", sep = "")

  # Exactly the seed and data-generation convention used in _1.
  set.seed(123L + as.integer(repetition))
  data1 <- generate_smallgroup_data(
    n1, t = 1L, subgroup_probabilities = specification$q1
  )
  set.seed(100123L + as.integer(repetition))
  data2 <- generate_smallgroup_data(
    n2, t = 2L, subgroup_probabilities = specification$q2
  )
  splits <- make_smallgroup_splits(data1, repetition)

  data1_y <- data1[, c(var_x, var_y, var_g), drop = FALSE]
  current <- data2[, c(var_x, var_g), drop = FALSE]
  target_y <- target[, c(var_x, var_y, var_g), drop = FALSE]
  training <- data1_y[splits$marginal_train, , drop = FALSE]
  calibration <- data1_y[splits$marginal_cal, , drop = FALSE]
  sample_counts <- sample_count_row_2(
    repetition, setting, data1_y, training, calibration, current
  )

  base_fit <- safe_compute(estimate_base_weight_2(
    training, calibration, current, specification,
    R = R, sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv
  ))
  score_fits <- setNames(lapply(c("AR", "CQR"), function(score) {
    safe_compute(fit_score_components_2(
      training, calibration, target_y, score, n_alpha
    ))
  }), c("AR", "CQR"))

  configurations <- make_configuration_table_2(c_choices, lambda_choices)
  diagnostic_rows <- vector("list", nrow(configurations))
  performance_rows <- vector("list", nrow(configurations) * 2L)
  performance_index <- 0L
  for (configuration_index in seq_len(nrow(configurations))) {
    configuration_row <- configurations[configuration_index, , drop = FALSE]
    weight_fit <- derive_weight_fit_2(
      base_fit, configuration_row$procedure, configuration_row$parameter
    )
    diagnostic <- make_weight_diagnostic_row_2(
      repetition, setting, configuration_row, weight_fit,
      calibration$G, sample_counts
    )
    diagnostic_rows[[configuration_index]] <- diagnostic
    for (score in c("AR", "CQR")) {
      performance_index <- performance_index + 1L
      performance_rows[[performance_index]] <- make_performance_row_2(
        repetition, setting, configuration_row, score,
        weight_fit, score_fits[[score]], target_y,
        diagnostic, sample_counts, alpha
      )
    }
  }

  list(
    results = bind_rows(performance_rows),
    weight_diagnostics = bind_rows(diagnostic_rows),
    sample_counts = sample_counts,
    rvec_original = if (base_fit$ok) base_fit$value$r_hat else NULL
  )
}

failed_setting_result_2 <- function(repetition, setting,
                                    c_choices, lambda_choices,
                                    reason) {
  configurations <- make_configuration_table_2(c_choices, lambda_choices)
  sample_counts <- empty_sample_count_row_2(repetition, setting)
  diagnostic_rows <- lapply(seq_len(nrow(configurations)), function(index) {
    make_weight_diagnostic_row_2(
      repetition, setting, configurations[index, , drop = FALSE],
      failed_safe_result_2(reason), integer(0), sample_counts
    )
  })
  diagnostics <- bind_rows(diagnostic_rows)
  performance_rows <- lapply(seq_len(nrow(configurations)), function(index) {
    bind_rows(lapply(c("AR", "CQR"), function(score) {
      make_performance_row_2(
        repetition, setting, configurations[index, , drop = FALSE], score,
        failed_safe_result_2(reason), failed_safe_result_2(reason),
        data.frame(Y = numeric(0)), diagnostics[index, , drop = FALSE],
        sample_counts, alpha = 0.2
      )
    }))
  })
  list(
    results = bind_rows(performance_rows),
    weight_diagnostics = diagnostics,
    sample_counts = sample_counts,
    rvec_original = NULL
  )
}

run_smallgroup_2_replicate <- function(
    repetition, n1, n2, n_test,
    alpha, n_alpha, R, sig_wo_cv, kliep_cv,
    c_choices, lambda_choices) {
  cat("------------------------------------------------------\n")
  cat("  Two-group stabilization - replicate ", repetition, "\n", sep = "")
  cat("------------------------------------------------------\n")

  set.seed(200123L + as.integer(repetition))
  target <- generate_smallgroup_data(n_test, t = 2L, test = TRUE)
  setting_results <- lapply(names(smallgroup_2_settings), function(setting) {
    tryCatch(
      run_smallgroup_2_setting_replicate(
        repetition = repetition,
        setting = setting,
        specification = smallgroup_2_settings[[setting]],
        target = target,
        n1 = n1, n2 = n2,
        alpha = alpha, n_alpha = n_alpha, R = R,
        sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
        c_choices = c_choices, lambda_choices = lambda_choices
      ),
      error = function(e) failed_setting_result_2(
        repetition, setting, c_choices, lambda_choices,
        paste("Setting-level failure:", conditionMessage(e))
      )
    )
  })
  names(setting_results) <- names(smallgroup_2_settings)
  setting_results
}


###############################################################################
# Summaries and user-facing outputs
###############################################################################

median_or_na_2 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else stats::median(x)
}

summarize_sample_counts_2 <- function(sample_counts) {
  setting <- unique(as.character(sample_counts$setting))
  if (length(setting) != 1L) {
    stop("Sample-count summary requires exactly one setting.")
  }
  rows <- list()
  row_index <- 0L
  for (sample in c("Historical", "Training", "Calibration", "Current")) {
    prefix <- tolower(sample)
    for (g in 1:2) {
      row_index <- row_index + 1L
      values <- sample_counts[[paste0(prefix, "_g", g, "_count")]]
      zero_values <- sample_counts[[paste0(prefix, "_g", g, "_zero")]]
      rows[[row_index]] <- data.frame(
        setting = setting, sample = sample, group = paste0("G=", g),
        mean = mean_or_na(values), sd = sd_or_na(values),
        min = if (all(!is.finite(values))) NA_real_ else min(values, na.rm = TRUE),
        median = median_or_na_2(values),
        max = if (all(!is.finite(values))) NA_real_ else max(values, na.rm = TRUE),
        zero_frequency = if (all(is.na(zero_values))) {
          NA_real_
        } else mean(zero_values, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(rows)
}

summarize_smallgroup_2_results <- function(results, weight_diagnostics,
                                           sample_counts) {
  performance <- results %>%
    group_by(
      setting, procedure, parameter, configuration,
      configuration_order, score
    ) %>%
    summarise(
      mean_empirical_coverage = mean_or_na(empirical_coverage),
      sd_empirical_coverage = sd_or_na(empirical_coverage),
      mean_interval_length = mean_or_na(average_interval_length),
      sd_interval_length = sd_or_na(average_interval_length),
      n_valid_coverage = sum(
        computation_ok & is.finite(empirical_coverage) &
          is.finite(average_interval_length)
      ),
      n_failed = sum(!computation_ok),
      failure_rate = mean(!computation_ok),
      .groups = "drop"
    )

  weight <- weight_diagnostics %>%
    group_by(
      setting, procedure, parameter, configuration, configuration_order
    ) %>%
    summarise(
      mean_weight_mae_g1 = mean_or_na(weight_mae_g1),
      sd_weight_mae_g1 = sd_or_na(weight_mae_g1),
      median_weight_mae_g1 = median_or_na_2(weight_mae_g1),
      mean_weight_medae_g1 = mean_or_na(weight_medae_g1),
      sd_weight_medae_g1 = sd_or_na(weight_medae_g1),
      median_weight_medae_g1 = median_or_na_2(weight_medae_g1),
      mean_weight_mae_g2 = mean_or_na(weight_mae_g2),
      sd_weight_mae_g2 = sd_or_na(weight_mae_g2),
      median_weight_mae_g2 = median_or_na_2(weight_mae_g2),
      mean_weight_medae_g2 = mean_or_na(weight_medae_g2),
      sd_weight_medae_g2 = sd_or_na(weight_medae_g2),
      median_weight_medae_g2 = median_or_na_2(weight_medae_g2),
      mean_weight_mae_all = mean_or_na(weight_mae_all),
      sd_weight_mae_all = sd_or_na(weight_mae_all),
      median_weight_mae_all = median_or_na_2(weight_mae_all),
      mean_weight_medae_all = mean_or_na(weight_medae_all),
      sd_weight_medae_all = sd_or_na(weight_medae_all),
      median_weight_medae_all = median_or_na_2(weight_medae_all),
      n_valid_weight_mae_g1 = sum(is.finite(weight_mae_g1)),
      n_valid_weight_medae_g1 = sum(is.finite(weight_medae_g1)),
      n_valid_weight_mae_g2 = sum(is.finite(weight_mae_g2)),
      n_valid_weight_medae_g2 = sum(is.finite(weight_medae_g2)),
      n_valid_weight_mae_all = sum(is.finite(weight_mae_all)),
      n_valid_weight_medae_all = sum(is.finite(weight_medae_all)),
      n_failed_weight_estimation = sum(!computation_ok),
      weight_failure_rate = mean(!computation_ok),
      .groups = "drop"
    )
  combined <- left_join(
    performance, weight,
    by = c(
      "setting", "procedure", "parameter",
      "configuration", "configuration_order"
    )
  )
  list(
    performance = arrange(performance, configuration_order, score),
    weight = arrange(weight, configuration_order),
    sample_counts = summarize_sample_counts_2(sample_counts),
    combined = arrange(combined, configuration_order, score)
  )
}

configuration_labels_2 <- function(configuration_table) {
  ifelse(
    configuration_table$procedure == "Original",
    "Original",
    ifelse(
      configuration_table$procedure == "Truncation",
      paste0("c=", vapply(
        configuration_table$parameter, format_parameter_2, character(1)
      )),
      paste0("lambda=", vapply(
        configuration_table$parameter, format_parameter_2, character(1)
      ))
    )
  )
}

make_smallgroup_2_coverage_figures <- function(results, setting, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  plot_data <- results[results$computation_ok, , drop = FALSE]
  plot_data <- plot_data[is.finite(plot_data$empirical_coverage), , drop = FALSE]
  plot_data$score <- factor(plot_data$score, levels = c("AR", "CQR"))
  setting_slug <- switch(
    setting, Balanced = "balanced", Small = "small",
    stop("Unknown small-group setting: ", setting)
  )
  coverage_y_limits <- switch(
    setting, Balanced = c(0.6, 0.9), Small = c(0, 1)
  )
  coverage_y_breaks <- switch(
    setting, Balanced = seq(0.6, 0.9, by = 0.1),
    Small = seq(0, 1, by = 0.25)
  )
  common_coverage_theme <- theme_test(base_size = 18) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(2, 2, 2, 2, unit = "pt"),
      axis.text.x = element_text(margin = margin(t = 2, unit = "pt")),
      axis.title.x = element_text(margin = margin(t = 3, unit = "pt")),
      axis.title.y = element_text(margin = margin(r = 3, unit = "pt")),
      legend.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.box.spacing = grid::unit(1, "pt"),
      legend.spacing.x = grid::unit(4, "pt")
    )

  truncation_data <- plot_data[
    plot_data$procedure %in% c("Truncation", "Original"), , drop = FALSE
  ]
  truncation_data$plot_position <- ifelse(
    truncation_data$procedure == "Original",
    length(truncation_c_vec_2) + 2,
    match(truncation_data$parameter, truncation_c_vec_2)
  )
  truncation_labels <- c(
    paste0("c=", vapply(
      truncation_c_vec_2, format_parameter_2, character(1)
    )),
    "No truncation"
  )
  truncation_breaks <- c(
    seq_along(truncation_c_vec_2), length(truncation_c_vec_2) + 2
  )
  truncation_plot <- ggplot(
    truncation_data,
    aes(
      x = plot_position, y = empirical_coverage, fill = score,
      group = interaction(plot_position, score)
    )
  ) +
    geom_boxplot(
      outlier.size = 0.6, width = 0.58,
      position = position_dodge(width = 0.7)
    ) +
    geom_hline(
      yintercept = 0.8, color = "red",
      linetype = "dashed", linewidth = 0.7
    ) +
    geom_vline(
      xintercept = length(truncation_c_vec_2) + 1,
      color = "grey45", linetype = "dashed", linewidth = 0.6
    ) +
    scale_x_continuous(
      breaks = truncation_breaks, labels = truncation_labels,
      limits = c(0.4, length(truncation_c_vec_2) + 2.6)
    ) +
    scale_y_continuous(breaks = coverage_y_breaks) +
    coord_cartesian(ylim = coverage_y_limits) +
    scale_fill_manual(values = c(AR = "#9ecae1", CQR = "#3182bd")) +
    labs(
      x = "Truncation threshold", y = "Empirical coverage",
      fill = "Score"
    ) +
    common_coverage_theme

  regularization_data <- plot_data[
    plot_data$procedure == "Regularization", , drop = FALSE
  ]
  regularization_data$plot_position <- match(
    round(regularization_data$parameter, 10),
    round(regularization_lambda_vec_2, 10)
  )
  regularization_labels <- as.expression(lapply(
    regularization_lambda_vec_2,
    function(value) {
      if (value == 0) {
        bquote(atop(lambda == 0, "(No regularization)"))
      } else {
        bquote(.(sprintf("%.1f", value)))
      }
    }
  ))
  regularization_plot <- ggplot(
    regularization_data,
    aes(
      x = plot_position, y = empirical_coverage, fill = score,
      group = interaction(plot_position, score)
    )
  ) +
    geom_boxplot(
      outlier.size = 0.6, width = 0.58,
      position = position_dodge(width = 0.7)
    ) +
    geom_hline(
      yintercept = 0.8, color = "red",
      linetype = "dashed", linewidth = 0.7
    ) +
    scale_x_continuous(
      breaks = seq_along(regularization_lambda_vec_2),
      labels = regularization_labels,
      limits = c(0.4, length(regularization_lambda_vec_2) + 0.6)
    ) +
    scale_y_continuous(breaks = coverage_y_breaks) +
    coord_cartesian(ylim = coverage_y_limits) +
    scale_fill_manual(values = c(AR = "#9ecae1", CQR = "#3182bd")) +
    labs(
      x = expression(lambda), y = "Empirical coverage",
      fill = "Score"
    ) +
    common_coverage_theme

  truncation_path <- file.path(
    figures_dir,
    paste0(
      "simulation_smallgroup_2_", setting_slug,
      "_truncation_coverage.pdf"
    )
  )
  regularization_path <- file.path(
    figures_dir,
    paste0(
      "simulation_smallgroup_2_", setting_slug,
      "_regularization_coverage.pdf"
    )
  )
  save_pdf_plot(truncation_plot, truncation_path, width = 6, height = 5)
  save_pdf_plot(regularization_plot, regularization_path, width = 8, height = 5)
  c(
    truncation_coverage_pdf = truncation_path,
    regularization_coverage_pdf = regularization_path
  )
}

build_smallgroup_2_all_information <- function(simulation_result) {
  results <- simulation_result$results
  diagnostics <- simulation_result$weight_diagnostics
  summary <- simulation_result$summary
  diagnostic_keys <- c(
    "repetition", "setting", "procedure", "parameter",
    "configuration", "configuration_order"
  )
  diagnostic_extra <- setdiff(names(diagnostics), names(results))
  diagnostics_for_join <- diagnostics[, c(
    diagnostic_keys, diagnostic_extra
  ), drop = FALSE]
  results %>%
    left_join(diagnostics_for_join, by = diagnostic_keys) %>%
    left_join(
      summary,
      by = c(
        "setting", "procedure", "parameter", "configuration",
        "configuration_order", "score"
      )
    ) %>%
    arrange(repetition, configuration_order, score)
}

format_excel_parameter_2 <- function(procedure, parameter) {
  ifelse(
    procedure == "Original", "",
    vapply(parameter, format_parameter_2, character(1))
  )
}

make_smallgroup_2_excel_sheets <- function(simulation_result) {
  performance <- simulation_result$performance_summary %>%
    transmute(
      Procedure = procedure,
      Parameter = format_excel_parameter_2(procedure, parameter),
      Score = score,
      `Mean empirical coverage` = mean_empirical_coverage,
      `SD empirical coverage` = sd_empirical_coverage,
      `Mean interval length` = mean_interval_length,
      `SD interval length` = sd_interval_length,
      `Valid repetitions` = n_valid_coverage,
      `Failed repetitions` = n_failed,
      `Failure rate` = failure_rate
    )
  weight <- simulation_result$weight_summary %>%
    transmute(
      Procedure = procedure,
      Parameter = format_excel_parameter_2(procedure, parameter),
      `Mean G=1 weight MAE` = mean_weight_mae_g1,
      `SD G=1 weight MAE` = sd_weight_mae_g1,
      `Median G=1 weight MAE` = median_weight_mae_g1,
      `Mean G=1 weight MedAE` = mean_weight_medae_g1,
      `SD G=1 weight MedAE` = sd_weight_medae_g1,
      `Median G=1 weight MedAE` = median_weight_medae_g1,
      `Mean G=2 weight MAE` = mean_weight_mae_g2,
      `SD G=2 weight MAE` = sd_weight_mae_g2,
      `Median G=2 weight MAE` = median_weight_mae_g2,
      `Mean G=2 weight MedAE` = mean_weight_medae_g2,
      `SD G=2 weight MedAE` = sd_weight_medae_g2,
      `Median G=2 weight MedAE` = median_weight_medae_g2,
      `Mean overall weight MAE` = mean_weight_mae_all,
      `SD overall weight MAE` = sd_weight_mae_all,
      `Median overall weight MAE` = median_weight_mae_all,
      `Mean overall weight MedAE` = mean_weight_medae_all,
      `SD overall weight MedAE` = sd_weight_medae_all,
      `Median overall weight MedAE` = median_weight_medae_all,
      `Valid G=1 MAE` = n_valid_weight_mae_g1,
      `Valid G=1 MedAE` = n_valid_weight_medae_g1,
      `Valid G=2 MAE` = n_valid_weight_mae_g2,
      `Valid G=2 MedAE` = n_valid_weight_medae_g2,
      `Valid overall MAE` = n_valid_weight_mae_all,
      `Valid overall MedAE` = n_valid_weight_medae_all,
      `Failed weight estimates` = n_failed_weight_estimation
    )
  counts <- simulation_result$sample_count_summary %>%
    transmute(
      Sample = sample, Group = group,
      Mean = mean, SD = sd, Min = min,
      Median = median, Max = max,
      `Zero frequency` = zero_frequency
    )
  configuration <- simulation_result$configuration
  config_table <- data.frame(
    Field = c(
      "Setting", "q1", "q2", "target_g", "n1", "n2", "n_test", "M",
      "alpha", "n_alpha", "c_choices", "lambda_choices",
      "R", "sig_wo_cv", "kliep_cv", "seed convention",
      "train/calibration split"
    ),
    Value = c(
      configuration$setting,
      paste(configuration$q1, collapse = ", "),
      paste(configuration$q2, collapse = ", "),
      paste(configuration$target_g, collapse = ", "),
      configuration$n1, configuration$n2, configuration$n_test,
      configuration$M, configuration$alpha, configuration$n_alpha,
      paste(configuration$c_choices, collapse = ", "),
      paste(configuration$lambda_choices, collapse = ", "),
      configuration$R_kliep, configuration$sig_wo_cv,
      as.character(configuration$kliep_cv),
      configuration$seed_convention,
      configuration$train_calibration_split
    ),
    stringsAsFactors = FALSE
  )
  list(
    Performance = as.data.frame(performance),
    WeightDiagnostics = as.data.frame(weight),
    SampleCounts = as.data.frame(counts),
    Configuration = config_table
  )
}

xml_escape_2 <- function(value) {
  value <- gsub("&", "&amp;", as.character(value), fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  value <- gsub('"', "&quot;", value, fixed = TRUE)
  gsub("'", "&apos;", value, fixed = TRUE)
}

excel_column_2 <- function(index) {
  output <- character(0)
  while (index > 0L) {
    remainder <- (index - 1L) %% 26L
    output <- c(LETTERS[remainder + 1L], output)
    index <- (index - remainder - 1L) %/% 26L
  }
  paste(output, collapse = "")
}

sheet_xml_2 <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  character_data <- as.data.frame(
    lapply(data, as.character), stringsAsFactors = FALSE,
    check.names = FALSE
  )
  header <- as.data.frame(
    as.list(names(data)), stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(header) <- names(data)
  all_rows <- rbind(header, character_data)
  rows <- vector("character", nrow(all_rows))
  for (row_index in seq_len(nrow(all_rows))) {
    cells <- character(ncol(all_rows))
    for (column_index in seq_len(ncol(all_rows))) {
      reference <- paste0(excel_column_2(column_index), row_index)
      value <- all_rows[row_index, column_index]
      style <- if (row_index == 1L) ' s="1"' else ""
      if (is.na(value)) value <- ""
      cells[column_index] <- paste0(
        '<c r="', reference, '" t="inlineStr"', style,
        '><is><t xml:space="preserve">', xml_escape_2(value),
        "</t></is></c>"
      )
    }
    rows[row_index] <- paste0(
      '<row r="', row_index, '">', paste(cells, collapse = ""), "</row>"
    )
  }
  paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
    '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>',
    '<sheetData>', paste(rows, collapse = ""), "</sheetData></worksheet>"
  )
}

write_xlsx_base_2 <- function(sheets, path) {
  staging <- tempfile("smallgroup_2_xlsx_")
  dir.create(file.path(staging, "_rels"), recursive = TRUE)
  dir.create(file.path(staging, "xl", "_rels"), recursive = TRUE)
  dir.create(file.path(staging, "xl", "worksheets"), recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)

  write_xml <- function(relative_path, text) {
    writeLines(text, file.path(staging, relative_path), useBytes = TRUE)
  }
  sheet_overrides <- paste0(
    '<Override PartName="/xl/worksheets/sheet', seq_along(sheets),
    '.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
  )
  write_xml("[Content_Types].xml", paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
    '<Default Extension="xml" ContentType="application/xml"/>',
    '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
    '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
    paste(sheet_overrides, collapse = ""), "</Types>"
  ))
  write_xml("_rels/.rels", paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>',
    "</Relationships>"
  ))
  sheet_nodes <- paste0(
    '<sheet name="', xml_escape_2(names(sheets)), '" sheetId="',
    seq_along(sheets), '" r:id="rId', seq_along(sheets), '"/>'
  )
  write_xml("xl/workbook.xml", paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
    '<sheets>', paste(sheet_nodes, collapse = ""), "</sheets></workbook>"
  ))
  relationships <- paste0(
    '<Relationship Id="rId', seq_along(sheets),
    '" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet',
    seq_along(sheets), '.xml"/>'
  )
  write_xml("xl/_rels/workbook.xml.rels", paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    paste(relationships, collapse = ""),
    '<Relationship Id="rId', length(sheets) + 1L,
    '" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>',
    "</Relationships>"
  ))
  write_xml("xl/styles.xml", paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
    '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><color rgb="FFFFFFFF"/><sz val="11"/><name val="Calibri"/></font></fonts>',
    '<fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF1F4E78"/><bgColor indexed="64"/></patternFill></fill></fills>',
    '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>',
    '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>',
    '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/></cellXfs>',
    "</styleSheet>"
  ))
  for (index in seq_along(sheets)) {
    write_xml(
      file.path("xl", "worksheets", paste0("sheet", index, ".xml")),
      sheet_xml_2(sheets[[index]])
    )
  }

  path <- normalizePath(path, mustWork = FALSE)
  if (file.exists(path)) unlink(path)
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(path, list.files(staging, recursive = TRUE, full.names = TRUE),
              root = staging)
  } else if (nzchar(Sys.which("zip"))) {
    previous <- getwd()
    on.exit(setwd(previous), add = TRUE)
    setwd(staging)
    status <- system2("zip", c("-q", "-r", shQuote(path), "."))
    setwd(previous)
    if (!identical(status, 0L)) stop("The zip command failed for Excel output.")
  } else if (.Platform$OS.type == "windows") {
    powershell <- Sys.which("powershell.exe")
    if (!nzchar(powershell)) powershell <- Sys.which("pwsh.exe")
    if (!nzchar(powershell)) stop("No ZIP backend is available for Excel output.")
    quote_ps <- function(x) gsub("'", "''", x, fixed = TRUE)
    command <- paste0(
      "Add-Type -AssemblyName System.IO.Compression.FileSystem; ",
      "[System.IO.Compression.ZipFile]::CreateFromDirectory('",
      quote_ps(staging), "','", quote_ps(path),
      "',[System.IO.Compression.CompressionLevel]::Optimal,$false)"
    )
    output <- system2(
      powershell,
      c("-NoProfile", "-NonInteractive", "-WindowStyle", "Hidden",
        "-Command", shQuote(command)),
      stdout = TRUE, stderr = TRUE
    )
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) {
      stop("PowerShell Excel archive creation failed: ",
           paste(output, collapse = " "))
    }
  } else {
    stop("No ZIP backend is available for Excel output.")
  }
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Excel output was not created.")
  }
  path
}

write_smallgroup_2_excel <- function(simulation_result, path) {
  sheets <- make_smallgroup_2_excel_sheets(simulation_result)
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    workbook <- openxlsx::createWorkbook()
    header_style <- openxlsx::createStyle(
      fgFill = "#1F4E78", fontColour = "#FFFFFF",
      textDecoration = "bold"
    )
    for (sheet_name in names(sheets)) {
      openxlsx::addWorksheet(workbook, sheet_name, gridLines = FALSE)
      openxlsx::writeData(
        workbook, sheet_name, sheets[[sheet_name]],
        headerStyle = header_style, withFilter = TRUE
      )
      openxlsx::freezePane(workbook, sheet_name, firstRow = TRUE)
      openxlsx::setColWidths(
        workbook, sheet_name, cols = seq_len(ncol(sheets[[sheet_name]])),
        widths = "auto"
      )
    }
    openxlsx::saveWorkbook(workbook, path, overwrite = TRUE)
  } else if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheets, path)
  } else {
    write_xlsx_base_2(sheets, path)
  }
  path
}

write_smallgroup_2_outputs <- function(simulation_result,
                                       output_dir, figures_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  setting <- simulation_result$configuration$setting
  slug <- tolower(setting)
  csv_path <- file.path(
    output_dir, paste0("simulation_smallgroup_2_", slug, "_results.csv")
  )
  xlsx_path <- file.path(
    output_dir, paste0("simulation_smallgroup_2_", slug, "_summary.xlsx")
  )
  all_information <- build_smallgroup_2_all_information(simulation_result)
  write.csv(all_information, csv_path, row.names = FALSE)
  write_smallgroup_2_excel(simulation_result, xlsx_path)
  coverage_paths <- make_smallgroup_2_coverage_figures(
    simulation_result$results, setting, figures_dir
  )
  output_paths <- c(
    results_csv = csv_path, summary_xlsx = xlsx_path, coverage_paths
  )
  stopifnot(length(output_paths) == 4L, all(file.exists(output_paths)))
  output_paths
}


###############################################################################
# Sequential simulation driver
###############################################################################

validate_smallgroup_2_choices <- function(c_choices, lambda_choices) {
  if (!is.numeric(c_choices) || length(c_choices) == 0L ||
      any(!is.finite(c_choices)) || any(c_choices <= 1)) {
    stop("Every c_choices value must be finite and greater than 1.")
  }
  if (!is.numeric(lambda_choices) || length(lambda_choices) == 0L ||
      any(!is.finite(lambda_choices)) ||
      any(lambda_choices < 0 | lambda_choices >= 1)) {
    stop("Every lambda_choices value must be in [0, 1).")
  }
  list(
    c_choices = unique(as.numeric(c_choices)),
    lambda_choices = unique(as.numeric(lambda_choices))
  )
}

# Runs all stabilization choices sequentially for each small-group setting.
DA_WCP_simul_smallgroup_2 <- function(
    n1 = 5000, n2 = 5000, n_test = 5000,
    M = 100, alpha = 0.2, n_alpha = 5, R = 4,
    sig_wo_cv = 1.2, kliep_cv = FALSE,
    c_choices = truncation_c_vec_2,
    lambda_choices = regularization_lambda_vec_2,
    save_outputs = TRUE,
    output_dir = "../Figures",
    figures_dir = "../Figures") {
  n1 <- validate_fixed_n(n1, "n1")
  n2 <- validate_fixed_n(n2, "n2")
  n_test <- validate_fixed_n(n_test, "n_test")
  M <- validate_fixed_n(M, "M")
  if (!is.finite(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1).")
  }
  choices <- validate_smallgroup_2_choices(c_choices, lambda_choices)
  c_choices <- choices$c_choices
  lambda_choices <- choices$lambda_choices
  configurations <- make_configuration_table_2(c_choices, lambda_choices)

  tasks <- lapply(seq_len(M), function(repetition) {
    list(repetition = repetition)
  })
  common_args <- list(
    n1 = n1, n2 = n2, n_test = n_test,
    alpha = alpha, n_alpha = n_alpha, R = R,
    sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
    c_choices = c_choices, lambda_choices = lambda_choices
  )

    replicate_results <- lapply(tasks, function(task) {
      do.call(run_smallgroup_2_replicate, c(task, common_args))
    })
    execution_info <- list(
      backend = "sequential", cores = 1L,
      details = list(), threads = 1L
    )

  setting_results <- lapply(names(smallgroup_2_settings), function(setting) {
    repetitions <- lapply(replicate_results, `[[`, setting)
    results <- bind_rows(lapply(repetitions, `[[`, "results"))
    weight_diagnostics <- bind_rows(lapply(
      repetitions, `[[`, "weight_diagnostics"
    ))
    sample_counts <- bind_rows(lapply(repetitions, `[[`, "sample_counts"))
    rvec_original <- lapply(repetitions, `[[`, "rvec_original")
    stopifnot(
      nrow(results) == M * nrow(configurations) * 2L,
      nrow(weight_diagnostics) == M * nrow(configurations),
      nrow(sample_counts) == M
    )
    summaries <- summarize_smallgroup_2_results(
      results, weight_diagnostics, sample_counts
    )
    specification <- smallgroup_2_settings[[setting]]
    configuration <- list(
      scenario = "S1 two-group small-subgroup stabilization",
      setting = setting,
      q1 = specification$q1,
      q2 = specification$q2,
      target_g = target_g_2,
      original_subgroup_map = original_subgroup_map,
      n1 = n1, n2 = n2, n_test = n_test,
      fixed_sample_sizes = TRUE, poisson_sample_size = FALSE,
      M = M, H = 1L, alpha = alpha, n_alpha = n_alpha,
      R_kliep = R, sig_wo_cv = sig_wo_cv, kliep_cv = kliep_cv,
      c_choices = c_choices, lambda_choices = lambda_choices,
      configurations = configurations,
      seed_convention = paste(
        "target: 200123+rep; historical: 123+rep;",
        "current: 100123+rep; split: 100000+rep"
      ),
      train_calibration_split = "individual-level random 50:50",
      expected_historical_counts = n1 * specification$q1,
      expected_training_counts = n1 * specification$q1 / 2,
      expected_calibration_counts = n1 * specification$q1 / 2,
      true_subgroup_factors = target_g_2 / specification$q1,
      procedures = c("Original", "Truncation", "Regularization"),
      original_r_hat_reused = TRUE,
      stabilization_scope = "KLIEP density ratio only",
      subgroup_factor_modified = FALSE,
      expected_output_files = c(
        results_csv = file.path(
          output_dir,
          paste0(
            "simulation_smallgroup_2_", tolower(setting), "_results.csv"
          )
        ),
        summary_xlsx = file.path(
          output_dir,
          paste0(
            "simulation_smallgroup_2_", tolower(setting), "_summary.xlsx"
          )
        ),
        truncation_coverage_pdf = file.path(
          figures_dir,
          paste0(
            "simulation_smallgroup_2_", tolower(setting),
            "_truncation_coverage.pdf"
          )
        ),
        regularization_coverage_pdf = file.path(
          figures_dir,
          paste0(
            "simulation_smallgroup_2_", tolower(setting),
            "_regularization_coverage.pdf"
          )
        )
      ),
      execution_backend = execution_info$backend,
      execution_cores = execution_info$cores,
      threads = 1L
    )
    setting_result <- list(
      results = results,
      weight_diagnostics = weight_diagnostics,
      sample_counts = sample_counts,
      sample_count_summary = summaries$sample_counts,
      performance_summary = summaries$performance,
      weight_summary = summaries$weight,
      summary = summaries$combined,
      rvec_original = rvec_original,
      execution_info = execution_info,
      configuration = configuration
    )
    if (isTRUE(save_outputs)) {
      setting_result$output_paths <- write_smallgroup_2_outputs(
        setting_result, output_dir, figures_dir
      )
    }
    setting_result
  })
  names(setting_results) <- names(smallgroup_2_settings)
  setting_results
}


###############################################################################
# Full-n M=2 sanity check
###############################################################################

run_smallgroup_2_sanity_check <- function(
    M = 2, save_outputs = FALSE,
    c_choices = truncation_c_vec_2,
    lambda_choices = regularization_lambda_vec_2) {
  sanity <- DA_WCP_simul_smallgroup_2(
    n1 = 5000, n2 = 5000, n_test = 5000, M = M,
    c_choices = c_choices, lambda_choices = lambda_choices,
    save_outputs = save_outputs
  )
  balanced <- sanity$Balanced
  small <- sanity$Small
  diagnostics <- bind_rows(
    balanced$weight_diagnostics, small$weight_diagnostics
  )
  successful <- diagnostics[diagnostics$computation_ok, , drop = FALSE]
  allowed_procedures <- c("Original", "Truncation", "Regularization")

  # Re-run _1 with the already estimated ratios, so this comparison adds no
  # KLIEP refitting and checks the Original implementation on identical data.
  cached_rvec <- lapply(seq_len(M), function(repetition) {
    list(
      Balanced = balanced$rvec_original[[repetition]],
      Small = small$rvec_original[[repetition]]
    )
  })
  reference_1 <- DA_WCP_simul_smallgroup_1(
    n1 = 5000, n2 = 5000, n_test = 5000, M = M,
    rvec = TRUE, rvec_marg = cached_rvec,
    save_outputs = FALSE
  )
  original_2 <- bind_rows(
    balanced$results[balanced$results$procedure == "Original", ],
    small$results[small$results$procedure == "Original", ]
  ) %>%
    arrange(repetition, setting, score)
  original_1 <- reference_1$results %>%
    arrange(repetition, setting, score)
  original_matches_1 <- isTRUE(all.equal(
    original_2$empirical_coverage,
    original_1$empirical_coverage,
    tolerance = 1e-12, check.attributes = FALSE
  )) && isTRUE(all.equal(
    original_2$average_interval_length,
    original_1$average_prediction_length,
    tolerance = 1e-12, check.attributes = FALSE
  ))
  regularization_zero <- bind_rows(
    balanced$results[
      balanced$results$procedure == "Regularization" &
        balanced$results$parameter == 0, ],
    small$results[
      small$results$procedure == "Regularization" &
        small$results$parameter == 0, ]
  ) %>%
    arrange(repetition, setting, score)
  regularization_zero_matches_original <-
    nrow(regularization_zero) == nrow(original_2) &&
    isTRUE(all.equal(
      regularization_zero$empirical_coverage,
      original_2$empirical_coverage,
      tolerance = 0, check.attributes = FALSE
    )) && isTRUE(all.equal(
      regularization_zero$average_interval_length,
      original_2$average_interval_length,
      tolerance = 0, check.attributes = FALSE
    )) && isTRUE(all.equal(
      regularization_zero$weight_mae_all,
      original_2$weight_mae_all,
      tolerance = 0, check.attributes = FALSE
    ))

  artificial_base <- list(
    ok = TRUE,
    value = list(
      r_hat = c(0.1, 3), r_true = c(0.1, 3),
      q1_hat = c(0.5, 0.005),
      subgroup_factor_hat = c(1, 100),
      subgroup_factor_true = c(1, 100),
      calibration_G = c(1L, 2L),
      estimated_raw_weight = c(0.1, 300),
      true_raw_weight = c(0.1, 300)
    )
  )
  artificial_truncation <- derive_weight_fit_2(
    artificial_base, "Truncation", 2
  )

  successful_counts <- with(
    successful, calibration_g1_count + calibration_g2_count
  )
  pooled_mae <- with(
    successful,
    (calibration_g1_count * weight_mae_g1 +
       calibration_g2_count * weight_mae_g2) / successful_counts
  )
  expected_files <- c(
    "simulation_smallgroup_2_balanced_results.csv",
    "simulation_smallgroup_2_balanced_summary.xlsx",
    "simulation_smallgroup_2_balanced_truncation_coverage.pdf",
    "simulation_smallgroup_2_balanced_regularization_coverage.pdf",
    "simulation_smallgroup_2_small_results.csv",
    "simulation_smallgroup_2_small_summary.xlsx",
    "simulation_smallgroup_2_small_truncation_coverage.pdf",
    "simulation_smallgroup_2_small_regularization_coverage.pdf"
  )
  checks <- list(
    exactly_two_settings = identical(names(sanity), c("Balanced", "Small")),
    exactly_two_subgroups = identical(target_g_2, c(0.5, 0.5)) &&
      length(original_subgroup_map) == 2L,
    balanced_probabilities_exact = identical(
      balanced$configuration$q1, c(0.55, 0.45)
    ) && identical(
      balanced$configuration$q2, c(0.45, 0.55)
    ),
    small_probabilities_exact = identical(
      small$configuration$q1, c(0.995, 0.005)
    ) && identical(
      small$configuration$q2, c(0.99, 0.01)
    ),
    target_probabilities_exact = identical(
      balanced$configuration$target_g, c(0.5, 0.5)
    ) && identical(small$configuration$target_g, c(0.5, 0.5)),
    fixed_sample_sizes = balanced$configuration$n1 == 5000 &&
      balanced$configuration$n2 == 5000 &&
      balanced$configuration$n_test == 5000 &&
      small$configuration$n1 == 5000 &&
      small$configuration$n2 == 5000 &&
      small$configuration$n_test == 5000,
    dgp_and_split_reused_from_1 = identical(
      balanced$configuration$seed_convention,
      small$configuration$seed_convention
    ) && identical(
      balanced$configuration$train_calibration_split,
      "individual-level random 50:50"
    ),
    original_matches_1 = original_matches_1,
    regularization_zero_matches_original =
      regularization_zero_matches_original,
    only_allowed_procedures = setequal(
      unique(as.character(balanced$results$procedure)), allowed_procedures
    ) && setequal(
      unique(as.character(small$results$procedure)), allowed_procedures
    ),
    same_stabilization_grid = identical(
      balanced$configuration$configurations,
      small$configuration$configurations
    ),
    truncation_changes_only_r = nrow(successful) > 0L && all(
      successful$transformation_scope == "KLIEP density ratio only"
    ) && artificial_truncation$ok &&
      max(artificial_truncation$value$r_used) <= 2 &&
      max(artificial_truncation$value$estimated_raw_weight) > 2,
    transformations_verified = nrow(successful) > 0L &&
      all(successful$transformation_verified) &&
      all(successful$weight_formula_verified),
    subgroup_factor_unchanged = nrow(successful) > 0L &&
      all(successful$subgroup_factor_unchanged),
    base_r_hat_reused = nrow(successful) > 0L &&
      all(successful$base_r_hat_reused),
    expected_historical_counts = identical(
      balanced$configuration$expected_historical_counts,
      c(2750, 2250)
    ) && identical(
      small$configuration$expected_historical_counts,
      c(4975, 25)
    ),
    expected_small_split_counts = identical(
      small$configuration$expected_training_counts, c(2487.5, 12.5)
    ) && identical(
      small$configuration$expected_calibration_counts, c(2487.5, 12.5)
    ),
    true_subgroup_factors = isTRUE(all.equal(
      balanced$configuration$true_subgroup_factors,
      c(0.5 / 0.55, 0.5 / 0.45), tolerance = 0
    )) && isTRUE(all.equal(
      small$configuration$true_subgroup_factors,
      c(0.5 / 0.995, 100), tolerance = 0
    )),
    mae_medae_nonnegative = nrow(successful) > 0L && all(
      is.finite(as.matrix(successful[, c(
        "weight_mae_g1", "weight_medae_g1",
        "weight_mae_g2", "weight_medae_g2",
        "weight_mae_all", "weight_medae_all"
      )]))
    ) && all(as.matrix(successful[, c(
      "weight_mae_g1", "weight_medae_g1",
      "weight_mae_g2", "weight_medae_g2",
      "weight_mae_all", "weight_medae_all"
    )]) >= 0),
    overall_mae_pools_all_rows = isTRUE(all.equal(
      successful$weight_mae_all, pooled_mae, tolerance = 1e-12
    )),
    overall_medae_is_direct = nrow(successful) > 0L &&
      all(successful$overall_medae_pools_all_rows),
    AR_CQR_share_weight_diagnostics = all(
      balanced$results$weight_mae_all[
        balanced$results$score == "AR"
      ] == balanced$results$weight_mae_all[
        balanced$results$score == "CQR"
      ], na.rm = TRUE
    ) && all(
      small$results$weight_mae_all[
        small$results$score == "AR"
      ] == small$results$weight_mae_all[
        small$results$score == "CQR"
      ], na.rm = TRUE
    ),
    no_KLIEP_error_diagnostics = !any(grepl(
      "kliep.*(mae|medae)", names(diagnostics), ignore.case = TRUE
    )),
    expected_output_names = identical(
      unname(c(
        basename(balanced$configuration$expected_output_files),
        basename(small$configuration$expected_output_files)
      )),
      expected_files
    ),
    fixed_result_dimensions = nrow(balanced$results) ==
      M * (1L + length(unique(c_choices)) +
             length(unique(lambda_choices))) * 2L &&
      nrow(small$results) ==
      M * (1L + length(unique(c_choices)) +
             length(unique(lambda_choices))) * 2L
  )
  failed_checks <- names(checks)[!vapply(checks, isTRUE, logical(1))]
  if (length(failed_checks) > 0L) {
    stop("Sanity checks failed: ", paste(failed_checks, collapse = ", "))
  }
  list(checks = checks, simulation = sanity, reference_1 = reference_1)
}


# Production example (not executed automatically):
# smallgroup_2_result <- DA_WCP_simul_smallgroup_2(
#   n1 = 5000, n2 = 5000, n_test = 5000, M = 100,
#   c_choices = c(2, 4, 6, 8, 10),
#   lambda_choices = seq(0, 0.9, by = 0.1)
# )

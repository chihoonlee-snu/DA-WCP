library(mvtnorm)
library(dplyr)
library(quantreg)
library(ggplot2)
library(tidyr)
library(nnet)
library(xgboost)
library(ranger)

source("KLIEP.R")
source("DA_WCP.R")


#####################################################################################
#
#                    Generating data & specifying distributions 
#
#####################################################################################

# ------------- Marginal distribution of Z ------------

prob1 <- matrix(nrow = 5, ncol = 5) # t = 1 data for 5 scenarios
prob2 <- matrix(nrow = 5, ncol = 5) # t = 2 data for 5 scenarios
target_g <- c(0.05, 0.1, 0.15, 0.3, 0.4) # t = 2 population

for(i in 1:5){
  if(i != 2) prob1[, i] <- rep(0.2, 5)
  else prob1[, i] <- target_g
}
for(i in 1:5){
  if(i != 3) prob2[, i] <- c(0.1, 0.15, 0.2, 0.25, 0.3)
  else prob2[, i] <- target_g
}

# --------- Conditional distribution X_t|Z_t ----------

# Base mean vectors 
mu1 <- function(k){
  return(c(k-3, sin(pi*k/2), 3-k, cos(pi*k/2)))
}

# Base covariance matrix
cov1 <- function(k, rho = 0.2){
  base_cov = matrix(ncol = 4, nrow = 4)
  for(i in 1:4){
    for(j in 1:4) base_cov[i, j] = rho^(abs(i-j))
  }
  return(base_cov + diag(rep(0.1*(k-1), 4)))
}

# Common shift vector, matrices
Delta_mu_shift <- c(0.3, -0.2, 0.4, 0.1)
Delta_Sigma_shift <- diag(c(0.05, -0.02, 0.03, 0.01))
Delta_mu_noshift <- c(0, 0, 0, 0)
Delta_Sigma_noshift <- diag(c(0, 0, 0, 0))


# ------- Outcome distribution Y_t|(X_t, Z_t) --------

beta <- c(5, -3, 2, 1)
gamma <- c(-10, -5, 0, 5, 10)

theta_hetero <- c(0.5, -0.4, 0.3, 0.2)
eta_hetero <- c(1.2, 1.6, 2.0, 2.4, 2.8)
theta_homo <- rep(0, 4)
eta_homo <- rep(log(5), 5)


# ----- Assigning parameters for each scenario -------

# parameter lists for scenarios (S1)-(S5)
theta_list = list(theta_hetero, theta_hetero, theta_hetero, theta_hetero, theta_homo)
eta_list = list(eta_hetero, eta_hetero, eta_hetero, eta_hetero, eta_homo)
Delta_mu_list = list(Delta_mu_shift, Delta_mu_shift, Delta_mu_shift, Delta_mu_noshift, Delta_mu_shift)
Delta_Sigma_list = list(Delta_Sigma_shift, Delta_Sigma_shift, Delta_Sigma_shift, Delta_Sigma_noshift, Delta_Sigma_shift)



# --------- Lists to save simulation results ---------

df_m = vector("list", 5) # list for 4 marginal methods
df_g = vector("list", 5) # list for 4 group-conditional methods
rvec_marg = vector("list", 5)
rvec_group = vector("list", 5)

# ------------ Data generation function -------------

# m : expected sample size
# t : time point (1 or 2)
# num : scenario number (1 to 5)
# test : whether the test dataset or not

generate_data <- function(m, t, num, test = F){
  if(test){
    n <- m
    Z <- sample(1:5, size = n, replace = TRUE, prob = target_g)
  }
  if(!test){
    n <- rpois(1, m)
    if(t == 1) Z <- sample(1:5, size = n, replace = TRUE, prob = prob1[, num])
    if(t == 2) Z <- sample(1:5, size = n, replace = TRUE, prob = prob2[, num])
  }
  X <- matrix(0, nrow = n, ncol = 4)
  Y <- numeric(n)
  Delta_mu = Delta_mu_list[[num]]
  Delta_Sigma = Delta_Sigma_list[[num]]
  theta = theta_list[[num]]
  eta = eta_list[[num]]
  
  for (k in 1:5) {
    idx <- which(Z == k)
    mu_k <- mu1(k) + if (t == 2) Delta_mu else 0
    Sigma_k <- cov1(k) + if (t == 2) Delta_Sigma else 0
    X[idx, ] <- MASS::mvrnorm(n = length(idx), mu = mu_k, Sigma = Sigma_k)
    for(i in idx) Y[i] <- X[i, ] %*% beta + gamma[k] + rnorm(1, sd = sqrt(exp(X[i, ] %*% theta + eta[k])))
  }
  
  df <- as.data.frame(X)
  colnames(df) <- paste0("X", 1:4)
  df$Z <- Z
  if(t == 1 | test) df$Y <- Y
  df$time <- t
  return(df)
}


#####################################################################################
#
#                                Simulation study
#
#####################################################################################


# ------------------------- Simulation function ------------------------------

# m1 : expected sample size at t = 1
# m2 : expected sample size at t = 2
# num : scenario number (1 ~ 5)
# var_group : names of subgroup variables
# var_conti : names of continuous predictor variables
# M : number of independent data generations
# K : number of subgroups
# R : number of folds in cross-validation in KLIEP (if kliep_cv = T)
# alpha : desired miscoverage level
# methods : list of methods to use (subset of "ar_m", "ar_g", "cqr_m", "cqr_g")
# sig_without_cv : sigma value for KLIEP method if kliep_cv = T
# kliep_cv : whether to use cross-validation in KLIEP procedure

DA_WCP_simul <- function(m1, m2, num, M = 100, K = 5, R = 4, alpha = 0.2, n_alpha = 5, rvec = F, 
                       rvec_marg = NA, rvec_group = NA, sig_wo_cv = 1.2, kliep_cv = F){
  
  # Simulation result
  ar_m_da_wcp <- list()
  ar_m_gwcp <- list()
  ar_g_da_wcp <- list()
  ar_g_gwcp <- list()
  cqr_m_da_wcp <- list()
  cqr_m_gwcp <- list()
  cqr_g_da_wcp <- list()
  cqr_g_gwcp <- list()
  
  if(!rvec) rvec_marg <- list()
  if(!rvec) rvec_group <- list()
  
  var_x = c("X1", "X2", "X3")
  var_y = c("Y")
  var_z = c()
  var_g = c("Z")
  var_h = c()
  
  for(m in 1:M){
    cat("------------------------------------------------------\n")
    cat("             (S", num, ") - ",  m, "th data progressing          \n", sep = "")
    cat("------------------------------------------------------\n")
    
    # Data generation
    set.seed(123+m)
    data1 <- generate_data(m1, t = 1, num)
    data2 <- generate_data(m2, t = 2, num)
    data_test <- generate_data(5000, t = 2, num, test = T)
    
    data1_y = data1[, c(var_x, var_y, var_z, var_g, var_h)] # t = 1 data
    data2_x = data2[, c(var_x, var_z, var_g, var_h)] # t = 2 data
    test <- data_test[, c(var_x, var_y, var_z, var_g, var_h)] # t = 2 test data
    n1 = nrow(data1_y)
    
    # Empirical coverage, average length
    in_intervals = vector(length = nrow(test))
    length_vec = vector(length = nrow(test))
    
    
    ######################### Marginal methods ##############################
    
    split_tra = sort(sample(n1, n1%/%2))
    split_cal = setdiff(1:n1, split_tra)
    train = data1_y[split_tra, ]
    cal = data1_y[split_cal, ]
    
    # -------------------------- AR_m_da_wcp ---------------------------- #
    if(!rvec) da_wcp_ar_m = DA_WCP(train, cal, train, data2_x, test, var_x, var_y, var_z, var_g, var_h,
                                   val_type = "marg", score_type = "ar", p2_vec = target_g, simul = T)
    if(rvec) da_wcp_ar_m = DA_WCP(train, cal, train, data2_x, test, var_x, var_y, var_z, var_g, var_h,
                                  val_type = "marg", score_type = "ar", p2_vec = target_g, simul = T,
                                  rvec = T, rvec_cal = rvec_marg[[m]])
    
    in_intervals = (test[[var_y]] >= da_wcp_ar_m$result_mat[, 1]) & (test[[var_y]] <= da_wcp_ar_m$result_mat[, 2])
    ar_m_da_wcp$eval.vec[m] = mean(in_intervals)
    ar_m_da_wcp$length.vec[m] = mean(da_wcp_ar_m$result_mat[, 2] - da_wcp_ar_m$result_mat[, 1])
    
    rvec_marg[[m]] = da_wcp_ar_m$rvec_cal
    
    # -------------------------- AR_m_gwcp ---------------------------- #
    gwcp_ar_m = DA_WCP(train, cal, train, data2_x, test, var_x, var_y, var_z, var_g, var_h,
                       gwcp = T, val_type = "marg", score_type = "ar", p2_vec = target_g, simul = T)
    
    in_intervals = (test[[var_y]] >= gwcp_ar_m$result_mat[, 1]) & (test[[var_y]] <= gwcp_ar_m$result_mat[, 2])
    ar_m_gwcp$eval.vec[m] = mean(in_intervals)
    ar_m_gwcp$length.vec[m] = mean(gwcp_ar_m$result_mat[, 2] - gwcp_ar_m$result_mat[, 1])

    # -------------------------- CQR_m_da_wcp ---------------------------- #
    da_wcp_cqr_m = DA_WCP(train, cal, train, data2_x, test, var_x, var_y, var_z, var_g, var_h,
                          val_type = "marg", score_type = "cqr", p2_vec = target_g, simul = T,
                          rvec = T, rvec_cal = rvec_marg[[m]])
    
    in_intervals = (test[[var_y]] >= da_wcp_cqr_m$result_mat[, 1]) & (test[[var_y]] <= da_wcp_cqr_m$result_mat[, 2])
    cqr_m_da_wcp$eval.vec[m] = mean(in_intervals)
    cqr_m_da_wcp$length.vec[m] = mean(da_wcp_cqr_m$result_mat[, 2] - da_wcp_cqr_m$result_mat[, 1])
    
    # -------------------------- CQR_m_gwcp ---------------------------- #
    gwcp_cqr_m = DA_WCP(train, cal, train, data2_x, test, var_x, var_y, var_z, var_g, var_h,
                        gwcp = T, val_type = "marg", score_type = "cqr", p2_vec = target_g, simul = T)
    
    in_intervals = (test[[var_y]] >= gwcp_cqr_m$result_mat[, 1]) & (test[[var_y]] <= gwcp_cqr_m$result_mat[, 2])
    cqr_m_gwcp$eval.vec[m] = mean(in_intervals)
    cqr_m_gwcp$length.vec[m] = mean(gwcp_cqr_m$result_mat[, 2] - gwcp_cqr_m$result_mat[, 1])
    
    
    ######################### Group-conditional methods ############################
    
    if(m == 1){
      ar_g_da_wcp$eval.mat = matrix(nrow = M, ncol = K)
      ar_g_gwcp$eval.mat = matrix(nrow = M, ncol = K)
      cqr_g_da_wcp$eval.mat = matrix(nrow = M, ncol = K)
      cqr_g_gwcp$eval.mat = matrix(nrow = M, ncol = K)
      
      ar_g_da_wcp$length.mat = matrix(nrow = M, ncol = K)
      ar_g_gwcp$length.mat = matrix(nrow = M, ncol = K)
      cqr_g_da_wcp$length.mat = matrix(nrow = M, ncol = K)
      cqr_g_gwcp$length.mat = matrix(nrow = M, ncol = K)
    }
    
    train = c()
    cal = c()
    
    for(k in 1:K){
      data1_y_g = data1_y[data1_y[[var_g]] == k, ]
      split_tra = sort(sample(nrow(data1_y_g), nrow(data1_y_g)%/%2))
      split_cal = setdiff(1:nrow(data1_y_g), split_tra)
      train = rbind(train, data1_y_g[split_tra, ])
      cal = rbind(cal, data1_y_g[split_cal, ])
    }
    
    # -------------------------- AR_g_da_wcp ---------------------------- #
    if(!rvec) da_wcp_ar_g = DA_WCP(train, cal, train, data2_x, test, var_x, var_y, var_z, var_g, var_h,
                                   val_type = "group", score_type = "ar", p2_vec = target_g, simul = T)
    if(rvec) da_wcp_ar_g = DA_WCP(train, cal, train, data2_x, test, var_x, var_y, var_z, var_g, var_h,
                                  val_type = "group", score_type = "ar", p2_vec = target_g, simul = T,
                                  rvec = T, rvec_cal = rvec_group[[m]])
    
    for(k in 1:K){
      in_intervals = (test[test[[var_g]] == k, var_y] >= da_wcp_ar_g$result_mat[test[[var_g]] == k, 1]) & 
        (test[test[[var_g]] == k, var_y] <= da_wcp_ar_g$result_mat[test[[var_g]] == k, 2])
      ar_g_da_wcp$eval.mat[m, k] = mean(in_intervals)
      ar_g_da_wcp$length.mat[m, k] = mean(da_wcp_ar_g$result_mat[test[[var_g]] == k, 2] -
                                            da_wcp_ar_g$result_mat[test[[var_g]] == k, 1])
    }
    
    rvec_group[[m]] = da_wcp_ar_g$rvec_cal
    
    # -------------------------- AR_g_gwcp ---------------------------- #
    gwcp_ar_g = DA_WCP(train, cal, train, data2_x, test, var_x, var_y, var_z, var_g, var_h,
                       gwcp = T, val_type = "group", score_type = "ar", p2_vec = target_g, simul = T)
    
    for(k in 1:K){
      in_intervals = (test[test[[var_g]] == k, var_y] >= gwcp_ar_g$result_mat[test[[var_g]] == k, 1]) & 
        (test[test[[var_g]] == k, var_y] <= gwcp_ar_g$result_mat[test[[var_g]] == k, 2])
      ar_g_gwcp$eval.mat[m, k] = mean(in_intervals)
      ar_g_gwcp$length.mat[m, k] = mean(gwcp_ar_g$result_mat[test[[var_g]] == k, 2] - 
                                          gwcp_ar_g$result_mat[test[[var_g]] == k, 1])
    }
    
    # -------------------------- CQR_g_da_wcp ---------------------------- #
    da_wcp_cqr_g = DA_WCP(train, cal, train, data2_x, test, var_x, var_y, var_z, var_g, var_h,
                          val_type = "group", score_type = "cqr", p2_vec = target_g, simul = T,
                          rvec = T, rvec_cal = rvec_group[[m]])
    
    for(k in 1:K){
      in_intervals = (test[test[[var_g]] == k, var_y] >= da_wcp_cqr_g$result_mat[test[[var_g]] == k, 1]) & 
        (test[test[[var_g]] == k, var_y] <= da_wcp_cqr_g$result_mat[test[[var_g]] == k, 2])
      cqr_g_da_wcp$eval.mat[m, k] = mean(in_intervals)
      cqr_g_da_wcp$length.mat[m, k] = mean(da_wcp_cqr_g$result_mat[test[[var_g]] == k, 2] - 
                                             da_wcp_cqr_g$result_mat[test[[var_g]] == k, 1])
    }
    
    # -------------------------- CQR_g_gwcp ---------------------------- #
    gwcp_cqr_g = DA_WCP(train, cal, train, data2_x, test, var_x, var_y, var_z, var_g, var_h,
                       gwcp = T, val_type = "group", score_type = "cqr", p2_vec = target_g, simul = T)
    
    for(k in 1:K){
      in_intervals = (test[test[[var_g]] == k, var_y] >= gwcp_cqr_g$result_mat[test[[var_g]] == k, 1]) & 
        (test[test[[var_g]] == k, var_y] <= gwcp_cqr_g$result_mat[test[[var_g]] == k, 2])
      cqr_g_gwcp$eval.mat[m, k] = mean(in_intervals)
      cqr_g_gwcp$length.mat[m, k] = mean(gwcp_cqr_g$result_mat[test[[var_g]] == k, 2] - 
                                           gwcp_cqr_g$result_mat[test[[var_g]] == k, 1])
    }

    cat("------------------ Empirical Coverage -------------------\n")
    cat("DA_WCP_AR_m :", round(mean(ar_m_da_wcp$eval.vec[1:m]), 3), "\n")
    cat("GWCP_AR_m :", round(mean(ar_m_gwcp$eval.vec[1:m]), 3), "\n")
    cat("DA_WCP_CQR_m :", round(mean(cqr_m_da_wcp$eval.vec[1:m]), 3), "\n")
    cat("GWCP_CQR_m :", round(mean(cqr_m_gwcp$eval.vec[1:m]), 3), "\n")
    
    cat("--------------- Average Interval Length ----------------\n")
    cat("DA_WCP_AR_m :", round(mean(ar_m_da_wcp$length.vec[1:m]), 3), "\n")
    cat("GWCP_AR_m :", round(mean(ar_m_gwcp$length.vec[1:m]), 3), "\n")
    cat("DA_WCP_CQR_m :", round(mean(cqr_m_da_wcp$length.vec[1:m]), 3), "\n")
    cat("GWCP_CQR_m :", round(mean(cqr_m_gwcp$length.vec[1:m]), 3), "\n")
  }
  
  eval.vec <- list(
    ar_m_da_wcp$eval.vec,
    ar_m_gwcp$eval.vec,
    cqr_m_da_wcp$eval.vec,
    cqr_m_gwcp$eval.vec
  )
  
  length.vec <- list(
    ar_m_da_wcp$length.vec,
    ar_m_gwcp$length.vec,
    cqr_m_da_wcp$length.vec,
    cqr_m_gwcp$length.vec
  )
  
  eval.mat <- list(
    ar_g_da_wcp$eval.mat, 
    ar_g_gwcp$eval.mat, 
    cqr_g_da_wcp$eval.mat, 
    cqr_g_gwcp$eval.mat
  )
  
  length.mat <- list(
    ar_g_da_wcp$length.mat, 
    ar_g_gwcp$length.mat, 
    cqr_g_da_wcp$length.mat, 
    cqr_g_gwcp$length.mat
  )

  df_m <- data.frame(
    coverage = do.call(base::c, eval.vec),
    length = do.call(base::c, length.vec),
    method = factor(rep(c("DA_WCP_AR_m", "GWCP_AR_m", "DA_WCP_CQR_m", "GWCP_CQR_m"), each = M))
  )
  
  df_g <- data.frame(
    coverage = do.call(base::c, eval.mat),
    length = do.call(base::c, length.mat),
    method = factor(rep(c("DA_WCP_AR_g", "GWCP_AR_g", "DA_WCP_CQR_g", "GWCP_CQR_g"), each = M*K))
  )
  
  return(list(df_m = df_m, df_g = df_g, rvec_marg = rvec_marg, rvec_group = rvec_group))
}

M = 100
m1 = 5000
m2 = 5000

for(num in 1:5){
  sim_num <- DA_WCP_simul(m1, m2, num, M = M, rvec = F)
  df_m[[num]] = sim_num$df_m
  df_g[[num]] = sim_num$df_g
  rvec_marg[[num]] = sim_num$rvec_marg
  rvec_group[[num]] = sim_num$rvec_group
  
  setwd("C://Users//lchk0//OneDrive//바탕 화면//Research//DA-WCP//Code")
  save.image(file = "DA_WCP_main_simul.RData")
}


#####################################################################################
#
#                         Visualization of the results
#
#####################################################################################

# Add new column
df_m[[1]]$sc <- "(S1)"
df_m[[2]]$sc <- "(S2)"
df_m[[3]]$sc <- "(S3)"
df_m[[4]]$sc <- "(S4)"
df_m[[5]]$sc <- "(S5)"

df_g[[1]]$sc <- "(S1)"
df_g[[2]]$sc <- "(S2)"
df_g[[3]]$sc <- "(S3)"
df_g[[4]]$sc <- "(S4)"
df_g[[5]]$sc <- "(S5)" 


# Combine into a single dataset
df_all_m <- bind_rows(df_m[[1]], df_m[[2]], df_m[[3]], df_m[[4]], df_m[[5]])

df_all_g <- bind_rows(df_g[[1]], df_g[[2]], df_g[[3]], df_g[[4]], df_g[[5]])


df_all_m$method <- factor(df_all_m$method, levels = c("GWCP_CQR_m", "GWCP_AR_m", "DA_WCP_CQR_m", "DA_WCP_AR_m"))


setwd("C://Users//lchk0//OneDrive//바탕 화면//Research//DA-WCP//Figures")
pdf(file = "simulation_coverage_m.pdf", width = 13, height = 3.5)
ggplot(df_all_m, aes(x = coverage, y = method, fill = method)) +
  geom_boxplot(outlier.size = 0.8) +
  geom_vline(xintercept = 0.8, color = "red", linewidth = 0.7) +
  facet_wrap(~ sc, nrow = 1) +
  scale_x_continuous(limits = c(0.6, 0.9)) +
  xlab("Empirical coverage") +
  ylab("") +
  scale_fill_manual(
    values = c(
      "DA_WCP_AR_m"   = "#9ecae1",  
      "GWCP_AR_m"   = "#fdae6b",  
      "DA_WCP_CQR_m"  = "#3182bd",
      "GWCP_CQR_m"  = "#fd8d3c"   
    )
  ) +
  theme_test(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "gray90"),
    axis.text.y = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  )
dev.off()


pdf(file = "simulation_length_m.pdf", width = 13, height = 3.5)
ggplot(df_all_m, aes(x = length, y = method, fill = method)) +
  geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~sc, nrow = 1) +
  scale_x_continuous(limits = c(4, 14)) +
  xlab("Average interval length") +
  ylab("") +
  scale_fill_manual(
    values = c(
      "DA_WCP_AR_m"   = "#9ecae1",  
      "GWCP_AR_m"   = "#fdae6b",  
      "DA_WCP_CQR_m"  = "#3182bd",
      "GWCP_CQR_m"  = "#fd8d3c"   
    )
  ) +
  theme_test(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "gray90"),
    axis.text.y = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  )
dev.off()




df_all_g$sc <- factor(df_all_g$sc, levels = c("(S1)", "(S2)", "(S3)", "(S4)", "(S5)"))
df_all_g$Subgroup <- as.factor(rep(c(rep('1', M), rep('2', M), rep('3', M), rep('4', M), rep('5', M)), 20))


pdf(file = "simulation_coverage_g.pdf", width = 11, height = 6)
ggplot(df_all_g, aes(x = Subgroup, y = coverage, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.size = 0.8) +
  facet_wrap(~sc, nrow = 2) +
  geom_hline(yintercept = 0.8,  color = "red", linewidth = 0.6) +
  ylim(0.55, 0.95) +
  labs(x = "Group", y = "Group-wise coverage") +
  scale_fill_manual(
    values = c(
      "DA_WCP_AR_g" = "#9ecae1",   # 연파랑
      "DA_WCP_CQR_g" = "#08519c",  # 진파랑
      "GWCP_AR_g" = "#fdae6b",       # 연주황
      "GWCP_CQR_g" = "#e6550d"       # 진주황
    ),
    name = "Method"
  ) +
  theme_test(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "gray95"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),  
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "left"
  )
dev.off()
















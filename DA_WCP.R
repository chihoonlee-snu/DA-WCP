library(nnet)
library(dplyr)
library(xgboost)
library(ranger)
library(quantreg)


#################### Code for DA-WCP (Design-adapted WCP) ####################

############# 1. DA_WCP (main function) ##############

# A function to compute WCP-type threshold
return_cutoff <- function(wtvec, scvec, a){
  ord = order(scvec)
  svec = sort(scvec)
  wvec = wtvec[ord]
  tot = 0
  for(i in 1:length(scvec)){
    tot = tot + wvec[i]
    if(tot >= 1-a) break
  }
  return(svec[i])
}

# KLIEP implementation function
setwd("C://Users//lchk0//OneDrive//바탕 화면//Research//DA-WCP//Code")
source("KLIEP.R")

# For categorical outcome, group-conditional methods
matrix_adjust <- function(mat, ref = factor(0:4), colmat){
  ret = matrix(0, nrow = nrow(mat), ncol = length(ref))
  colnames(ret) = ref
  ret[, ref %in% colmat] = mat
  return(ret)
}


# For continuous outcome, returns nrow(dtest)*2 (left, right endpoints)
# For categorical outcome, returns nrow(dtest)*M (0 or 1 for each category)
# val_type : marg, group
# score_type : ar, cqr (for continuous outcomes) / rf, nn, xgb (for categorical outcomes)

DA_WCP <- function(dtra, dcal, dttra, d2, dtest, var_x, var_y, var_z, var_g, var_h, val_type, score_type, 
                   gwcp = F, var_wt_x = NA, var_wt_y = NA, alpha = 0.2, p2_vec = NA, n_alpha = 5, 
                   rvec = F, rvec_cal = NA, R = 4, sig_wo_cv = 1.2, kliep_cv = F, simul = F, eval = T){
  
  # Define new variable (z_tilde)
  var_zt = c(var_z, var_g)
  
  # (alpha_l, alpha_h) grid for CQR score
  alpha_l = seq(from = 0, to = 0.5, length.out = (n_alpha+2))[-c(1, (n_alpha+2))]
  alpha_h = seq(from = 0.5, to = 1.0, length.out = (n_alpha+2))[-c(1, (n_alpha+2))]
  alpha_grid = as.matrix(expand.grid(alpha_l, alpha_h))
  
  ## -------------------- counts --------------------
  K_zt <- if (length(var_zt) == 0) 1L else nrow(dtra %>% dplyr::count(dplyr::across(dplyr::all_of(var_zt))))
  K_z  <- if (length(var_z)  == 0) 1L else nrow(dtra %>% dplyr::count(dplyr::across(dplyr::all_of(var_z))))
  K_g  <- if (length(var_g)  == 0) 1L else nrow(dtra %>% dplyr::count(dplyr::across(dplyr::all_of(var_g))))
  H    <- if (length(var_h)  == 0) 1L else nrow(dtra %>% dplyr::count(dplyr::across(dplyr::all_of(var_h))))
  
  ## -------------------- encode zt --------------------
  if (length(var_zt) == 0 || K_zt <= 1) {
    dtra$zt  <- 1L; dcal$zt  <- 1L; dttra$zt <- 1L; d2$zt <- 1L; dtest$zt <- 1L
    K_zt <- 1L
  } else {
    lev_zt <- levels(interaction(dtra[, var_zt, drop = FALSE], drop = TRUE))
    dtra$zt  <- as.integer(factor(interaction(dtra[,  var_zt, drop = FALSE], drop = TRUE), levels = lev_zt))
    dcal$zt  <- as.integer(factor(interaction(dcal[,  var_zt, drop = FALSE], drop = TRUE), levels = lev_zt))
    dttra$zt <- as.integer(factor(interaction(dttra[, var_zt, drop = FALSE], drop = TRUE), levels = lev_zt))
    d2$zt    <- as.integer(factor(interaction(d2[,    var_zt, drop = FALSE], drop = TRUE), levels = lev_zt))
    dtest$zt <- as.integer(factor(interaction(dtest[, var_zt, drop = FALSE], drop = TRUE), levels = lev_zt))
  }
  
  ## -------------------- encode z --------------------
  if (length(var_z) == 0 || K_z <= 1) {
    dtra$z  <- 1L; dcal$z  <- 1L; dttra$z <- 1L; d2$z <- 1L; dtest$z <- 1L
    K_z <- 1L
  } else {
    lev_z <- levels(interaction(dtra[, var_z, drop = FALSE], drop = TRUE))
    dtra$z  <- as.integer(factor(interaction(dtra[,  var_z, drop = FALSE], drop = TRUE), levels = lev_z))
    dcal$z  <- as.integer(factor(interaction(dcal[,  var_z, drop = FALSE], drop = TRUE), levels = lev_z))
    dttra$z <- as.integer(factor(interaction(dttra[, var_z, drop = FALSE], drop = TRUE), levels = lev_z))
    d2$z    <- as.integer(factor(interaction(d2[,    var_z, drop = FALSE], drop = TRUE), levels = lev_z))
    dtest$z <- as.integer(factor(interaction(dtest[, var_z, drop = FALSE], drop = TRUE), levels = lev_z))
  }
  
  ## -------------------- encode g --------------------
  if (length(var_g) == 0 || K_g <= 1) {
    dtra$g  <- 1L; dcal$g  <- 1L; dttra$g <- 1L; d2$g <- 1L; dtest$g <- 1L
    K_g <- 1L
  } else {
    lev_g <- levels(interaction(dtra[, var_g, drop = FALSE], drop = TRUE))
    dtra$g  <- as.integer(factor(interaction(dtra[,  var_g, drop = FALSE], drop = TRUE), levels = lev_g))
    dcal$g  <- as.integer(factor(interaction(dcal[,  var_g, drop = FALSE], drop = TRUE), levels = lev_g))
    dttra$g <- as.integer(factor(interaction(dttra[, var_g, drop = FALSE], drop = TRUE), levels = lev_g))
    d2$g    <- as.integer(factor(interaction(d2[,    var_g, drop = FALSE], drop = TRUE), levels = lev_g))
    dtest$g <- as.integer(factor(interaction(dtest[, var_g, drop = FALSE], drop = TRUE), levels = lev_g))
  }
  
  ## -------------------- encode st --------------------
  if (length(var_h) == 0 || H <= 1) {
    dtra$st  <- 1L; dcal$st  <- 1L; dttra$st <- 1L; d2$st <- 1L; dtest$st <- 1L
    H <- 1L
  } else {
    lev_st <- levels(interaction(dtra[, var_h, drop = FALSE], drop = TRUE))
    dtra$st  <- as.integer(factor(interaction(dtra[,  var_h, drop = FALSE], drop = TRUE), levels = lev_st))
    dcal$st  <- as.integer(factor(interaction(dcal[,  var_h, drop = FALSE], drop = TRUE), levels = lev_st))
    dttra$st <- as.integer(factor(interaction(dttra[, var_h, drop = FALSE], drop = TRUE), levels = lev_st))
    d2$st    <- as.integer(factor(interaction(d2[,    var_h, drop = FALSE], drop = TRUE), levels = lev_st))
    dtest$st <- as.integer(factor(interaction(dtest[, var_h, drop = FALSE], drop = TRUE), levels = lev_st))
  }
  
  
  # Regression formula
  formula_m <- as.formula(paste(var_y, "~", paste(c(var_x, var_z, var_g), collapse = " + ")))
  formula_g <- as.formula(paste(var_y, "~", paste(c(var_x, var_z), collapse = " + ")))
  
  ### Estimating subgroup-proportions
  q1h_mat = matrix(nrow = H, ncol = K_g) # Training demographic proportion for each stratum
  q2h_mat = matrix(nrow = H, ncol = K_g) # t = 2 dataset demographic proportion for each stratum (for evaluation)
  if(!simul) p2_vec = vector(length = K_g) # Target demographic proportion
  v1_mat = matrix(nrow = K_g, ncol = K_z)
  v2_mat = matrix(nrow = K_g, ncol = K_z)
  
  if(!simul) for(k_g in 1:K_g) p2_vec[k_g] = sum(d2[[var_wt_x]][d2$g == k_g])/sum(d2[[var_wt_x]])
  
  # Only for marginal method
  if(val_type == "marg"){
    for(h in 1:H){
      if(!simul & eval){
        d2_y_h = d2[!is.na(d2[[var_y]]) & d2$st == h, ]
        q2h_mat[h, ] = proportions(table(factor(d2_y_h$g, levels = 1:K_g)))
      }
      if(!gwcp){
        dtra_h = dtra[dtra$st == h, ]
        q1h_mat[h, ] = proportions(table(factor(dtra_h$g, levels = 1:K_g)))
      }
    }
  }

  if(val_type == "marg" & gwcp) q1h_mat[1, ] = proportions(table(factor(dcal$g, levels = 1:K_g)))
  
  # For both marginal, group-conditional methods, only DA-WCP
  if(!gwcp){
    for(k_g in 1:K_g){
      dttra_g = dttra[dttra$g == k_g, ]
      d2_g = d2[d2$g == k_g, ]
      
      v1_mat[k_g, ] = proportions(table(factor(dttra_g$z, levels = 1:K_z)))
      v2_mat[k_g, ] = proportions(table(factor(d2_g$z, levels = 1:K_z)))
    }
  }
  
  ### Estimating subgroup-specific density-ratio functions (KLIEP)
  if(!rvec & !gwcp){
    rvec_cal = vector(length = nrow(dcal))
    for(k_zt in 1:K_zt){
      dttra_zt = dttra[dttra$zt == k_zt, var_x, drop = F]
      dcal_zt = dcal[dcal$zt == k_zt, var_x, drop = F]
      d2_zt = d2[d2$zt == k_zt, var_x, drop = F]
      
      cat(k_zt, "th group KLIEP progressing (n1 = ", nrow(dttra_zt), ", n2 = ", nrow(d2_zt), ") \n", sep = "")
      
      mu_vec = colMeans(d2_zt)
      sd_vec = apply(d2_zt, 2, sd)
      
      dttra_normalized = as.matrix(sweep(sweep(dttra_zt, 2, mu_vec, "-"), 2, sd_vec, "/"))
      dcal_normalized = as.matrix(sweep(sweep(dcal_zt, 2, mu_vec, "-"), 2, sd_vec, "/"))
      d2_normalized = as.matrix(sweep(sweep(d2_zt, 2, mu_vec, "-"), 2, sd_vec, "/"))
      rvec_cal[dcal$zt == k_zt] = rvec_fun(dcal_normalized, dttra_normalized, d2_normalized, R = R, sigma = sig_wo_cv, cv = kliep_cv)
    }
  }
  
  ### Assigning estimated WCP weights to calibration data
  weight_vec = rep(1, length = nrow(dcal))
  
  if(val_type == "marg" & !gwcp){
    for(k_g in 1:K_g){
      for(k_z in 1:K_z){
        for(h in 1:H){
          temp = rvec_cal[dcal$g == k_g & dcal$z == k_z & dcal$st == h]
          weight_vec[dcal$g == k_g & dcal$z == k_z & dcal$st == h] = 
            temp*(v2_mat[k_g, k_z]/v1_mat[k_g, k_z])*(p2_vec[k_g]/q1h_mat[h, k_g])
        }
      }
    }
    weight_vec = proportions(weight_vec)
  }
  
  if(val_type == "marg" & gwcp){
    for(k_g in 1:K_g) weight_vec[dcal$g == k_g] = p2_vec[k_g]/q1h_mat[1, k_g]
    weight_vec = proportions(weight_vec)
  }
  
  if(val_type == "group" & !gwcp){
    for(k_g in 1:K_g){
      for(k_z in 1:K_z){
        temp = rvec_cal[dcal$g == k_g & dcal$z == k_z]
        weight_vec[dcal$g == k_g & dcal$z == k_z] = temp*v2_mat[k_g, k_z]/v1_mat[k_g, k_z]
      }
      weight_vec[dcal$g == k_g] = proportions(weight_vec[dcal$g == k_g])
    }
  }
  
  if(val_type == "group" & gwcp){
    for(k_g in 1:K_g) weight_vec[dcal$g == k_g] = proportions(weight_vec[dcal$g == k_g])
  }

  ### Run DA-WCP (or GWCP)
  if(val_type == "marg"){
    
    ### Continuous outcome
    
    # (1) Absolute residual score
    if(score_type == "ar"){
      if(simul) wls = lm(formula_m, data = dtra)
      if(!simul) wls = lm(formula_m, data = dtra, weights = dtra[[var_wt_y]])
      pred_cal = as.vector(predict(wls, dcal))
      pred_test = as.vector(predict(wls, dtest))
      score_vec = abs(pred_cal - dcal[[var_y]])
      
      # WCP-type threshold
      t_alpha = return_cutoff(weight_vec, score_vec, alpha)
      
      result_mat = as.data.frame(cbind(pred_test - t_alpha, pred_test + t_alpha))
      colnames(result_mat) = c("left_end", "right_end")
    }
    
    # (2) CQR score
    if(score_type == "cqr"){
      score_mat = matrix(nrow = nrow(dcal), ncol = n_alpha**2)
      pred_test_l = matrix(nrow = nrow(dtest), ncol = n_alpha**2)
      pred_test_h = matrix(nrow = nrow(dtest), ncol = n_alpha**2)
      
      for(ind in 1:n_alpha**2){
        if(simul){
          wls_l = rq(formula_m, tau = alpha_grid[ind, 1], data = dtra)
          wls_h = rq(formula_m, tau = alpha_grid[ind, 2], data = dtra)
        }
        if(!simul){
          wls_l = rq(formula_m, tau = alpha_grid[ind, 1], data = dtra, weights = dtra[[var_wt_y]])
          wls_h = rq(formula_m, tau = alpha_grid[ind, 2], data = dtra, weights = dtra[[var_wt_y]])
        }
        pred_l = as.vector(predict(wls_l, dcal))
        pred_h = as.vector(predict(wls_h, dcal))
        pred_test_l[, ind] = as.vector(predict(wls_l, dtest))
        pred_test_h[, ind] = as.vector(predict(wls_h, dtest))
        score_mat[, ind] = apply(cbind(pred_l - dcal[[var_y]], dcal[[var_y]] - pred_h), 1, max)
      }
      
      # WCP-type threshold for each (alpha_l, alpha_h)
      length_mat = matrix(nrow = nrow(dtest), ncol = n_alpha**2)
      length_vec = vector(length = n_alpha**2)
      t_alpha_vec = vector(length = n_alpha**2)
      for(ind in 1:n_alpha**2){
        t_alpha_vec[ind] = return_cutoff(weight_vec, score_mat[, ind], alpha)
        length_mat[, ind] = 2*t_alpha_vec[ind] + pred_test_h[, ind] - pred_test_l[, ind]
        length_vec[ind] = mean(length_mat[, ind])
      }
      
      # (alpha_l, alpha_h) that yields shortest average interval length
      ind_cqr = which.min(length_vec)
      selected_alpha = c(alpha_grid[ind_cqr, 1], alpha_grid[ind_cqr, 2])
      
      # WCP-type threshold
      t_alpha = t_alpha_vec[ind_cqr]
      
      result_mat = as.data.frame(cbind(pred_test_l[, ind_cqr] - t_alpha, pred_test_h[, ind_cqr] + t_alpha))
      colnames(result_mat) = c("left_end", "right_end")
    }
    
    if(score_type %in% c("ar", "cqr") & !simul & eval){
      in_interval = (dtest[[var_y]] <= result_mat[, 2]) & (dtest[[var_y]] >= result_mat[, 1])
      eval_weight = rep(1, nrow(dtest))
      for(k_g in 1:K_g){
        for(h in 1:H) eval_weight[dtest$g == k_g & dtest$st == h] = p2_vec[k_g]/q2h_mat[h, k_g]
      }
      coverage = mean(eval_weight*in_interval)
      avglen = mean(result_mat[, 2] - result_mat[, 1])
      return(list(result_mat = result_mat, coverage = coverage, avglen = avglen, rvec_cal = rvec_cal))
    }
    
    ### Categorical outcome
    if(score_type %in% c("xgb", "rf", "nn")){
      scaled_dtra_x <- scale(dtra[, var_x])
      
      scaled_dtra = dtra
      scaled_dtra[, var_x] <- scaled_dtra_x
      
      scaled_dcal = dcal
      scaled_dcal[, var_x] <- scale(dcal[, var_x],
                                    center = attr(scaled_dtra_x, "scaled:center"),
                                    scale = attr(scaled_dtra_x, "scaled:scale"))
      
      scaled_dtest = dtest
      scaled_dtest[, var_x] <- scale(dtest[, var_x],
                                     center = attr(scaled_dtra_x, "scaled:center"),
                                     scale = attr(scaled_dtra_x, "scaled:scale"))
      M = max(as.numeric(dtra[[var_y]]))
    }
    
    # (3) XGBoost score
    if(score_type == "xgb"){
      xgb_model_tra <- model.matrix(as.formula(paste(var_y, "~ . - 1")), scaled_dtra[, c(var_x, var_y, var_z, var_g)])
      xgb_y <- as.numeric(scaled_dtra[[var_y]]) - 1
      if(simul) xgb_dtra <- xgb.DMatrix(data = xgb_model_tra, label = xgb_y)
      if(!simul) xgb_dtra <- xgb.DMatrix(data = xgb_model_tra, label = xgb_y, weight = scaled_dtra[[var_wt_y]])
      fit <- xgboost(data = xgb_dtra, objective = "multi:softprob", num_class = M, nrounds = 1000, verbose = 0)
      
      xgb_model_cal <- model.matrix(as.formula(paste(var_y, "~ . - 1")), scaled_dcal[, c(var_x, var_y, var_z, var_g)])
      pred_cal <- predict(fit, newdata = xgb_model_cal)
      pred_mat <- matrix(pred_cal, ncol = M, byrow = TRUE)
      
      xgb_model_test <- model.matrix(as.formula(paste(var_y, "~ . - 1")), scaled_dtest[, c(var_x, var_y, var_z, var_g)])
      pred_test <- predict(fit, newdata = xgb_model_test)
      test_mat <- matrix(pred_test, ncol = M, byrow = TRUE)
    }
    
    # (4) Random forest score
    if(score_type == "rf"){
      if(!simul) fit <- ranger(formula = formula_m, data = scaled_dtra, probability = TRUE, 
                               classification = TRUE, verbose = FALSE, case.weights = scaled_dtra[[var_wt_y]])
      if(simul) fit <- ranger(formula = formula_m, data = scaled_dtra, probability = TRUE, 
                              classification = TRUE, verbose = FALSE)
      pred_mat <- predict(fit, scaled_dcal)$predictions
      test_mat <- predict(fit, scaled_dtest)$predictions
    }
    
    # (5) Neural net score
    if(score_type == "nn"){
      if(!simul) fit <- nnet(formula = formula_m, data = scaled_dtra, weights = scaled_dtra[[var_wt_y]], 
                             decay = 0.01, size = M, maxit = 10000, trace = F)
      if(simul) fit <- nnet(formula = formula_m, data = scaled_dtra, decay = 0.01, 
                             size = M, maxit = 10000, trace = F)
      pred_mat <- predict(fit, scaled_dcal, type = "raw")
      test_mat <- predict(fit, scaled_dtest, type = "raw")
    }
    
    # WCP-type threshold for categorical outcomes
    if(score_type %in% c("rf", "xgb", "nn")){
      score_vec = rep(1, nrow(scaled_dcal))
      result_mat = matrix(0, nrow = nrow(scaled_dtest), ncol = M)
      
      for (i in 1:nrow(scaled_dcal)) score_vec[i] = 1 - pred_mat[i, scaled_dcal[[var_y]][i]]
      t_alpha = return_cutoff(weight_vec, score_vec, alpha)
      
      for (j in 1:M) result_mat[, j] = as.numeric(test_mat[, j] >= 1 - t_alpha)
      
      if(!simul & eval){
        in_interval = result_mat[cbind(seq_len(nrow(result_mat)), as.numeric(dtest[[var_y]]))]
        eval_weight = rep(1, nrow(dtest))
        for(k_g in 1:K_g){
          for(h in 1:H) eval_weight[dtest$g == k_g & dtest$st == h] = p2_vec[k_g]/q2h_mat[h, k_g]
        }
        coverage = mean(eval_weight*in_interval)
        return(list(result_mat = result_mat, coverage = coverage, rvec_cal = rvec_cal))
      }
    }
  }
  
  if(val_type == "group"){
    
    if(score_type %in% c("ar", "cqr")) result_mat = matrix(nrow = nrow(dtest), ncol = 2)
    if(score_type %in% c("xgb", "rf", "nn")){
      M = max(as.numeric(dtra[[var_y]]))
      result_mat = matrix(0, nrow = nrow(dtest), ncol = M)
    }
    
    for(k_g in 1:K_g){
      dtra_g = dtra[dtra$g == k_g, ]
      dcal_g = dcal[dcal$g == k_g, ]
      dtest_g = dtest[dtest$g == k_g, ]
      weight_vec_g = weight_vec[dcal$g == k_g]
      
      ### Continuous outcome
      
      # (1) Absolute residual score
      if(score_type == "ar"){
        if(simul) wls_g = lm(formula_g, data = dtra_g)
        if(!simul) wls_g = lm(formula_g, data = dtra_g, weights = dtra_g[[var_wt_y]])
        pred_cal_g = as.vector(predict(wls_g, dcal_g))
        pred_test_g = as.vector(predict(wls_g, dtest_g))
        score_vec_g = abs(pred_cal_g - dcal_g[[var_y]])
        
        # WCP-type threshold
        t_alpha_g = return_cutoff(weight_vec_g, score_vec_g, alpha)
        
        result_mat[dtest$g == k_g, ] = cbind(pred_test_g - t_alpha_g, pred_test_g + t_alpha_g)
      }
      
      # (2) CQR score
      if(score_type == "cqr"){
        score_mat_g = matrix(nrow = nrow(dcal_g), ncol = n_alpha**2)
        pred_test_l_g = matrix(nrow = nrow(dtest_g), ncol = n_alpha**2)
        pred_test_h_g = matrix(nrow = nrow(dtest_g), ncol = n_alpha**2)
        
        for(ind in 1:n_alpha**2){
          if(simul){
            wls_l_g = rq(formula_g, tau = alpha_grid[ind, 1], data = dtra_g)
            wls_h_g = rq(formula_g, tau = alpha_grid[ind, 2], data = dtra_g)
          }
          if(!simul){
            wls_l_g = rq(formula_g, tau = alpha_grid[ind, 1], data = dtra_g, weights = dtra_g[[var_wt_y]])
            wls_h_g = rq(formula_g, tau = alpha_grid[ind, 2], data = dtra_g, weights = dtra_g[[var_wt_y]])
          }
          pred_l_g = as.vector(predict(wls_l_g, dcal_g))
          pred_h_g = as.vector(predict(wls_h_g, dcal_g))
          pred_test_l_g[, ind] = as.vector(predict(wls_l_g, dtest_g))
          pred_test_h_g[, ind] = as.vector(predict(wls_h_g, dtest_g))
          score_mat_g[, ind] = apply(cbind(pred_l_g - dcal_g[[var_y]], dcal_g[[var_y]] - pred_h_g), 1, max)
        }
        
        # WCP-type threshold for each (alpha_l, alpha_h)
        length_mat_g = matrix(nrow = nrow(dtest_g), ncol = n_alpha**2)
        length_vec_g = vector(length = n_alpha**2)
        t_alpha_vec = vector(length = n_alpha**2)
        for(ind in 1:n_alpha**2){
          t_alpha_vec[ind] = return_cutoff(weight_vec_g, score_mat_g[, ind], alpha)
          length_mat_g[, ind] = 2*t_alpha_vec[ind] + pred_test_h_g[, ind] - pred_test_l_g[, ind]
          length_vec_g[ind] = mean(length_mat_g[, ind])
        }
        
        # (alpha_l, alpha_h) that yields shortest average interval length
        ind_cqr_g = which.min(length_vec_g)
        selected_alpha_g = c(alpha_grid[ind_cqr_g, 1], alpha_grid[ind_cqr_g, 2])
        
        # WCP-type threshold
        t_alpha_g = t_alpha_vec[ind_cqr_g]
        
        result_mat[dtest$g == k_g, ] = cbind(pred_test_l_g[, ind_cqr_g] - t_alpha_g, 
                                             pred_test_h_g[, ind_cqr_g] + t_alpha_g)
      }
      
      ### Categorical outcome
      if(score_type %in% c("xgb", "rf", "nn")){
        scaled_dtra_x_g <- scale(dtra_g[, var_x])
        
        scaled_dtra_g = dtra_g
        scaled_dtra_g[, var_x] <- scaled_dtra_x_g
        
        scaled_dcal_g = dcal_g
        scaled_dcal_g[, var_x] <- scale(dcal_g[, var_x],
                                        center = attr(scaled_dtra_x_g, "scaled:center"),
                                        scale = attr(scaled_dtra_x_g, "scaled:scale"))
        
        scaled_dtest_g = dtest_g
        scaled_dtest_g[, var_x] <- scale(dtest_g[, var_x],
                                         center = attr(scaled_dtra_x_g, "scaled:center"),
                                         scale = attr(scaled_dtra_x_g, "scaled:scale"))
        M = max(as.numeric(dtra[[var_y]]))
      }
      
      # (3) XGBoost score
      if(score_type == "xgb"){
        xgb_model_tra_g <- model.matrix(as.formula(paste(var_y, "~ . - 1")), scaled_dtra_g[, c(var_x, var_y, var_z)])
        xgb_y_g <- as.numeric(scaled_dtra_g[[var_y]]) - 1
        if(simul) xgb_dtra_g <- xgb.DMatrix(data = xgb_model_tra_g, label = xgb_y_g)
        if(!simul) xgb_dtra_g <- xgb.DMatrix(data = xgb_model_tra_g, label = xgb_y_g, weight = scaled_dtra_g[[var_wt_y]])
        fit <- xgboost(data = xgb_dtra_g, objective = "multi:softprob", num_class = M, nrounds = 1000, verbose = 0)
        
        xgb_model_cal_g <- model.matrix(as.formula(paste(var_y, "~ . - 1")), scaled_dcal_g[, c(var_x, var_y, var_z)])
        pred_cal_g <- predict(fit, newdata = xgb_model_cal_g)
        pred_mat_g <- matrix(pred_cal_g, ncol = M, byrow = TRUE)
        
        xgb_model_test_g <- model.matrix(as.formula(paste(var_y, "~ . - 1")), scaled_dtest_g[, c(var_x, var_y, var_z)])
        pred_test_g <- predict(fit, newdata = xgb_model_test_g)
        test_mat_g <- matrix(pred_test_g, ncol = M, byrow = TRUE)
      }
      
      # (4) Random forest score
      if(score_type == "rf"){
        colmat = unique(scaled_dtra_g[[var_y]])
        if(!simul) fit <- ranger(formula = formula_g, data = scaled_dtra_g, probability = TRUE, 
                                 classification = TRUE, verbose = FALSE, case.weights = scaled_dtra_g[[var_wt_y]])
        if(simul) fit <- ranger(formula = formula_g, data = scaled_dtra_g, probability = TRUE, 
                                classification = TRUE, verbose = FALSE)
        pred_mat_g <- matrix_adjust(predict(fit, scaled_dcal_g)$predictions, colmat = colmat)
        test_mat_g <- matrix_adjust(predict(fit, scaled_dtest_g)$predictions, colmat = colmat)
      }
      
      # (5) Neural net score
      if(score_type == "nn"){
        colmat = unique(scaled_dtra_g[[var_y]])
        if(!simul) fit <- nnet(formula = formula_g, data = scaled_dtra_g, weights = scaled_dtra_g[[var_wt_y]], 
                               decay = 0.01, size = M, maxit = 10000, trace = F)
        if(simul) fit <- nnet(formula = formula_g, data = scaled_dtra_g, decay = 0.01, 
                              size = M, maxit = 10000, trace = F)
        pred_mat_g <- matrix_adjust(predict(fit, scaled_dcal_g, type = "raw"), colmat = colmat)
        test_mat_g <- matrix_adjust(predict(fit, scaled_dtest_g, type = "raw"), colmat = colmat)
      }
      
      # WCP-type threshold for categorical outcomes
      if(score_type %in% c("rf", "xgb", "nn")){
        score_vec_g = rep(1, nrow(scaled_dcal_g))
        for(i in 1:nrow(scaled_dcal_g)) score_vec_g[i] = 1 - pred_mat_g[i, scaled_dcal_g[[var_y]][i]]
        t_alpha_g = return_cutoff(weight_vec_g, score_vec_g, alpha)
        for(j in 1:M) result_mat[dtest$g == k_g, j] = as.numeric(test_mat_g[, j] >= 1 - t_alpha_g)
      }
      
    }
    
    if(score_type %in% c("ar", "cqr") & !simul & eval){
      coverage = vector(length = K_g)
      avglen = vector(length = K_g)
      in_interval = (dtest[[var_y]] <= result_mat[, 2]) & (dtest[[var_y]] >= result_mat[, 1])
      for(k_g in 1:K_g){
        coverage[k_g] = mean(in_interval[dtest$g == k_g])
        avglen[k_g] = mean(result_mat[dtest$g == k_g, 2] - result_mat[dtest$g == k_g, 1])
      }
      return(list(result_mat = result_mat, coverage = coverage, avglen = avglen, rvec_cal = rvec_cal))
    }
    
    if(score_type %in% c("xgb", "rf", "nn") & !simul & eval){
      in_interval = result_mat[cbind(seq_len(nrow(result_mat)), as.numeric(dtest[[var_y]]))]
      coverage = vector(length = K_g)
      for(k_g in 1:K_g) coverage[k_g] = mean(in_interval[dtest$g == k_g])
      return(list(result_mat = result_mat, coverage = coverage, rvec_cal = rvec_cal))
    }
  }
  
  return(list(result_mat = result_mat, rvec_cal = rvec_cal))
}







### Direct importance estimation (KLIEP)

library(mvtnorm)

g_b_fun <- function(sigma, Xtra, Xtest, Xcenter){
  b = nrow(Xcenter)
  gamma = vector(length = b)
  beta = matrix(nrow = nrow(Xtest), ncol = b)
  for(l in 1:b){
    gamma[l] = sum(apply(Xtra, 1, function(z) return(k_fun(z, Xcenter[l, ], sigma))))
    beta[, l] = apply(Xtest, 1, function(z) return(k_fun(z, Xcenter[l, ], sigma)))
  }
  return(list(gamma = gamma, beta = beta))
}

k_fun <- function(x, y, sigma){
  return(exp(-sum((x-y)**2)/(2*sigma**2)))
}

obj_fun <- function(alpha, sigma, beta){
  ret = sum(log(as.vector(beta%*%alpha)))
  return(ret/nrow(beta))
}

obj_fun_cv <- function(alpha, sigma, Xtestr, Xcenter){
  beta = matrix(nrow = nrow(Xtestr), ncol = length(alpha))
  for(l in 1:length(alpha)) beta[, l] = 
      apply(Xtestr, 1, function(z) return(k_fun(z, Xcenter[l, ], sigma)))
  ret =  sum(log(as.vector(beta%*%alpha)))
  return(ret/nrow(beta))
}

initial_point <- function(sigma, beta, gamma, N, ntrar){
  set.seed(12345)
  b = length(gamma)
  x = runif(N*b, 0, 2)
  alpha_mat = matrix(x, nrow = N)
  for(i in 1:N) alpha_mat[i,] = alpha_mat[i,]*ntrar/sum(alpha_mat[i,]*gamma)
  obj_vec = vector(length = N)
  for(i in 1:N) obj_vec[i] = obj_fun(alpha_mat[i,], sigma, beta)
  return(alpha_mat[which.max(obj_vec), ])
}

check_constraint <- function(gamma, alpha, ntrar){
  if(abs(sum(gamma*alpha) - ntrar) < 1e-6 & min(alpha) >= 0) return("OK")
  else return("NO")
}

update_alpha <- function(eps, alpha, beta, gamma){
  b = length(gamma)
  grad_vec = rep(0, b)
  strata = alpha > 0
  for(j in 1:nrow(beta)) grad_vec = grad_vec + beta[j, ]/sum(alpha*beta[j, ])
  grad_vec = grad_vec/nrow(beta)
  
  gamma_strata = gamma[strata]
  alpha_strata = alpha[strata]
  grad_vec_pos = grad_vec[strata]
  grad_vec_pos = grad_vec_pos - sum(grad_vec_pos*gamma_strata)*gamma_strata/sum(gamma_strata**2)
  grad_vec[strata] = grad_vec_pos
  grad_vec[!strata] = 0
  
  detect_change = vector(length = length(grad_vec_pos))
  neg_ind = grad_vec_pos < 0
  if(sum(neg_ind) == 0 | min(alpha + eps*grad_vec) >= 0) return(alpha + eps*grad_vec)
  else{
    detect_change[!neg_ind] = Inf
    detect_change[neg_ind] = -alpha[strata][neg_ind]/grad_vec_pos[neg_ind]
    min_ind = ((1:b)[strata])[which.min(detect_change)]
    eps_prime = min(detect_change)
    alpha_new = alpha + eps_prime*grad_vec
    alpha_new[min_ind] = 0
    return(alpha_new)
  }
}

opt_importance <- function(Xtra, Xtest, eps = 1e-2, tol = 1e-6, sigma = 1, N = 100){
  set.seed(12345)
  Xcenter = Xtest[sample(1:nrow(Xtest), min(100, nrow(Xtest))), ]
  g_b = g_b_fun(sigma, Xtra, Xtest, Xcenter)
  gamma = g_b$gamma
  beta = g_b$beta
  
  alpha0 = initial_point(sigma, beta, gamma, N, nrow(Xtra))
  alpha_old = alpha0 - 1
  alpha_new = alpha0
  iter = 1
  obj_new = obj_fun(alpha_new, sigma, beta)
  gamma_sq_sum = sum(gamma**2)
  eps_prime = eps
  alpha_diff = 100
  
  while(alpha_diff > tol & eps_prime > 1e-5){
    alpha_old = alpha_new
    alpha_new = update_alpha(eps_prime, alpha_old, beta, gamma)
    obj_old = obj_new
    obj_new = obj_fun(alpha_new, sigma, beta)
    
    if(iter%%200 == 1 | obj_new < obj_old){
      alpha_diff = norm(alpha_new - alpha_old, '2')
      if(obj_new < obj_old){
        # cat("bad situation....\n")
        # cat("object fun :", obj_new, "\n")
        # cat("eps prime  :", eps_prime, "\n")
        alpha_new = alpha_old
        alpha_diff = 100
        eps_prime = eps_prime/2
        obj_new = obj_old
      }
      else{
        # cat("object fun :", obj_new, "\n")
        # cat("alpha diff :", alpha_diff, "\n")
        # cat("constraint :", check_constraint(gamma, alpha_new, nrow(Xtra)), "\n")
        # cat("eps prime  :", eps_prime, "\n")
        eps_prime = min(2*eps_prime, eps)
      }
      # cat("--------------------------------------\n")
      # print(alpha_new)
    }
    iter = iter + 1
    if(iter > 10000){
      # cat("too many iterations...\n")
      break
    }
  }
  return(list(alpha = alpha_new, obj = obj_new, Xcenter = Xcenter))
}



ratio_vec <- function(Xtra, Xcenter, alpha, sigma){
  gamma_mat = matrix(nrow = nrow(Xtra), ncol = length(alpha))
  for(l in 1:length(alpha)) gamma_mat[, l] = apply(Xtra, 1, function(z) return(k_fun(z, Xcenter[l, ], sigma)))
  return(apply(gamma_mat, 1, function(z) return(sum(z*alpha))))
}


assign_balanced_values <- function(n, k) {
  set.seed(12345)
  base_count <- n %/% k
  remainder <- n %% k     
  counts <- rep(base_count, k)
  
  if (remainder > 0) counts[1:remainder] <- counts[1:remainder] + 1
  values <- unlist(mapply(rep, 1:k, counts))
  sample(values)
}

KLIEP_CV <- function(Xtra, Xtest, R, sigma_vec, eps = 100, tol = 1e-5, N = 100){
  cv_assign = assign_balanced_values(nrow(Xtest), R)
  score_vec = rep(0, length(sigma_vec))
  for(i in 1:length(sigma_vec)){
    sigma = sigma_vec[i]
    for(r in 1:R){
      cat("sigma :", sigma, ", r :", r, "\n")
      Xtest_r = Xtest[cv_assign == r, ]
      Xtest_nr = Xtest[cv_assign != r, ]
      opt = opt_importance(Xtra, Xtest_nr, eps, tol, sigma, N)
      score_vec[i] = score_vec[i] + obj_fun_cv(opt$alpha, sigma, Xtest_r, opt$Xcenter)/R
    }
  }
  sigma = sigma_vec[which.max(score_vec)]
  return(list(sigma = sigma, score_vec = score_vec))
}

# b = KLIEP_CV(X_tra, X_test, R = 5, sigma_vec = seq(from = 0.5, to = 2.5, by = 0.5), eps = 200, tol = 1e-4)
# d = opt_importance(X_tra, X_test, eps = 200, tol = 1e-5, sigma = b$sigma)
# rvec = ratio_vec(X_tra, d$Xcenter, d$alpha, b$sigma)

rvec_fun <- function(Xratio, Xtra, Xtest, R = 5, eps = 200, tol = 1e-4, sigma_vec = seq(from = 1, to = 2, by = 0.25), sigma = 1.2, cv = F){
  if(cv){
    b = KLIEP_CV(Xtra, Xtest, R, sigma_vec, eps, tol)
    d = opt_importance(Xtra, Xtest, eps, tol, b$sigma)
    rvec = ratio_vec(Xratio, d$Xcenter, d$alpha, b$sigma)
  }
  else{
    d = opt_importance(Xtra, Xtest, eps, tol, sigma)
    rvec = ratio_vec(Xratio, d$Xcenter, d$alpha, sigma)
  }
  return(rvec)
}




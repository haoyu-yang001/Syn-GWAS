score_test_ObsG_power_sim_one_g_to_one_Y <- function(mydf) {
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  GRM <- mydf$GRM
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  
  n_obs <- length(obs_protein_index)
  
  data_obsG <- data.frame(
    y = test_df$Y_obs,
    x = test_df$X
  )%>% na.omit()
  
  model_lm1 <- lm(y ~ x, data = data_obsG)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  
  a1 <- as.numeric(t(model_lm1$residuals) %*% GRM_o %*% model_lm1$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm1$residuals * model_lm1$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  
  # 计算 Σ11^{-1}
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df$Y_obs[obs_protein_index]
  
  X_all <- data.frame(intercept = rep(1,n),test_df$X)
  X_obs <- X_all[obs_protein_index,]
  
  X_obs <- as.matrix(X_obs)
  X_all <- as.matrix(X_all)
  
  G_all <- mydf$G_all
  
  SX <- inv_Sigma11 %*% X_obs
  XtSX_inv <- solve(crossprod(X_obs, SX))
  beta_hat <- XtSX_inv %*% (t(SX) %*% Y)
  residual <- Y - X_obs %*% beta_hat
  invS11_res <- inv_Sigma11 %*% residual
  
  ## for obs test statistic
  G_obs <- G_all[obs_protein_index]
  SU <- as.numeric(crossprod(G_obs, invS11_res))
  A <- crossprod(SX, G_obs)
  VU <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_ObsG","negative_log10_pval_ObsG")
  
  return(results)
}

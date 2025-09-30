SynSurrG_typeIerror_step1 <- function(mydf) {
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  
  # r <- rank(test_df %>% pull(Y_all))
  # test_df$Y_all_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of y_all
  # test_df <- INT(test_df, "Y_obs") 
  # r <- rank(test_df %>% pull(yhat))
  # test_df$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of yhat
  # 
  GRM <- mydf$GRM
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  GRM_obs_full <- GRM[obs_protein_index,1:n] # GRM子矩阵，维度 n_obs × n
  
  n_obs <- length(obs_protein_index)
  
  reml_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df$X
  )%>% na.omit()
  
  lm_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df$X
  )
  
  model_lm1 <- lm(y ~ x, data = reml_data_SynSurrG)
  model_lm2 <- lm(haty ~ x, data = lm_data_SynSurrG)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  GRM_oall <- GRM
  diag(GRM_oall) <- 0
  
  a1 <- as.numeric(t(model_lm1$residuals) %*% GRM_o %*% model_lm1$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm1$residuals * model_lm1$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  a1 <- as.numeric(t(model_lm2$residuals) %*% GRM_oall %*% model_lm2$residuals)
  hat_tau_S2 <- a1 / sum((GRM_oall)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm2$residuals * model_lm2$residuals)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM)))/nrow(GRM)
  
  a1 <- as.numeric(t(model_lm2$residuals[obs_protein_index]) %*% GRM_o %*% model_lm1$residuals)
  hat_tau_TS <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm2$residuals[obs_protein_index] * model_lm1$residuals)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  
  # obtain Σ12: n_obs × n
  I_values <- rep(hat_sigma_TS, min(n_obs, n))
  I_matrix <- sparseMatrix(
    i = 1:min(n_obs, n),j = 1:min(n_obs, n),
    x = I_values,dims = c(n_obs, n)
  )
  Sigma12 <- hat_tau_TS * GRM_obs_full + I_matrix
  
  # obtain Σ22: n × n
  Sigma22 <- hat_tau_S2 * GRM + hat_sigma_S2 * Diagonal(n)
  
  # 计算 Σ11^{-1}
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  # obtain Σ22^{-1} block wise inverse 
  inv_Sigma22 <- matrix_inv_block(wait_matrix=Sigma22)
  
  #### compute v11 v12 v22
  # Sigma11, Sigma12, Sigma22, and inv_Sigma22 (the block‐wise inverse from merged_blocks)
  # Note that Sigma21 = t(Sigma12)
  # Compute the product Σ12 * Σ22⁻¹
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  # Compute A = Σ11 - Σ12 * Σ22⁻¹ * Σ21
  # Here, Σ21 = t(Sigma12)
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  
  V11 <- matrix_inv_Amatrix(Amatrix = A)
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df$Y_obs[obs_protein_index]
  hatY <- test_df$yhat
  
  X_all <- data.frame(intercept = rep(1,n),test_df$X)
  X_obs <- X_all[obs_protein_index,]
  
  X_obs <- as.matrix(X_obs)
  X_all <- as.matrix(X_all)
  
  ## the following is the value that do not need to compute for each snp
  
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))
  
  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
  B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  
  bt2 <- Sigma12_Sigma22inv %*% X_all
  
  params <- c(
    "obs_protein_index","V11", "Y", "hatY", "X_obs", "Att",
    "A22", "Atb", "B1", "B_mat",
    "B2_1", "B2", "bt2", "Btt1",
    "inv_Sigma22", "Sigma11",
    "Sigma12_Sigma22inv", "Sigma12"
  )
  
  SynSurrG_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(SynSurrG_typeIerror_step1_pars)
}

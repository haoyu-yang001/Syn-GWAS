SynSurr_typeIerror_step1 <- function(mydf,independent_indices) {
  
  test_df_independent <- data.frame(X = mydf$X_all[independent_indices], 
                                    Y_all = mydf$Y[independent_indices], 
                                    Y_obs = mydf$Y_obs[independent_indices], 
                                    yhat = mydf$S[independent_indices])
  
  n <- nrow(test_df_independent)
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))
  
  n_obs <- length(obs_protein_index)
  
  reml_data_SynSurr <- data.frame(
    y = test_df_independent$Y_obs,
    haty = test_df_independent$yhat,
    x = test_df_independent$X
  )%>% na.omit()
  
  lm_data_SynSurr <- data.frame(
    y = test_df_independent$Y_obs,
    haty = test_df_independent$yhat,
    x = test_df_independent$X
  )
  
  model_lm1 <- lm(y ~ x, data = reml_data_SynSurr)
  model_lm2 <- lm(haty ~ x, data = lm_data_SynSurr)
  
  hat_sigma_T2 <- var(model_lm1$residuals)
  hat_sigma_S2 <- var(model_lm2$residuals)
  hat_sigma_TS <- cov(model_lm2$residuals[obs_protein_index],model_lm1$residuals)
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_sigma_T2 * Diagonal(n_obs)
  
  # obtain Σ12: n_obs × n
  I_values <- rep(hat_sigma_TS, min(n_obs, n))
  I_matrix <- sparseMatrix(
    i = 1:min(n_obs, n),j = 1:min(n_obs, n),
    x = I_values,dims = c(n_obs, n)
  )
  Sigma12 <- I_matrix
  
  # obtain Σ22: n × n
  Sigma22 <- hat_sigma_S2 * Diagonal(n)
  
  # 计算 Σ11^{-1}
  chol_Sigma11 <- chol(Sigma11)         # Cholesky
  inv_Sigma11 <- chol2inv(chol_Sigma11) 
  
  chol_Sigma22 <- chol(Sigma22)         # Cholesky
  inv_Sigma22 <- chol2inv(chol_Sigma22) 
  
  #### compute v11 v12 v22
  # Sigma11, Sigma12, Sigma22, and inv_Sigma22 (the block‐wise inverse from merged_blocks)
  # Note that Sigma21 = t(Sigma12)
  # Compute the product Σ12 * Σ22⁻¹
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  # Compute A = Σ11 - Σ12 * Σ22⁻¹ * Σ21
  # Here, Σ21 = t(Sigma12)
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  chol_A <- chol(A)         # Cholesky
  V11 <- chol2inv(chol_A) 
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df_independent$Y_obs[obs_protein_index]
  hatY <- test_df_independent$yhat
  
  X_all <- data.frame(intercept = rep(1,n),test_df_independent$X)
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
  
  SynSurr_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(SynSurr_typeIerror_step1_pars)
}

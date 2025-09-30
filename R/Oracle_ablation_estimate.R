Oracle_ablation_estimate <- function(mydf,independent_indices) {
  
  test_df_independent <- data.frame(X = mydf$X_all[independent_indices,], 
                                    Y_all = mydf$Y[independent_indices], 
                                    Y_obs = mydf$Y_obs[independent_indices], 
                                    yhat = mydf$S[independent_indices])
  
  n <- nrow(test_df_independent)
  
  data_oracle <- data.frame(
    y = test_df_independent$Y_all,
    x = test_df_independent[, grepl("^X\\.", names(test_df_independent))]
  )
  
  preds <- grep("^x\\.X\\.", names(data_oracle), value = TRUE)
  fml <- reformulate(preds, response = "y")
  
  model_lm <- lm(fml, data = data_oracle)
  
  hat_sigma_T2 <- var(model_lm$residuals)
  
  Sigma11 <- hat_sigma_T2 * Diagonal(n)
  
  chol_Sigma11 <- chol(Sigma11)         # Cholesky
  inv_Sigma11 <- chol2inv(chol_Sigma11) 
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_oracle$y
  X <- as.matrix(cbind(rep(1,n),test_df_independent[, grepl("^X\\.", names(test_df_independent))]))
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  beta_hat <- XtSX_inv %*% (t(SX) %*% Y)
  residual <- Y - X %*% beta_hat
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "invS11_res","SX","inv_Sigma11","XtSX_inv"
  )
  
  Oracle_ablation_est_pars <- mget(params, envir = environment())
  
  return(Oracle_ablation_est_pars)
  
}

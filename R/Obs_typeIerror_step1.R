Obs_typeIerror_step1 <- function(mydf,independent_indices) {
  
  test_df_independent <- data.frame(X = mydf$X_all[independent_indices], 
                                    Y_all = mydf$Y[independent_indices], 
                                    Y_obs = mydf$Y_obs[independent_indices], 
                                    yhat = mydf$S[independent_indices])
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))
  
  data_obs <- data.frame(
    y = test_df_independent$Y_obs,
    x = test_df_independent$X
  )%>% na.omit()
  
  n_obs <- nrow(data_obs)
  
  model_lm <- lm(y ~ x, data = data_obs)
  
  hat_sigma_T2 <- var(model_lm$residuals)
  Sigma11 <- hat_sigma_T2 * Diagonal(n_obs)
  
  chol_Sigma11 <- chol(Sigma11)         # Cholesky
  inv_Sigma11 <- chol2inv(chol_Sigma11) 
  
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_obs$y
  X <- cbind(rep(1,length(Y)),data_obs$x)
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  beta_hat <- XtSX_inv %*% (t(SX) %*% Y)
  residual <- Y - X %*% beta_hat
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "obs_protein_index","invS11_res","SX","inv_Sigma11","XtSX_inv"
  )
  
  Obs_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(Obs_typeIerror_step1_pars)
  
}

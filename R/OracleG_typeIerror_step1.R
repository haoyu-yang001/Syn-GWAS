OracleG_typeIerror_step1 <- function(mydf) {
  
  GRM <- mydf$GRM
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  # r <- rank(test_df %>% pull(Y_all))
  # test_df$Y_all_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of y_all
  # test_df <- INT(test_df, "Y_obs") 
  # r <- rank(test_df %>% pull(yhat))
  # test_df$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of yhat
  
  reml_data_oracle <- data.frame(
    y = test_df$Y_all,
    x = test_df$X
  )
  
  model_lm <- lm(y ~ x, data = reml_data_oracle)
  
  GRM_oall <- GRM
  diag(GRM_oall) <- 0
  
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_oall %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_oall)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM)))/nrow(GRM)
  
  Sigma11 <- hat_tau_T2 * GRM + hat_sigma_T2 * Diagonal(n)
  
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- reml_data_oracle$y
  X <- cbind(rep(1,n),reml_data_oracle$x)
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  beta_hat <- XtSX_inv %*% (t(SX) %*% Y)
  residual <- Y - X %*% beta_hat
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "invS11_res","SX","inv_Sigma11","XtSX_inv"
  )
  
  OracleG_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(OracleG_typeIerror_step1_pars)
  
}

#' @description （可选）更长说明
#' @export

DGP_GRMr_gmodel_step3 <- function(n_obs,miss,tauT = 0.7, sigmaT = 0.7, pve_x = 0, pve_g = 0.0001,
                                  maf, kinsetting, Genotype_each=NULL, pre_vars_step2, GRM) {
  
  list2env(pre_vars_step2, envir = environment())
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  # Genotype and covariate.
  if(!is.na(maf)){
    G_all <- generate_gmatrix_kinship(n_obs,miss,maf,kinsetting)
  }
  if(!is.null(Genotype_each)){
    G_all <- Genotype_each
  }
  X_all <- stats::rnorm(n = n_all)
  Z_all <- cbind(G_all,X_all)
  
  # 计算随机效应和误差项的总变异
  r <- tauT^2 + sigmaT^2
  # 计算总表型变异
  V_total <- r / (1 - pve_g - pve_x)
  
  # 根据期望的贡献比例设置固定效应系数
  beta_G <- sqrt(pve_g * V_total /var(G_all))
  beta_X <- sqrt(pve_x * V_total/var(X_all)) 
  
  # 固定效应部分
  mu_all <- Z_all %*% c(beta_G, beta_X)
  
  # 生成hatY
  S <- as.numeric(mu_all + t(chol_Sigma22) %*% rnorm(n_all))
  
  # 生成Y
  mu_cond <- mu_all[1:n_obs] + Sigma12_Sigma22inv %*% (S - mu_all)
  epsl <- rnorm(n_all)
  Y_obs <- rep(NA,n_all)
  Y_obs[1:n_obs] <- as.numeric(mu_cond + t(L_cond) %*% epsl[1:n_obs])
  
  ### generate oracle Y 
  mu_cond_oracle <- mu_all + Sigma12_oracle_Sigma22inv %*% (S - mu_all)
  Y <- as.numeric(mu_cond_oracle + t(L_cond_oracle) %*% epsl)
  
  out <- list(
    G_all = G_all,
    X_all = X_all,
    S = S,
    Y = Y,
    Y_obs = Y_obs,
    GRM = GRM
  )
  return(out)
}

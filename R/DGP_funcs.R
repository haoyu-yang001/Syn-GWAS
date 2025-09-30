DGP <- function(n_obs,miss,rho,maf = 0.25,tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5, 
                pve_g = 0, pve_x = 0.10) {
  
  ## using rho to generate the tauTS and sigmaTS
  tauTS <- sqrt((tauT^2+sigmaT^2)*(tauS^2+sigmaS^2))*rho*3/5
  sigmaTS <- sqrt((tauT^2+sigmaT^2)*(tauS^2+sigmaS^2))*rho*2/5
  
  # Subjects with observed target outcomes.
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  # Genotype and covariate.
  G_all <- stats::rbinom(n = n_all, size = 2, prob = maf)
  X_all <- stats::rnorm(n = n_all)
  Z_all <- cbind(G_all, X_all)
  
  # 计算随机效应和误差项的总变异
  r <- tauT^2 + sigmaT^2
  # 计算总表型变异
  V_total <- r / (1 - pve_g - pve_x)
  
  # 根据期望的贡献比例设置固定效应系数
  beta_G <- sqrt(pve_g * V_total /var(G_all))
  beta_X <- sqrt(pve_x * V_total/var(X_all)) 
  
  # generate the block-wise GRM matrix
  blocks <- seq(1, n_all-1, by = 2)
  rows <- rep(blocks, each = 4) + c(0, 0, 1, 1)  # 行索引
  cols <- rep(blocks, each = 4) + c(0, 1, 0, 1)  # 列索引
  values <- rep(c(1, 0.5, 0.5, 1), times = length(blocks))  # 元素值
  
  GRM <- sparseMatrix(i = rows,j = cols,x = values, dims = c(n_all, n_all))
  
  GRM_obs <- GRM[1:n_obs,1:n_obs]
  GRM_obs_all <- GRM[1:n_obs,]
  colnames(GRM_obs) <- 1:n_obs
  rownames(GRM_obs) <- 1:n_obs
  
  colnames(GRM) <- 1:n_all
  rownames(GRM) <- 1:n_all
  
  # 固定效应部分
  mu_all <- Z_all %*% c(beta_G, beta_X)
  
  ## generate observed Y_obs and suggrate S
  Sigma11 <- tauT^2 * GRM_obs + sigmaT^2 * Diagonal(n = n_obs, x = 1)
  Sigma22 <- tauS^2 * GRM + sigmaS^2 * Diagonal(n = n_all, x = 1)
  
  I_matrix <- sparseMatrix(
    i = 1:n_obs,j = 1:n_obs,x = sigmaTS,
    dims = c(n_obs, n_all)
  )
  
  Sigma12 <- tauTS * GRM_obs_all + I_matrix
  Sigma21 <- t(Sigma12)  # 对称
  
  # Sigma22 的逆
  chol_Sigma22 <- chol(Sigma22)         # Cholesky
  inv_Sigma22 <- chol2inv(chol_Sigma22) 
  
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  
  S <- as.numeric(mu_all + t(chol_Sigma22) %*% rnorm(n_all))
  
  mu_cond <- mu_all[1:n_obs] + Sigma12 %*% inv_Sigma22 %*% (S - mu_all)
  Sigma_cond <- Sigma11 - Sigma12 %*% inv_Sigma22 %*% Sigma21
  #L_cond <- chol(Sigma_cond)
  #Sigma_cond <- as.matrix(nearPD(Sigma_cond)$mat)
  L_cond <- chol(Sigma_cond)
  epsl <- rnorm(n_all)
  Y_obs <- rep(NA,n_all)
  Y_obs[1:n_obs] <- as.numeric(mu_cond + t(L_cond) %*% epsl[1:n_obs])
  
  ### generate oracle Y 
  Sigma11_oracle <- tauT^2 * GRM + sigmaT^2 * Diagonal(n = n_all, x = 1)
  I_matrix_oracle <- sparseMatrix(
    i = 1:n_all, j = 1:n_all, x = sigmaTS, dims = c(n_all, n_all)
  )
  
  Sigma12_oracle <- tauTS * GRM + I_matrix_oracle
  Sigma21_oracle <- t(Sigma12_oracle)  # 对称
  
  Sigma12_oracle_Sigma22inv <- Sigma12_oracle %*% inv_Sigma22
  
  mu_cond_oracle <- mu_all + Sigma12_oracle %*% inv_Sigma22 %*% (S - mu_all)
  Sigma_cond_oracle <- Sigma11_oracle - Sigma12_oracle %*% inv_Sigma22 %*% Sigma21_oracle
  L_cond_oracle <- chol(Sigma_cond_oracle)
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

DGP_GRMr <- function(n_obs,miss,GRMr){
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  # generate the block-wise GRM matrix
  blocks <- seq(1, n_all-1, by = 2)
  rows <- rep(blocks, each = 4) + c(0, 0, 1, 1)  # 行索引
  cols <- rep(blocks, each = 4) + c(0, 1, 0, 1)  # 列索引
  values <- rep(c(1, GRMr, GRMr, 1), times = length(blocks))  # 元素值
  
  GRM <- sparseMatrix(i = rows,j = cols,x = values, dims = c(n_all, n_all))
  return(GRM)
}

DGP_GRMr_nullmodel_step1 <- function(n_obs,miss,tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5,GRM){
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  GRM_obs <- GRM[1:n_obs,1:n_obs]
  GRM_obs_all <- GRM[1:n_obs,]
  colnames(GRM_obs) <- 1:n_obs
  rownames(GRM_obs) <- 1:n_obs
  
  colnames(GRM) <- 1:n_all
  rownames(GRM) <- 1:n_all
  
  ## 预先计算一些量，避免重复
  Sigma11 <- tauT^2 * GRM_obs + sigmaT^2 * Diagonal(n = n_obs, x = 1)
  Sigma22 <- tauS^2 * GRM + sigmaS^2 * Diagonal(n = n_all, x = 1)
  
  # Sigma22 的逆
  chol_Sigma22 <- chol(Sigma22)         # Cholesky
  inv_Sigma22 <- chol2inv(chol_Sigma22) 
  
  Sigma11_oracle <- tauT^2 * GRM + sigmaT^2 * Diagonal(n = n_all, x = 1)
  
  params <- c("Sigma11","Sigma22","inv_Sigma22","Sigma11_oracle","chol_Sigma22")
  pre_vars_step1 <- mget(params, envir = environment())
  
  return(pre_vars_step1)
}

DGP_GRMr_nullmodel_step2 <- function(n_obs,miss,tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5,
                                     rho,GRM,pre_vars_step1){
  
  list2env(pre_vars_step1, envir = environment())
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  GRM_obs <- GRM[1:n_obs,1:n_obs]
  GRM_obs_all <- GRM[1:n_obs,]
  colnames(GRM_obs) <- 1:n_obs
  rownames(GRM_obs) <- 1:n_obs
  
  colnames(GRM) <- 1:n_all
  rownames(GRM) <- 1:n_all
  
  ## using rho to generate the tauTS and sigmaTS
  tauTS <- sqrt((tauT^2+sigmaT^2)*(tauS^2+sigmaS^2))*rho*3/5
  sigmaTS <- sqrt((tauT^2+sigmaT^2)*(tauS^2+sigmaS^2))*rho*2/5
  
  I_matrix <- sparseMatrix(
    i = 1:n_obs,j = 1:n_obs,x = sigmaTS,
    dims = c(n_obs, n_all)
  )
  
  Sigma12 <- tauTS * GRM_obs_all + I_matrix
  Sigma21 <- t(Sigma12)  # 对称
  
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  
  Sigma_cond <- Sigma11 - Sigma12_Sigma22inv %*% Sigma21
  #L_cond <- chol(Sigma_cond)
  #Sigma_cond <- as.matrix(nearPD(Sigma_cond)$mat)
  L_cond <- chol(Sigma_cond)
  
  I_matrix_oracle <- sparseMatrix(
    i = 1:n_all, j = 1:n_all, x = sigmaTS, dims = c(n_all, n_all)
  )
  
  Sigma12_oracle <- tauTS * GRM + I_matrix_oracle
  Sigma21_oracle <- t(Sigma12_oracle)  # 对称
  
  Sigma12_oracle_Sigma22inv <- Sigma12_oracle %*% inv_Sigma22
  
  Sigma_cond_oracle <- Sigma11_oracle - Sigma12_oracle_Sigma22inv %*% Sigma21_oracle
  L_cond_oracle <- chol(Sigma_cond_oracle)
  
  params <- c("chol_Sigma22","Sigma12_Sigma22inv","L_cond",
              "Sigma12_oracle_Sigma22inv","L_cond_oracle")
  pre_vars_step2 <- mget(params, envir = environment())
  
  return(pre_vars_step2)
}

DGP_GRMr_nullmodel_step3 <- function(n_obs,miss,tauT = 0.7, sigmaT = 0.7, pve_x = 0.10, 
                                     pre_vars_step2, GRM) {
  
  list2env(pre_vars_step2, envir = environment())
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  # Genotype 0 and covariate.
  X_all <- stats::rnorm(n = n_all)
  Z_all <- cbind(rep(0,n_all),X_all)
  
  # 计算随机效应和误差项的总变异
  r <- tauT^2 + sigmaT^2
  # 计算总表型变异
  V_total <- r / (1 - 0 - pve_x)
  
  # 根据期望的贡献比例设置固定效应系数
  beta_G <- 0
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
    X_all = X_all,
    S = S,
    Y = Y,
    Y_obs = Y_obs,
    GRM = GRM
  )
  return(out)
}

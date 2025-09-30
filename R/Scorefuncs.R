compute_score <- function(Ug, Vu) {
  score <- as.numeric(Ug^2 / Vu)
  negative_log10_pval <- -pchisq(score, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  list(score = score, negative_log10_pval = negative_log10_pval)
}

score_test_SynSurrG_multiply <- function(g_matrix,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  
  n_obs <- length(obs_protein_index)
  n <- nrow(g_matrix)
  
  snpindex <- 1:ncol(g_matrix)
  
  G_all <- g_matrix #[,snpindex]
  G_obs <- as.matrix(G_all[obs_protein_index,])
  
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  #Att <- inv_Sigma22 %*% X_all
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X
  #A22 <- crossprod(X_all, Att)
  
  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
  
  bt1 <- Sigma12_Sigma22inv %*% G_all
  
  Ah <- t(temA) %*% t(Sigma12)
  
  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)
  
  V1 <- rowSums((VV1 %*% Sigma11) * VV1)
  
  # 定义单个 SNP 计算的函数
  compute_snp <- function(snp_d) {
    
    ## 1. estimate alpha 
    # 构造 A 部分
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)
    # 构造 alpha 右端项 
    Ar0 <- rbind(t(temA[,snp_d]),t(Att)) # Z %*% Sigma_22^{-1}
    # 得到 alpha 估计值
    alpha <- solve(A_mat, Ar0 %*% hatY)
    
    
    ## 2. estimate beta
    #Btt1 <- V11 %*% X_obs
    #Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
    #B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    #B1 <- B_mat %*% Y
    
    
    #bt2 <- Sigma12_Sigma22inv %*% X_all
    b1 <- cbind(bt1[,snp_d],bt2)
    #b3 <- Sigma12_Sigma22inv %*% hatY
    
    #Sigma12_Sigma22inv_rs_alpha <- b3 - b1 %*% alpha # Sigma12_Sigma22inv (hatY-Z alpha)
    #beta <- B1 - B_mat %*% Sigma12_Sigma22inv_rs_alpha
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)
    
    ## 3. calculate the score
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))
    
    ## 4. calculate the variance of the score
    
    #the first element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV1 <- t(Stt1[snp_d,] - VV1_2[snp_d,]) # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    
    
    #the second element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV2
    #VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
    Arr <- rbind(Ah[snp_d,],Atb)
    Ar <- solve(A_mat, Arr)
    
    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar
    
    V2 <- VV2 %*% VV1[snp_d,]
    
    #V1 <- VV1 %*% Sigma11  %*% t(VV1)
    #V2 <- VV2 %*% t(Sigma12)  %*% t(VV1)
    #V3 <- VV2 %*% Sigma22  %*% t(VV2)
    
    VU <- as.numeric(V1[snp_d] + V2)
    
    inv_beta <- Stt1[snp_d,] %*% G_obs[,snp_d]
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)
    
    list(#alpha = alpha, beta = beta, 
      S_syn=S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }
  
  # 利用 lapply 对所有 SNP 进行计算
  results <- lapply(seq_along(snpindex), compute_snp)
  
  # 将各个结果整合成矩阵和向量
  #hat_alpha   <- do.call(cbind, lapply(results, `[[`, "alpha"))
  #hat_beta    <- do.call(cbind, lapply(results, `[[`, "beta"))
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")
  
  #S_obs <- as.numeric(crossprod(G_obs, invS11_res))
  #A <- crossprod(SX, G_obs)
  #VU_obs <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * XtSX_inv %*% A))
  
  surr_result <- compute_score(S_syn, VU_syn)
  #obs_result <- compute_score(S_obs, VU_obs)
  
  list(
    T_score_SynSurrG = surr_result$score,
    negative_log10_pval_SynSurrG = surr_result$negative_log10_pval,
    hat_beta_SynSurrG = hat_beta_syn,
    var_hat_beta_SynSurrG = var_hat_beta_syn
  )
}

score_test_SynSurr_multiply <- function(g_matrix,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  g_matrix <- as.matrix(g_matrix[independent_indices,])
  
  n_obs <- length(obs_protein_index)
  n <- nrow(g_matrix)
  
  snpindex <- 1:ncol(g_matrix)
  
  G_all <- g_matrix #[,snpindex]
  G_obs <- as.matrix(G_all[obs_protein_index,])
  
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  #Att <- inv_Sigma22 %*% X_all
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X
  #A22 <- crossprod(X_all, Att)
  
  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
  
  bt1 <- Sigma12_Sigma22inv %*% G_all
  
  Ah <- t(temA) %*% t(Sigma12)
  
  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)
  
  V1 <- rowSums((VV1 %*% Sigma11) * VV1)
  
  # 定义单个 SNP 计算的函数
  compute_snp <- function(snp_d) {
    
    ## 1. estimate alpha 
    # 构造 A 部分
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)
    # 构造 alpha 右端项 
    Ar0 <- rbind(t(temA[,snp_d]),t(Att)) # Z %*% Sigma_22^{-1}
    # 得到 alpha 估计值
    alpha <- solve(A_mat, Ar0 %*% hatY)
    
    
    ## 2. estimate beta
    #Btt1 <- V11 %*% X_obs
    #Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
    #B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    #B1 <- B_mat %*% Y
    
    
    #bt2 <- Sigma12_Sigma22inv %*% X_all
    b1 <- cbind(bt1[,snp_d],bt2)
    #b3 <- Sigma12_Sigma22inv %*% hatY
    
    #Sigma12_Sigma22inv_rs_alpha <- b3 - b1 %*% alpha # Sigma12_Sigma22inv (hatY-Z alpha)
    #beta <- B1 - B_mat %*% Sigma12_Sigma22inv_rs_alpha
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)
    
    ## 3. calculate the score
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))
    
    ## 4. calculate the variance of the score
    
    #the first element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV1 <- t(Stt1[snp_d,] - VV1_2[snp_d,]) # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    
    
    #the second element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV2
    #VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
    Arr <- rbind(Ah[snp_d,],Atb)
    Ar <- solve(A_mat, Arr)
    
    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar
    
    V2 <- VV2 %*% VV1[snp_d,]
    
    #V1 <- VV1 %*% Sigma11  %*% t(VV1)
    #V2 <- VV2 %*% t(Sigma12)  %*% t(VV1)
    #V3 <- VV2 %*% Sigma22  %*% t(VV2)
    
    VU <- as.numeric(V1[snp_d] + V2)
    
    inv_beta <- Stt1[snp_d,] %*% G_obs[,snp_d]
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)
    
    list(#alpha = alpha, beta = beta, 
      S_syn=S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }
  
  # 利用 lapply 对所有 SNP 进行计算
  results <- lapply(seq_along(snpindex), compute_snp)
  
  # 将各个结果整合成矩阵和向量
  #hat_alpha   <- do.call(cbind, lapply(results, `[[`, "alpha"))
  #hat_beta    <- do.call(cbind, lapply(results, `[[`, "beta"))
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")
  
  #S_obs <- as.numeric(crossprod(G_obs, invS11_res))
  #A <- crossprod(SX, G_obs)
  #VU_obs <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * XtSX_inv %*% A))
  
  surr_result <- compute_score(S_syn, VU_syn)
  #obs_result <- compute_score(S_obs, VU_obs)
  
  list(
    T_score_SynSurr = surr_result$score,
    negative_log10_pval_SynSurr = surr_result$negative_log10_pval,
    hat_beta_SynSurr = hat_beta_syn,
    var_hat_beta_SynSurr = var_hat_beta_syn
  )
}

score_test_OracleG_multiply <- function(g_matrix,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  
  ## for oracle test statistic
  SU <- as.numeric(crossprod(g_matrix, invS11_res))
  A <- crossprod(SX, g_matrix)
  VU <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_OracleG","negative_log10_pval_OracleG")
  
  inv_beta <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_OracleG <- hat_beta
  results$var_hat_beta_OracleG <- var_hat_beta
  
  return(results)
}

score_test_Oracle_multiply <- function(g_matrix,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  g_matrix <- as.matrix(g_matrix[independent_indices,])
  
  ## for oracle test statistic
  SU <- as.numeric(crossprod(g_matrix, invS11_res))
  A <- crossprod(SX, g_matrix)
  VU <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_Oracle","negative_log10_pval_Oracle")
  
  inv_beta <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_Oracle <- hat_beta
  results$var_hat_beta_Oracle <- var_hat_beta
  
  return(results)
}

score_test_ObsG_multiply <- function(g_matrix,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  g_matrix_sub <- g_matrix[obs_protein_index,]
  
  ## for obs test statistic
  SU <- as.numeric(crossprod(g_matrix_sub, invS11_res))
  A <- crossprod(SX, g_matrix_sub)
  VU <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_ObsG","negative_log10_pval_ObsG")
  
  inv_beta <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_ObsG <- hat_beta
  results$var_hat_beta_ObsG <- var_hat_beta
  
  return(results)
}

score_test_Obs_multiply <- function(g_matrix,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  g_matrix_all_sub <- as.matrix(g_matrix[independent_indices,])
  
  g_matrix_sub <- g_matrix_all_sub[obs_protein_index,]
  
  ## for obs test statistic
  SU <- as.numeric(crossprod(g_matrix_sub, invS11_res))
  A <- crossprod(SX, g_matrix_sub)
  VU <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_Obs","negative_log10_pval_Obs")
  
  inv_beta <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_Obs <- hat_beta
  results$var_hat_beta_Obs <- var_hat_beta
  
  return(results)
}

score_test_SynSurrG_single <- function(G_all,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  
  n_obs <- length(obs_protein_index)
  n <- length(G_all)
  
  G_obs <- G_all[obs_protein_index]
  
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  #Att <- inv_Sigma22 %*% X_all
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X
  #A22 <- crossprod(X_all, Att)
  
  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
  
  bt1 <- Sigma12_Sigma22inv %*% G_all
  
  Ah <- t(temA) %*% t(Sigma12)
  
  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)
  
  V1 <- rowSums((VV1 %*% Sigma11) * VV1)
  
  # 定义单个 SNP 计算的函数
  compute_snp <- function(snp_d) {
    
    ## 1. estimate alpha 
    # 构造 A 部分
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)
    # 构造 alpha 右端项 
    Ar0 <- rbind(t(temA[,snp_d]),t(Att)) # Z %*% Sigma_22^{-1}
    # 得到 alpha 估计值
    alpha <- solve(A_mat, Ar0 %*% hatY)
    
    
    ## 2. estimate beta
    #Btt1 <- V11 %*% X_obs
    #Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
    #B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    #B1 <- B_mat %*% Y
    
    
    #bt2 <- Sigma12_Sigma22inv %*% X_all
    b1 <- cbind(bt1[,snp_d],bt2)
    #b3 <- Sigma12_Sigma22inv %*% hatY
    
    #Sigma12_Sigma22inv_rs_alpha <- b3 - b1 %*% alpha # Sigma12_Sigma22inv (hatY-Z alpha)
    #beta <- B1 - B_mat %*% Sigma12_Sigma22inv_rs_alpha
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)
    
    ## 3. calculate the score
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))
    
    ## 4. calculate the variance of the score
    
    #the first element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV1 <- t(Stt1[snp_d,] - VV1_2[snp_d,]) # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    
    
    #the second element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV2
    #VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
    Arr <- rbind(Ah[snp_d,],Atb)
    Ar <- solve(A_mat, Arr)
    
    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar
    
    V2 <- VV2 %*% VV1[snp_d,]
    
    #V1 <- VV1 %*% Sigma11  %*% t(VV1)
    #V2 <- VV2 %*% t(Sigma12)  %*% t(VV1)
    #V3 <- VV2 %*% Sigma22  %*% t(VV2)
    
    VU <- as.numeric(V1[snp_d] + V2)
    
    inv_beta <- Stt1 %*% G_obs
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)
    
    list(#alpha = alpha, beta = beta, 
      S_syn=S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }
  
  # 利用 lapply 对所有 SNP 进行计算
  results <- lapply(seq_along(1), compute_snp)
  
  # 将各个结果整合成矩阵和向量
  #hat_alpha   <- do.call(cbind, lapply(results, `[[`, "alpha"))
  #hat_beta    <- do.call(cbind, lapply(results, `[[`, "beta"))
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")
  
  #S_obs <- as.numeric(crossprod(G_obs, invS11_res))
  #A <- crossprod(SX, G_obs)
  #VU_obs <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * XtSX_inv %*% A))
  
  surr_result <- compute_score(S_syn, VU_syn)
  #obs_result <- compute_score(S_obs, VU_obs)
  
  list(
    T_score_SynSurrG = surr_result$score,
    negative_log10_pval_SynSurrG = surr_result$negative_log10_pval,
    hat_beta_SynSurrG = hat_beta_syn,
    var_hat_beta_SynSurrG = var_hat_beta_syn
  )
}

score_test_SynSurr_single <- function(G_all,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  G_all <- G_all[independent_indices]
  
  n_obs <- length(obs_protein_index)
  n <- length(G_all)
  
  G_obs <- G_all[obs_protein_index]
  
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  #Att <- inv_Sigma22 %*% X_all
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X
  #A22 <- crossprod(X_all, Att)
  
  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
  
  bt1 <- Sigma12_Sigma22inv %*% G_all
  
  Ah <- t(temA) %*% t(Sigma12)
  
  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)
  
  V1 <- rowSums((VV1 %*% Sigma11) * VV1)
  
  # 定义单个 SNP 计算的函数
  compute_snp <- function(snp_d) {
    
    ## 1. estimate alpha 
    # 构造 A 部分
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)
    # 构造 alpha 右端项 
    Ar0 <- rbind(t(temA[,snp_d]),t(Att)) # Z %*% Sigma_22^{-1}
    # 得到 alpha 估计值
    alpha <- solve(A_mat, Ar0 %*% hatY)
    
    
    ## 2. estimate beta
    #Btt1 <- V11 %*% X_obs
    #Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
    #B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    #B1 <- B_mat %*% Y
    
    
    #bt2 <- Sigma12_Sigma22inv %*% X_all
    b1 <- cbind(bt1[,snp_d],bt2)
    #b3 <- Sigma12_Sigma22inv %*% hatY
    
    #Sigma12_Sigma22inv_rs_alpha <- b3 - b1 %*% alpha # Sigma12_Sigma22inv (hatY-Z alpha)
    #beta <- B1 - B_mat %*% Sigma12_Sigma22inv_rs_alpha
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)
    
    ## 3. calculate the score
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))
    
    ## 4. calculate the variance of the score
    
    #the first element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV1 <- t(Stt1[snp_d,] - VV1_2[snp_d,]) # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    
    
    #the second element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV2
    #VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
    Arr <- rbind(Ah[snp_d,],Atb)
    Ar <- solve(A_mat, Arr)
    
    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar
    
    V2 <- VV2 %*% VV1[snp_d,]
    
    #V1 <- VV1 %*% Sigma11  %*% t(VV1)
    #V2 <- VV2 %*% t(Sigma12)  %*% t(VV1)
    #V3 <- VV2 %*% Sigma22  %*% t(VV2)
    
    VU <- as.numeric(V1[snp_d] + V2)
    
    inv_beta <- Stt1 %*% G_obs
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)
    
    list(#alpha = alpha, beta = beta, 
      S_syn=S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }
  
  # 利用 lapply 对所有 SNP 进行计算
  results <- lapply(seq_along(1), compute_snp)
  
  # 将各个结果整合成矩阵和向量
  #hat_alpha   <- do.call(cbind, lapply(results, `[[`, "alpha"))
  #hat_beta    <- do.call(cbind, lapply(results, `[[`, "beta"))
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")
  
  #S_obs <- as.numeric(crossprod(G_obs, invS11_res))
  #A <- crossprod(SX, G_obs)
  #VU_obs <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * XtSX_inv %*% A))
  
  surr_result <- compute_score(S_syn, VU_syn)
  #obs_result <- compute_score(S_obs, VU_obs)
  
  list(
    T_score_SynSurr = surr_result$score,
    negative_log10_pval_SynSurr = surr_result$negative_log10_pval,
    hat_beta_SynSurr = hat_beta_syn,
    var_hat_beta_SynSurr = var_hat_beta_syn
  )
}

score_test_OracleG_single <- function(G_all,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  
  ## for oracle test statistic
  SU <- as.numeric(crossprod(G_all, invS11_res))
  A <- crossprod(SX, G_all)
  VU <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_OracleG","negative_log10_pval_OracleG")
  
  inv_beta <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_OracleG <- hat_beta
  results$var_hat_beta_OracleG <- var_hat_beta

  return(results)
}

score_test_Oracle_single <- function(G_all,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  G_all <- G_all[independent_indices]
  
  ## for oracle test statistic
  SU <- as.numeric(crossprod(G_all, invS11_res))
  A <- crossprod(SX, G_all)
  VU <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_Oracle","negative_log10_pval_Oracle")
  
  inv_beta <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_Oracle <- hat_beta
  results$var_hat_beta_Oracle <- var_hat_beta
  
  return(results)
}

score_test_ObsG_single <- function(G_all,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  G_all_sub <- G_all[obs_protein_index]
  
  ## for obs test statistic
  SU <- as.numeric(crossprod(G_all_sub, invS11_res))
  A <- crossprod(SX, G_all_sub)
  VU <- as.numeric(colSums(G_all_sub * (inv_Sigma11 %*% G_all_sub)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_ObsG","negative_log10_pval_ObsG")
  
  inv_beta <- as.numeric(colSums(G_all_sub * (inv_Sigma11 %*% G_all_sub)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_ObsG <- hat_beta
  results$var_hat_beta_ObsG <- var_hat_beta
  
  return(results)
}

score_test_Obs_single <- function(G_all,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  G_all_sub <- G_all[independent_indices]
  G_sub <- G_all_sub[obs_protein_index]
  
  ## for obs test statistic
  SU <- as.numeric(crossprod(G_sub, invS11_res))
  A <- crossprod(SX, G_sub)
  VU <- as.numeric(colSums(G_sub * (inv_Sigma11 %*% G_sub)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_Obs","negative_log10_pval_Obs")
  
  inv_beta <- as.numeric(colSums(G_sub * (inv_Sigma11 %*% G_sub)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_Obs <- hat_beta
  results$var_hat_beta_Obs <- var_hat_beta
  
  return(results)
}

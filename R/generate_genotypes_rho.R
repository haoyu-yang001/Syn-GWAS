generate_genotypes_rho <- function(n_obs, miss, n_snps, maf = 0.25, rho = 0.5) {
  # n_families: 家庭数量(5000个家庭=10000个人)
  # n_snps: SNP数量
  # maf: 次要等位基因频率
  # rho: 家庭成员间的基因型相关性
  
  # 1. 生成独立的标准正态分布变量(家庭间独立)
  # 每个SNP生成2n_families个独立样本
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  n_families <- n_all/2
  Z <- matrix(rnorm(n_snps * 2 * n_families), nrow = n_snps)
  
  # 2. 创建家庭成员间的相关性
  # 对于每个家庭中的两个人，他们的基因型是Z1和rho*Z1 + sqrt(1-rho^2)*Z2
  # 其中Z1和Z2是独立的
  Z_family <- matrix(0, nrow = n_snps, ncol = 2 * n_families)
  
  for (i in 1:n_families) {
    # 家庭中第一个人的基因型(完全独立)
    Z_family[, 2*i - 1] <- Z[, 2*i - 1]
    
    # 家庭中第二个人的基因型(与第一个人相关)
    Z_family[, 2*i] <- rho * Z[, 2*i - 1] + sqrt(1 - rho^2) * Z[, 2*i]
  }
  
  # 3. 转换为二项分布以控制MAF
  # 计算分位数阈值
  threshold <- qnorm(1 - maf)
  
  # 初始化基因型矩阵(0,1,2)
  genotypes <- matrix(0, nrow = 2 * n_families, ncol = n_snps)
  
  # 转换为基因型
  for (i in 1:n_snps) {
    # 对于每个SNP，将正态分布转换为基因型
    # 0: 低于阈值(主要等位基因纯合子)
    # 1: 在阈值和对称阈值之间(杂合子)
    # 2: 高于对称阈值(次要等位基因纯合子)
    
    # 使用对称的双侧阈值
    upper_threshold <- threshold
    lower_threshold <- -threshold
    
    # 转换为基因型
    gt <- ifelse(Z_family[i,] < lower_threshold, 0, 
                 ifelse(Z_family[i,] > upper_threshold, 2, 1))
    
    genotypes[, i] <- gt
  }
  
  # 4. 标准化基因型(均值为0，方差为1)
  # 首先中心化
  genotypes_centered <- scale(genotypes, center = TRUE, scale = FALSE)
  
  # 然后标准化方差
  genotypes_scaled <- genotypes_centered / sqrt(2 * maf * (1 - maf))
  
  # 转置矩阵，使行为个体，列为SNP
  genotypes_final <- t(genotypes_scaled)
  
  return(t(genotypes_final))
}

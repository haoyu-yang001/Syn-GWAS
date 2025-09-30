find_blocks_vectorized <- function(A) {
  n <- nrow(A)
  
  sm <- summary(A)

  max_j <- rep(0, n)
  tmp <- tapply(sm$j, sm$i, max)
  max_j[as.integer(names(tmp))] <- tmp
  
  zero_rows <- which(max_j == 0)
  if(length(zero_rows) > 0){
    max_j[zero_rows] <- zero_rows
  }
  
  cum_max <- cummax(max_j)
  
  block_boundaries <- which(seq_len(n) == cum_max)
  
  blocks <- vector("list", length(block_boundaries))
  start <- 1
  for (i in seq_along(block_boundaries)) {
    end <- block_boundaries[i]
    blocks[[i]] <- start:end
    start <- end + 1
  }
  
  return(blocks)
}

find_blocks_vectorized_threshold <- function(A, threshold = 0.005) {
  n <- nrow(A)
  # 提取非零元素的三元组信息，并仅保留绝对值大于等于阈值的元素
  sm <- summary(A)
  sm <- sm[abs(sm$x) >= threshold, ]
  
  # 计算每一行的最大非零列索引（若某行全零，则为0）
  max_j <- rep(0, n)
  if(nrow(sm) > 0) {
    tmp <- tapply(sm$j, sm$i, max)
    max_j[as.integer(names(tmp))] <- tmp
  }
  
  # 对于全零行，将其设置为行号（表示该行仅与自身有关，不影响后续行）
  zero_rows <- which(max_j == 0)
  if (length(zero_rows) > 0) {
    max_j[zero_rows] <- zero_rows
  }
  
  # 计算每一行的“最远影响范围”
  cum_max <- cummax(max_j)
  
  # 找到所有满足行号等于累计最大值的位置，作为 block 的边界
  block_boundaries <- which(seq_len(n) == cum_max)
  
  # 根据 block 边界划分 block
  blocks <- vector("list", length(block_boundaries))
  start <- 1
  for (i in seq_along(block_boundaries)) {
    end <- block_boundaries[i]
    blocks[[i]] <- start:end
    start <- end + 1
  }
  
  return(blocks)
}

block_matrix <- function(tl, tr, bl, br) {
  rbind(cbind(tl, tr), cbind(bl, br))
}

matrix_inv_block <- function(wait_matrix,threshold = 1000){
  blocks <- find_blocks_vectorized(wait_matrix)
  merged_blocks <- list()
  current_block <- blocks[[1]]
  for (i in 2:length(blocks)) {
    if (length(current_block) + length(blocks[[i]]) < threshold) {
      current_block <- c(current_block, blocks[[i]])
    } else {
      # Otherwise, add the current block to the merged list and start a new one.
      merged_blocks <- c(merged_blocks, list(current_block))
      current_block <- blocks[[i]]
    }
  }
  # Append the last accumulated block
  merged_blocks <- c(merged_blocks, list(current_block))
  
  ncores <- 2
  block_inv_list <- mclapply(merged_blocks, function(block) {
    subSigma <- wait_matrix[block, block]
    chol_subSigma <- chol(subSigma)
    chol2inv(chol_subSigma)
  }, mc.cores = ncores)
  
  inv_wait_matrix <- bdiag(block_inv_list)
  #cat("Completed block-wise inversion Sigma11 using parallel computation.\n")
  return(inv_wait_matrix)
}

matrix_inv_Amatrix <- function(Amatrix,thr=0.006,threshold=1000){
  blocks_A <- find_blocks_vectorized_threshold(Amatrix,threshold = thr) #set the thr for block detaction
  ja <- max(sapply(blocks_A, length))
  while(ja > 1000){
    thr <- thr + 0.001
    blocks_A <- find_blocks_vectorized_threshold(Amatrix,threshold = thr) #set the thr for block detaction
    ja <- max(sapply(blocks_A, length))
  }
  
  merged_blocks <- list()
  current_block <- blocks_A[[1]]
  if(length(blocks_A) > 1 ){
    for (i in 2:length(blocks_A)) {
      if (length(current_block) + length(blocks_A[[i]]) < threshold) {
        current_block <- c(current_block, blocks_A[[i]])
      } else {
        # Otherwise, add the current block to the merged list and start a new one.
        merged_blocks <- c(merged_blocks, list(current_block))
        current_block <- blocks_A[[i]]
      }
    }
  }
  # Append the last accumulated block
  merged_blocks <- c(merged_blocks, list(current_block))
  #sapply(merged_blocks, length)
  
  ncores <- 2
  block_inv_list <- mclapply(merged_blocks, function(block) {
    subA <- Amatrix[block, block]
    #chol_subA <- chol(subA + diag(1e-5, nrow(subA)))
    #chol2inv(chol_subA)
    solve(subA)
  }, mc.cores = ncores)
  
  V11 <- bdiag(block_inv_list)
  #cat("Completed V11 block-wise inversion using parallel computation.\n")
  return(V11)
}

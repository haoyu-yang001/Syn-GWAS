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

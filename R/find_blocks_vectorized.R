find_blocks_vectorized <- function(A) {
  n <- nrow(A)
  # 提取非零元素的三元组信息
  sm <- summary(A)
  
  # 计算每一行的最大非零列索引（若某行全零，则为0）
  max_j <- rep(0, n)
  tmp <- tapply(sm$j, sm$i, max)
  max_j[as.integer(names(tmp))] <- tmp
  
  # 对于全零行，将其设置为行号（表示该行只影响自身）
  zero_rows <- which(max_j == 0)
  if(length(zero_rows) > 0){
    max_j[zero_rows] <- zero_rows
  }
  
  # 使用 cummax() 计算每一行的“最远影响范围”
  cum_max <- cummax(max_j)
  
  # 找到所有满足行号等于 cummax 的位置，即 block 的边界
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

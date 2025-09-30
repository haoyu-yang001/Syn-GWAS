find_blocks_vectorized <- function(A) {
  n <- nrow(A)

  # 统一拿 (i, j, x) 三元组
  if (inherits(A, "sparseMatrix")) {
    sm <- Matrix::summary(A)   # 确保用到 Matrix 的方法
  } else {
    nz <- which(A != 0, arr.ind = TRUE)
    sm <- data.frame(i = nz[,1], j = nz[,2], x = A[nz], row.names = NULL)
  }

  max_j <- integer(n)
  if (nrow(sm) > 0) {
    tmp <- tapply(sm$j, sm$i, max)
    max_j[as.integer(names(tmp))] <- tmp
  }
  zero_rows <- which(max_j == 0)
  if (length(zero_rows) > 0) max_j[zero_rows] <- zero_rows

  cum_max <- cummax(max_j)
  bnd <- which(seq_len(n) == cum_max)

  blocks <- vector("list", length(bnd))
  start <- 1L
  for (k in seq_along(bnd)) {
    end <- bnd[k]; blocks[[k]] <- start:end; start <- end + 1L
  }
  blocks
}

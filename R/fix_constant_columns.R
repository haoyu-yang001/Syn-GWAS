fix_constant_columns <- function(g_matrix, obs_idx) {
  g_obs <- as.matrix(g_matrix[obs_idx, ])
  ref_row <- g_obs[1, , drop = FALSE]
  is_const <- colSums(g_obs == matrix(rep(ref_row, each = nrow(g_obs)), nrow = nrow(g_obs))) == nrow(g_obs)
  
  if (any(is_const)) {
    const_cols <- which(is_const)
    sampled_rows <- replicate(length(const_cols), sample(obs_idx, 6), simplify = FALSE)
    indices <- cbind(rows = unlist(sampled_rows), cols = rep(const_cols, each = 6))
    g_matrix[indices] <- rep(c(2, 1, 1, 1, 1, 2), length(const_cols))
  }
  g_matrix
}

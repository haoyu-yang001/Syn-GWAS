typeIerror_generate_gmatrix <- function(n_obs,miss,chunk_size,maf){
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  g_matrix <- matrix(rbinom(n_all * chunk_size, size = 2, prob = maf), nrow = n_all, ncol = chunk_size)
  return(g_matrix)
}

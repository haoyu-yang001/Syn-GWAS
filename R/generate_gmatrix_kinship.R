generate_gmatrix_kinship <- function(n_obs,miss,maf,kinsetting){
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  draw_gamete <- function(g) rbinom(length(g), size = 1, prob = g / 2)
  
  g_grand_father <- rbinom(n_all/2, size = 2, prob = maf)
  g_grand_mother <- rbinom(n_all/2, size = 2, prob = maf)
  
  if(kinsetting == "unkin"){
    geno_vec <- c(rbind(g_grand_father, g_grand_mother))
  }
  if(kinsetting == "siblings"){
    g_father <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    g_uncle <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    geno_vec <- c(rbind(g_father, g_uncle))
  }
  if(kinsetting == "cousin"){
    g_father <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    g_uncle <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    g_mother <- rbinom(n_all/2, size = 2, prob = maf)
    g_me <- draw_gamete(g_father) + draw_gamete(g_mother)
    g_rand <- rbinom(n_all/2, size = 2, prob = maf)
    g_far_sister <- draw_gamete(g_uncle) + draw_gamete(g_rand)
    geno_vec <- c(rbind(g_me, g_far_sister))
  }
  if(kinsetting == "halfsiblings"){
    g_uncle1 <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    g_grand_mother2 <- rbinom(n_all/2, size = 2, prob = maf)
    g_uncle2 <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother2)
    geno_vec <- c(rbind(g_uncle1, g_uncle2))
  }
  return(geno_vec)
}

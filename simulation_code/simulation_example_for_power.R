library(devtools)
install_github("haoyu-yang001/Syn-GWAS",force = TRUE)
library(withr)
library(Matrix)
library(SynGWAS)

arraynum <- 200
num_snps_each_array <- 50

GRMr <- 0.125 #0.125 0.25 0.5
chunk_array <- 1
kinsetting <- "halfsiblings" #"cousin" "halfsiblings" "siblings"

n_obs <- 1000

snp_heritability_vec <- seq(0.0001, 0.001, 0.0001)

tauT = 0.7
tauS = 0.4
sigmaT = 0.7
sigmaS = 0.5

pve_x = 0
maf <- 0.25

miss <- 0.5


GRM <- DGP_GRMr(n_obs,miss,GRMr)

hh <- find_blocks_vectorized(GRM)
independent_indices <- with_seed(25324, {
  sapply(hh, function(x) x[sample.int(length(x), 1)])
})

pre_vars_step1 <- DGP_GRMr_nullmodel_step1(n_obs,miss,tauT, tauS, sigmaT, sigmaS, GRM)

rho <- 0.75


pre_vars_step2 <- DGP_GRMr_nullmodel_step2(n_obs,miss,tauT, tauS, sigmaT, sigmaS,rho,GRM,pre_vars_step1)

pve_g <- 0.0002

cat("Processing miss", miss, "rho", rho, "pve_g", pve_g, "chunk array:", chunk_array,"GRMr", GRMr, "\n")

results <- c()

for (ss in 1:num_snps_each_array) {
  
  mydf <- DGP_GRMr_gmodel_step3(n_obs,miss,tauT = tauT, sigmaT = sigmaT, pve_x = pve_x, pve_g = pve_g,
                                maf, kinsetting, Genotype_each=NULL, pre_vars_step2, GRM)
  
  G_all <- mydf$G_all
  
  SynSurrG_typeIerror_step1_pars <- SynSurrG_typeIerror_step1(mydf)
  OracleG_typeIerror_step1_pars <- OracleG_typeIerror_step1(mydf)
  ObsG_typeIerror_step1_pars <- ObsG_typeIerror_step1(mydf)
  
  SynSurr_typeIerror_step1_pars <- SynSurr_typeIerror_step1(mydf,independent_indices)
  Oracle_typeIerror_step1_pars <- Oracle_typeIerror_step1(mydf,independent_indices)
  Obs_typeIerror_step1_pars <- Obs_typeIerror_step1(mydf,independent_indices)
  
  SynSurrG_results <- score_test_SynSurrG_single(G_all,step1_pars = SynSurrG_typeIerror_step1_pars)
  OracleG_results <- score_test_OracleG_single(G_all,step1_pars = OracleG_typeIerror_step1_pars)
  ObsG_results <- score_test_ObsG_single(G_all,step1_pars = ObsG_typeIerror_step1_pars)
  
  SynSurr_results <- score_test_SynSurr_single(G_all,step1_pars = SynSurr_typeIerror_step1_pars,independent_indices)
  Oracle_results <- score_test_Oracle_single(G_all,step1_pars = Oracle_typeIerror_step1_pars,independent_indices)
  Obs_results <- score_test_Obs_single(G_all,step1_pars = Obs_typeIerror_step1_pars,independent_indices)
  
  results <- rbind(results,cbind(OracleG_results$negative_log10_pval_OracleG,
                                 OracleG_results$hat_beta_OracleG,
                                 OracleG_results$var_hat_beta_OracleG,
                                 ObsG_results$negative_log10_pval_ObsG,
                                 ObsG_results$hat_beta_ObsG,
                                 ObsG_results$var_hat_beta_ObsG,
                                 SynSurrG_results$negative_log10_pval_SynSurrG,
                                 SynSurrG_results$hat_beta_SynSurrG,
                                 SynSurrG_results$var_hat_beta_SynSurrG,
                                 Oracle_results$negative_log10_pval_Oracle,
                                 Oracle_results$hat_beta_Oracle,
                                 Oracle_results$var_hat_beta_Oracle,
                                 Obs_results$negative_log10_pval_Obs,
                                 Obs_results$hat_beta_Obs,
                                 Obs_results$var_hat_beta_Obs,
                                 SynSurr_results$negative_log10_pval_SynSurr,
                                 SynSurr_results$hat_beta_SynSurr,
                                 SynSurr_results$var_hat_beta_SynSurr))
  
  if (ss %% 2 == 0) print(paste0("computing to ss = ",ss))
}

colnames(results) <- c("OracleG_negative_log10_pval","OracleG_hat_beta","OracleG_var_hat_beta",
                       "ObsG_negative_log10_pval","ObsG_hat_beta","ObsG_var_hat_beta",
                       "SynSurrG_negative_log10_pval","SynSurrG_hat_beta","SynSurrG_var_hat_beta",
                       "Oracle_negative_log10_pval","Oracle_hat_beta","Oracle_var_hat_beta",
                       "Obs_negative_log10_pval","Obs_hat_beta","Obs_var_hat_beta",
                       "SynSurr_negative_log10_pval","SynSurr_hat_beta","SynSurr_var_hat_beta")

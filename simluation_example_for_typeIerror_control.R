library(devtools)
install_github("haoyu-yang001/Syn-GWAS",force = TRUE)
library(withr)
library(Matrix)
library(SynGWAS)

arraynum <- 500

GRMr <- 0.25
chunk_array <- 1
kinsetting <- "cousin"

n_obs <- 1000

tauT = 0.7
tauS = 0.4
sigmaT = 0.7
sigmaS = 0.5

pve_x = 0
maf <- 0.25

chunk_size <- 200  # Load 200 SNPs at a time
num_snps_each_array <- 400
num_chunks <- ceiling(num_snps_each_array / chunk_size)  # Total chunks

miss <- 0.5

GRM <- DGP_GRMr(n_obs,miss,GRMr)

hh <- find_blocks_vectorized(GRM)

independent_indices <- with_seed(25324, {
  sapply(hh, function(x) x[sample.int(length(x), 1)])
})

pre_vars_step1 <- DGP_GRMr_nullmodel_step1(n_obs,miss,tauT, tauS, sigmaT, sigmaS, GRM)

rho <- 0.75

cat("Processing GRMr", GRMr, "miss", miss, "rho", rho, "chunk array:", chunk_array, "\n")

pre_vars_step2 <- DGP_GRMr_nullmodel_step2(n_obs,miss,tauT, tauS, sigmaT, sigmaS,
                                           rho,GRM,pre_vars_step1)

mydf <- DGP_GRMr_nullmodel_step3(n_obs,miss,tauT, sigmaT, pve_x, pre_vars_step2, GRM)

SynSurrG_typeIerror_step1_pars <- SynSurrG_typeIerror_step1(mydf)
OracleG_typeIerror_step1_pars <- OracleG_typeIerror_step1(mydf)
ObsG_typeIerror_step1_pars <- ObsG_typeIerror_step1(mydf)

SynSurr_typeIerror_step1_pars <- SynSurr_typeIerror_step1(mydf,independent_indices)
Oracle_typeIerror_step1_pars <- Oracle_typeIerror_step1(mydf,independent_indices)
Obs_typeIerror_step1_pars <- Obs_typeIerror_step1(mydf,independent_indices)

SynSurrG_results <- list()
OracleG_results <- list()
ObsG_results <- list()

SynSurr_results <- list()
Oracle_results <- list()
Obs_results <- list()

for (chunk_index in seq_len(num_chunks)) {
  
  chunk_index <- 1
  geno_matrix <- replicate(
    chunk_size,generate_gmatrix_kinship(n_obs, miss, maf, kinsetting),
    simplify = "matrix")
  
  SynSurrG_results[[chunk_index]] <- score_test_SynSurrG_multiply(g_matrix = geno_matrix,step1_pars = SynSurrG_typeIerror_step1_pars)
  OracleG_results[[chunk_index]] <- score_test_OracleG_multiply(g_matrix = geno_matrix,step1_pars = OracleG_typeIerror_step1_pars)
  ObsG_results[[chunk_index]] <- score_test_ObsG_multiply(g_matrix = geno_matrix,step1_pars = ObsG_typeIerror_step1_pars)
  
  SynSurr_results[[chunk_index]] <- score_test_SynSurr_multiply(g_matrix = geno_matrix,step1_pars = SynSurr_typeIerror_step1_pars,independent_indices)
  Oracle_results[[chunk_index]] <- score_test_Oracle_multiply(g_matrix = geno_matrix,step1_pars = Oracle_typeIerror_step1_pars,independent_indices)
  Obs_results[[chunk_index]] <- score_test_Obs_multiply(g_matrix = geno_matrix,step1_pars = Obs_typeIerror_step1_pars,independent_indices)
  
  if (chunk_index %% 50 == 0) {
    print(chunk_index)
  }
}

merged_SynSurrG <- merge_results(SynSurrG_results)
merged_OracleG  <- merge_results(OracleG_results)
merged_ObsG     <- merge_results(ObsG_results)
merged_SynSurr  <- merge_results(SynSurr_results)
merged_Oracle   <- merge_results(Oracle_results)
merged_Obs      <- merge_results(Obs_results)

SynSurrG_negative_log10_pvalues <- merged_SynSurrG$negative_log10_pval_SynSurrG
#SynSurrG_hat_beta <- merged_SynSurrG$hat_beta_SynSurrG
#SynSurrG_var_hat_beta <- merged_SynSurrG$var_hat_beta_SynSurrG

OracleG_negative_log10_pvalues <- merged_OracleG$negative_log10_pval_OracleG
#OracleG_hat_beta <- merged_OracleG$hat_beta_OracleG
#OracleG_var_hat_beta <- merged_OracleG$var_hat_beta_OracleG

ObsG_negative_log10_pvalues <- merged_ObsG$negative_log10_pval_ObsG
#ObsG_hat_beta <- merged_ObsG$hat_beta_ObsG
#ObsG_var_hat_beta <- merged_ObsG$var_hat_beta_ObsG

SynSurr_negative_log10_pvalues <- merged_SynSurr$negative_log10_pval_SynSurr
#SynSurr_hat_beta <- merged_SynSurr$hat_beta_SynSurr
#SynSurr_var_hat_beta <- merged_SynSurr$var_hat_beta_SynSurr

Oracle_negative_log10_pvalues <- merged_Oracle$negative_log10_pval_Oracle
#Oracle_hat_beta <- merged_Oracle$hat_beta_Oracle
#Oracle_var_hat_beta <- merged_Oracle$var_hat_beta_Oracle

Obs_negative_log10_pvalues <- merged_Obs$negative_log10_pval_Obs
#Obs_hat_beta <- merged_Obs$hat_beta_Obs
#Obs_var_hat_beta <- merged_Obs$var_hat_beta_Obs

results <- cbind(OracleG_negative_log10_pvalues,
                 #OracleG_hat_beta,
                 #OracleG_var_hat_beta,
                 SynSurrG_negative_log10_pvalues,
                 #SynSurrG_hat_beta,
                 #SynSurrG_var_hat_beta,
                 ObsG_negative_log10_pvalues,
                 #ObsG_hat_beta,
                 #ObsG_var_hat_beta,
                 Oracle_negative_log10_pvalues,
                 #Oracle_hat_beta,
                 #Oracle_var_hat_beta,
                 SynSurr_negative_log10_pvalues,
                 #SynSurr_hat_beta,
                 #SynSurr_var_hat_beta,
                 Obs_negative_log10_pvalues)
#Obs_hat_beta,
#Obs_var_hat_beta)


library(data.table)
library(dplyr)
library(stringr)
library("mrbma")
library(tibble)

# clump across all instruments from the metabolites to get independent SNPs
isoleucine_PD <- fread("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/harmonised/meta_PD/f.x376_invnorm_GCST90200399_PD_harmonised.txt") %>% 
select(SNP,beta.exposure,se.exposure,pval.exposure)

leucine_PD <- fread("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/harmonised/meta_PD/f.x397_invnorm_GCST90200389_PD_harmonised.txt") %>% 
select(SNP,beta.exposure,se.exposure,pval.exposure)

all_snps <- rbind(isoleucine_PD,leucine_PD)

exposure_data_clumped <- all_snps %>%
       rename(rsid = SNP,
                  pval = pval.exposure) %>%
       ieugwasr::ld_clump(clump_r2 = 0.001,
                 clump_p = 5e-6,
                 clump_kb = 10000,
                 plink_bin = genetics.binaRies::get_plink_binary(), 
                 bfile = "/scratch/prj/proitsi/lachlan/PhD_project_2/reference_panel/EUR")

unique_snps <- unique(exposure_data_clumped$rsid)


# Extract these SNPs from the original GWAS files
X1 <- fread("/scratch/prj/proitsi/sumstats/metabolites_meta/chen_hysi_meta/f.x397_invnorm_GCST90200389_meta1.txt") %>% select(Effect,StdErr,MarkerName,`P-value`) %>% rename(SNP = MarkerName, betaX1 = Effect, seX1 = StdErr,p = `P-value`) %>%  filter(SNP %in% unique_snps)

X2 <- fread("/scratch/prj/proitsi/sumstats/metabolites_meta/chen_hysi_meta/f.x376_invnorm_GCST90200399_meta1.txt") %>% select(Effect,StdErr,MarkerName,`P-value`) %>% rename(SNP = MarkerName, betaX2 = Effect, seX2 = StdErr,p = `P-value`) %>%  filter(SNP %in% unique_snps)

X_all <- inner_join(X1,X2, by = "SNP") 

X_all <- X_all %>% relocate(SNP, .after = seX2)


# Extract from the outcome
Y_X1 <- fread("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/harmonised/meta_PD/f.x376_invnorm_GCST90200399_PD_harmonised.txt") %>% 
select(SNP,beta.outcome,se.outcome,pval.outcome) %>% filter(SNP %in% unique_snps)

Y_X2 <- fread("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/harmonised/meta_PD/f.x397_invnorm_GCST90200389_PD_harmonised.txt") %>% 
select(SNP,beta.outcome,se.outcome,pval.outcome) %>% filter(SNP %in% unique_snps)

Y_joined <- rbind(Y_X1,Y_X2)

Y_joined <- Y_joined %>% arrange(SNP) %>% rename(betaYG = beta.outcome, sebetaYG = se.outcome, pvalYG = pval.outcome)

Y_joined <- Y_joined %>% distinct()


# Create matrices for exposure_beta, exposure_se, and exposure_pval for BMA
exposure_beta <- as.matrix(X_all %>% select(SNP, betaX1, betaX2) %>% column_to_rownames(var = "SNP"))
exposure_se <- as.matrix(X_all %>% select(SNP, seX1, seX2) %>% column_to_rownames(var = "SNP"))
exposure_pval <- as.matrix(X_all %>% select(SNP, p.x, p.y) %>% column_to_rownames(var = "SNP"))
expname <- data.frame(id.exposure = c("leucine", "isoleucine"), exposure = c("Leucine", "Isoleucine"))

outcome_beta <- as.vector(Y_joined$betaYG)
outcome_se <- as.vector(Y_joined$sebetaYG)
outcome_pval <- as.vector(Y_joined$pvalYG)
outname <- data.frame(id.outcome = "PD", outcome = "PD")



# Create the harmonized data list for mr-bma
harmonized_data <- list(
  exposure_beta = exposure_beta,
  exposure_se = exposure_se,
  exposure_pval = exposure_pval,
  outcome_beta = outcome_beta,  # Assuming this is pre-defined or needs to be similarly structured
  outcome_se = outcome_se,      # Assuming this is pre-defined or needs to be similarly structured
  outcome_pval = outcome_pval,  # Assuming this is pre-defined or needs to be similarly structured
  expname = expname,
  outname = outname
)


# Run mr-bma
result <- mr_bma(harmonized_data,
                calculate_p = FALSE, 
                nrepeat = 10000, 
                remove_outliers = TRUE,
                remove_influential = TRUE,
                prior_prob = 0.1,
                prior_sigma = 0.25)

result$model_best

bma_results <- as.data.frame(result$mip_table[,1:2])

bma_results$outcome <- "PD"

bma_results <- bma_results %>% relocate(outcome, .after = exposure)

bma_results <- bma_results %>% rename(exposure = "rf")


# Write results
write.table(bma_results,"/scratch/prj/proitsi/lachlan/PhD_project_2/mrbma/out/PD_bma.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

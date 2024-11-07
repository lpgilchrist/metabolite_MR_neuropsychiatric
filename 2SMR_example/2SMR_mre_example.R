# module load r/4.1.3-gcc-10.3.0-withx-rmath-standalone-python3+-chk-version

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(stringr)
library(R.utils)
library(ieugwasr)
library(ggplot2)


#Set working directory
setwd("/scratch/prj/proitsi/sumstats/metabolites_meta/chen_hysi_meta/")


outcome_file<-fread("/scratch/prj/proitsi/lachlan/PhD_project_2/depression_sumstats/munged/Depression_MVP_UKB_PGC_Finngen1.txt")

freq <- fread("/scratch/prj/proitsi/lachlan/PhD_project_2/reference_panel/ref_freq.frq") %>% select(SNP,A1,A2,MAF)


# Extract the phenotype names
outcome_file$Allele1 <- toupper(outcome_file$Allele1)
outcome_file$Allele2 <- toupper(outcome_file$Allele2)

outcome_file$Phenotype <- "Depression"


# Set the exposure and outcome names 
outcome_name <- "Depression"

outcome_file$Z <- outcome_file$Effect/outcome_file$StdErr

  args <- commandArgs(trailingOnly = TRUE)

  file <- args[1] 

  print(file)

  exposure_name <- str_remove(file, pattern = "_meta1.txt")
 
  # read in GWAS file

# Read in the metabolites GWAS file
exposure_file <- fread(paste0("/scratch/prj/proitsi/sumstats/metabolites_meta/chen_hysi_meta/",file), header=T,data.table=F)

head(exposure_file)


# Set the exposure name
exposure_file$Phenotype <- exposure_name


# Set the path to the matched/proxy instrument file
use_instruments_loc <- paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/snappy_out/meta_proxies/",exposure_name,"_Depression_matched_proxies.txt")



# Read in the instrument rsids
use_instruments <- fread(use_instruments_loc, header = FALSE) %>% select(V2) %>% rename(SNP = V2)




# Extract thhe instruments from the exposure
exposure_data <- exposure_file %>% filter(MarkerName %in% use_instruments$SNP)


# Obtain the allele frequencies from 1KG for harmonisation step 
exposure_data$Allele1 <- toupper(exposure_data$Allele1)

exposure_data$Allele2 <- toupper(exposure_data$Allele2)



# Get allele frequencies from the reference panel as these are missing and are used for harmonisation and are missing from Hysi et al.
gwas_ref_freq_match<-merge(exposure_data,freq,by.x=c('MarkerName','Allele1','Allele2'),by.y=c('SNP','A1','A2'))

print(nrow(gwas_ref_freq_match))


gwas_ref_freq_switch<-merge(exposure_data,freq,by.x=c('MarkerName','Allele1','Allele2'),by.y=c('SNP','A2','A1'))

print(nrow(gwas_ref_freq_switch))


gwas_ref_freq_switch$MAF<-1-gwas_ref_freq_switch$MAF


gwas_ref_freq<-rbind(gwas_ref_freq_match, gwas_ref_freq_switch) 


names(gwas_ref_freq)[names(gwas_ref_freq) == 'MAF']<-'Freq1'


exposure_data <- gwas_ref_freq


nrow(exposure_data)



# Format the exposure data
exposure_data <- format_data(
  exposure_data,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "StdErr",
  eaf_col = "Freq1",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  samplesize_col = "N",
  min_pval = 1e-200,
  z_col = "Z",
  chr_col = "Chromosome",
  pos_col = "Position",
  log_pval = FALSE)


# Test instrument strength and exclude instruments with F < 10
exposure_data$F <- exposure_data$beta.exposure^2/exposure_data$se.exposure^2

print(exposure_data)

exposure_data <- exposure_data %>% filter(F >= 10)

print(exposure_data)

exposure_data <- exposure_data %>% select(-F)


# Format the outcome data
outcome_data <- format_data(
  outcome_file,
  type = "outcome",
  snps = exposure_data$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "StdErr",
  eaf_col = "Freq1",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  samplesize_col = "N",
  min_pval = 1e-200,
  z_col = "Z",
  chr_col = "Chromosome",
  pos_col = "Position",
  log_pval = FALSE)


# Harmonise the exposure and outcome, excluding pallindromic variants that cannot be aligned
harmonised_data <- harmonise_data(
    exposure_dat = exposure_data, 
    outcome_dat = outcome_data, action = 2)

harmonised_data <- harmonised_data %>% filter(mr_keep == "TRUE")

# Write the harmonised data for later use in sensitivity analyses
write.table(harmonised_data, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/harmonised/meta_depression/",exposure_name,"_Depression_harmonised.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

if (nrow(harmonised_data) < 5){
   sink("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/output/MR_mre/meta_depression/metab_mre_not_run.txt", append = TRUE)
    print(paste0("There are only ",nrow(harmonised_data)," instruments available for ",exposure_name," and Depression after harmonisation – did not run.."))
  sink()

}else{

nrow(harmonised_data)


# Run the IVW-MR (mre) analysis
mr_results <- mr(harmonised_data, method_list=c("mr_ivw_mre"))

print(mr_results)

# Obtain the F-statistic
mr_results$Mean_F <- mean(harmonised_data$beta.exposure^2/harmonised_data$se.exposure^2)
mr_results$Min_F <- min(harmonised_data$beta.exposure^2/harmonised_data$se.exposure^2)
mr_results$Max_F <- max(harmonised_data$beta.exposure^2/harmonised_data$se.exposure^2)

# Obtain I^2
mr_results$I_squared <- Isq(harmonised_data$beta.exposure,harmonised_data$se.exposure)


# Obtain 95% CI for beta
mr_results$lower_b <- (mr_results$b - 1.96*mr_results$se)

mr_results <- mr_results %>% relocate(lower_b, .before = b)

mr_results$upper_b <- (mr_results$b + 1.96*mr_results$se)

mr_results <- mr_results %>% relocate(upper_b, .after = b)


# Obtain OR and 95% CI
mr_results$lower_OR <- exp(mr_results$b - 1.96*mr_results$se)

mr_results$OR <- exp(mr_results$b)

mr_results$upper_OR <- exp(mr_results$b + 1.96*mr_results$se)

# Obtain the clumpng threshold used
if (max(harmonised_data$pval.exposure)<= 5e-8){
mr_results$clump_thres <- 5e-8

}else{mr_results$clump_thres <- 5e-6}

mr_results <- mr_results %>% relocate(clump_thres, .before =  nsnp)

mr_results <- mr_results %>% relocate(exposure, .before =  outcome)

mr_results <- mr_results %>% select(-id.exposure, -id.outcome)



write.table(mr_results, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/output/MR_mre/meta_depression/",exposure_name,"_depression_mre_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)}




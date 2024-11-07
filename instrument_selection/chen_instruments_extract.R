# module load r/4.1.3-gcc-10.3.0-withx-rmath-standalone-python3+-chk-version

library(remotes)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(R.utils)
library(ieugwasr)
library(TwoSampleMR)

# This script creates a a list of instruments in the risk factors and a list of snps that are shared between the depression GWAS and the risk factor
# We can use these output together in snappy to extract proxy instruments if they are missing from the outcome dataset

#Set working directory
setwd("/scratch/prj/proitsi/sumstats/Chen_metabolites/Chen_munged/")

ref_snps <- fread("/scratch/prj/proitsi/lachlan/PhD_project_2/reference_panel/ref_freq.frq") %>% select(SNP)

depression_snps <- fread("/scratch/prj/proitsi/lachlan/PhD_project_2/depression_sumstats/munged/Depression_MVP_UKB_PGC_Finngen1.txt") %>% select(MarkerName) %>% rename(SNP = MarkerName) 



exposure_files <- list.files(pattern=".txt.gz")

exposure_files 


for(i in exposure_files){
  # read in GWAS files

exposure_file<-fread(i)

exposure_name <- str_remove(i, pattern = "_munged.txt.gz")


# Extract the matching snps in order to obtain proxies later
exposure_file <- exposure_file %>% inner_join(ref_snps)

exposure_snps <- exposure_file %>% select(SNP)

print(paste0("The number of rows before joining is ",nrow(exposure_snps)))

matched_snps <- exposure_snps %>% inner_join(depression_snps)

print(paste0("The number of rows after joining is ",nrow(matched_snps)))

print(paste0("Writing SNPs..."))

write.table(matched_snps, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/snappy_snps/snappy_snps_chen/",exposure_name,"_Depression_snps.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


exposure_file <- exposure_file %>% filter(P<=5e-6)


print(paste0("Extracting IVs for ",exposure_name))



exposure_data <- format_data(
  exposure_file,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "FRQ",
  effect_allele_col = "A2",
  other_allele_col = "A1",
  pval_col = "P",
  samplesize_col = "N",
  min_pval = 1e-200,
  z_col = "Z",
  chr_col = "CHR",
  pos_col = "BP",
  log_pval = FALSE)


# Clump the exposure data
if (min(exposure_data$pval.exposure)<= 5e-8){
exposure_data_clumped <- exposure_data %>%
       rename(rsid = SNP,
                  pval = pval.exposure) %>%
       ieugwasr::ld_clump(clump_r2 = 0.001,
                 clump_p = 5e-8,
                 clump_kb = 10000,
                 plink_bin = genetics.binaRies::get_plink_binary(), 
                 bfile = "/scratch/prj/proitsi/lachlan/PhD_project_2/reference_panel/EUR")

}else{

  sink("/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/chen_instruments/chen_5e-6.txt", append = TRUE)
  print(paste0(exposure_name," has no IVs at 5e-8, using 5e-6"))
  sink()

  exposure_data_clumped <- exposure_data %>%
       rename(rsid = SNP,
                  pval = pval.exposure) %>%
       ieugwasr::ld_clump(clump_r2 = 0.001,
                 clump_p = 5e-6,
                 clump_kb = 10000,
                 plink_bin = genetics.binaRies::get_plink_binary(), 
                 bfile = "/scratch/prj/proitsi/lachlan/PhD_project_2/reference_panel/EUR")}


exposure_data_clumped <- exposure_data_clumped %>% rename(SNP = rsid, pval.exposure = pval)


print(nrow(exposure_data_clumped))

# Check for sufficient insturments
if (nrow(exposure_data_clumped)< 5){

  sink("/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/chen_instruments/chen_5e-6.txt", append = TRUE)
  print(paste0(exposure_name," has only ",nrow(exposure_data_clumped)," IVs at 5e-8, using 5e-6"))
  sink()

  exposure_data_clumped <- exposure_data %>%
       rename(rsid = SNP,
                  pval = pval.exposure) %>%
       ieugwasr::ld_clump(clump_r2 = 0.001,
                 clump_p = 5e-6,
                 clump_kb = 10000,
                 plink_bin = genetics.binaRies::get_plink_binary(), 
                 bfile = "/scratch/prj/proitsi/lachlan/PhD_project_2/reference_panel/EUR")


exposure_data_clumped <- exposure_data_clumped %>% rename(SNP = rsid, pval.exposure = pval)

}


print(nrow(exposure_data_clumped))


exposure_data_clumped <- exposure_data_clumped %>% select(SNP)

print(paste0("Writing instruments for snappy for ",exposure_name))

write.table(exposure_data_clumped, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/chen_instruments/",exposure_name,"_instruments.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

}



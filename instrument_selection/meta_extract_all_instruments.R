# module load r
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(R.utils)
library(ieugwasr)
library(TwoSampleMR)

# This script creates a a list of instruments
# We can use these output together in snappy to extract proxy instruments if they are missing from the outcome dataset

#Set working directory
setwd("/scratch/prj/proitsi/sumstats/metabolites_meta/chen_hysi_meta/")

exposure_files <- list.files(pattern="meta1")


exposure_files 


for(i in exposure_files){
  # read in GWAS files

exposure_file<-fread(i)

exposure_name <- str_remove(i, pattern = "_meta1.txt")

exposure_file$Allele1 <- toupper(exposure_file$Allele1)

exposure_file$Allele2 <- toupper(exposure_file$Allele2)

exposure_file <- exposure_file %>% filter(N == 17108)

exposure_file <- exposure_file %>% filter(`P-value`<=5e-6)



print(paste0("Extracting IVs for ",exposure_name))



exposure_data <- format_data(
  exposure_file,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "StdErr",
  eaf_col = "FRQ",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  samplesize_col = "N",
  min_pval = 1e-200,
  chr_col = "Chromosome",
  pos_col = "Position",
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

  sink("/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/meta_instruments/metab_meta_5e-6.txt", append = TRUE)
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

  sink("/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/meta_instruments/metab_meta_5e-6.txt", append = TRUE)
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

write.table(exposure_data_clumped, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/meta_instruments/",exposure_name,"_instruments.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

}



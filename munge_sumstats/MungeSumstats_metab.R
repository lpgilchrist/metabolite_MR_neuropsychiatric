## Clean summary statistics using the MungeSumstats package in R ##

library(MungeSumstats)
library(data.table)
library(dplyr)
library(stringr)
library(R.utils)

# Read in 1000 Genomes MAFs
#freq <- fread("/scratch/prj/ukbiobank/usr/Lachlan/polygenic_paper/COLOC_reporter/ld_reference/all_chr_EUR/plink.frq") %>% select(SNP,A1,A2,MAF)

# Set names for renaming to avoid issues

possible_colnames <- c(CHR = "chromosome",
  BP = "base_pair_location",
  SNP = "variant_id",
  EFFECT_ALLELE = "effect_allele",
  NON_EFFECT_ALLELE = "other_allele",
  FREQ = "effect_allele_frequency",
  BETA = "beta",
  SE = "sebeta", SE = "standard_error",
  P = "p_value")


#Set working directory

setwd("/scratch/prj/proitsi/sumstats/Chen_metabolites/")

data("sumstatsColHeaders")

args <- commandArgs(trailingOnly = TRUE)


  file <- args[1] 

  print(file)

# Extract trait name

  trait_name <- str_remove(file, pattern = "_buildGRCh38.tsv.gz")



# Read in the metabolites GWAS file
df<-fread(paste0("/scratch/prj/proitsi/sumstats/Chen_metabolites/",file), header=T,data.table=F)

# df <- fread("GCST90199625_buildGRCh38.tsv.gz")

print(paste0(file," has been read!"))


# Rename columns
print("Renaming columns...")
df <- df %>%
   dplyr::rename(any_of(possible_colnames))
print("Columns renamed!")

df$N <- 8299 

# Select required columns
print("Selecting columns...")
df <- df %>% select(SNP, CHR, BP, EFFECT_ALLELE, NON_EFFECT_ALLELE, FREQ, BETA, SE, P, N)
print(head(df))

# Remove all NAs as this will cause things to crash etc.
print(paste0("Removing NAs..."))
  df <- df[!is.na(df$BETA), ]
  df <- df[!is.na(df$SE), ]
print(paste0("NAs removed!!"))

# Filter allele freq
df <- df %>%
    filter(FREQ >= 0.01 & FREQ <= 0.99)


#Create Z column
  df$Z <- df$BETA/df$SE


  format_sumstats(
  path = df,
  ref_genome = "GRCh38", # Set reference genome build for the sumstats
  dbSNP = 144,
  convert_ref_genome = "GRCh37", # Set genome build for conversion should it not match (default = NULL)
  convert_small_p = TRUE, # For p-values outside of R's range, these are converted to 0
  convert_large_p = TRUE, # For p-value over 1, these are converted to 1
  convert_neg_p = TRUE, # For p-values less that 0 (i.e. negative) these are converted to 0
  compute_z = FALSE, # Whether to convert beta and se to z-score (not required, set as FALSE)
  force_new_z = FALSE, # For if z column already exists (not required, set as FALSE)
  compute_n = 0L, # Compute missing N for SNPs (not required, set as 0)
  convert_n_int = TRUE, # If N is not an integer, this is rounded
  impute_beta = FALSE, # Impute missing effect value if beta column missing (not required)
  impute_se = FALSE, # Impute missing effect standard error value if column missing (not required) 
  analysis_trait = NULL, # For if multiple traits are contained within the sumstats (not required)
  INFO_filter = 0.7, # If INFO column is present, filter on specified value
  FRQ_filter = 0, # No filtering applied if set to 0. Instead filter in MAF column seperately.
  pos_se = TRUE, # Check all standard errors are positive, if not remove those that arent
  effect_columns_nonzero = TRUE, # Removes SNPs with beta effects at 0, set to FALSE if some SNPs may legitimately have zero effect
  N_std = 5, # Remove SNPs with >5 SDs above mean N
  N_dropNA = FALSE, # Drop rows where N is missing
  rmv_chrPrefix = TRUE, # Drop chr or CHR prefix is present in CHR column (i.e. chr1 etc.)
  on_ref_genome = TRUE, # Checks all SNPs are on the reference genome and imputes if missing
  strand_ambig_filter = FALSE, # Remove strand ambigous SNPs
  allele_flip_check = TRUE, # Check reference allele (A1 here) requires flipping
  allele_flip_drop = TRUE, # If neither allele matches reference genome, they are dropped
  allele_flip_z = TRUE, # Flip the Z score column along with the allele
  allele_flip_frq = TRUE, # Flip the effect allele column
  bi_allelic_filter = TRUE, # Remove non-bialleic SNPs
  snp_ids_are_rs_ids = TRUE, # Defines whether SNPs IDs should be inferred as RSIDs
  remove_multi_rs_snp = TRUE, # For SNPs with more than 1 RSID, keep only the first (i.e. if rs123_rs456 keep only rs123)
  frq_is_maf = TRUE, # Set as true, stops renaming of column if major allele freq is inferred
  indels = TRUE, # Exclude indels
  sort_coordinates = TRUE, # Sort the coordinates of the resulting summary statistics
  nThread = 1, # Number of threads for analysis
  save_path = paste0("/scratch/prj/proitsi/sumstats/Chen_metabolites/Chen_munged/",trait_name,"_munged.txt.gz"),
  write_vcf = FALSE, # Write vcf format file (set to FALSE)
  tabix_index = FALSE, # Related to vcf command above, ignore
  return_data = TRUE, # Return data in a specified format (see next command)
  return_format = "data.table", # Return format as data table
  ldsc_format = FALSE, # If requiring format for use straight away in ldsc
  log_folder_ind = FALSE, # For if log of filtered out snps are required
  log_mungesumstats_msgs = TRUE, # For a log of all messages from the above commands
  log_folder = paste0("/scratch/prj/proitsi/sumstats/Chen_metabolites/Chen_munged/munged_logs/",trait_name,"_logs"),  
  imputation_ind = FALSE, # If columns should be added to the new sumstats specifying filter steps
  force_new = FALSE, # 
  mapping_file = sumstatsColHeaders
)

print("Done")



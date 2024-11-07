#module load test_switch_kcl
#source test_switch

#module load r/4.1.3-gcc-10.3.0-withx-rmath-standalone-python3+-chk-version

# Load the required packages
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(stringr)
library(R.utils)
library(ieugwasr)
library(ggplot2)
library(MRcML)

mre_meta_results = data.frame()
mre_meta_dep_het = data.frame()
mre_meta_dep_pleio = data.frame()
all_single_results = data.frame()
all_leave_out_results = data.frame()

setwd("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/output/MR_sensitivity/all_sensitivity/")


outcomes <- c("depression","AD","ALS","ANX","BP","SCZ","MS","PD")

# Before running this analysis, a file containing the FDR significant metabolites from IVW-MR has been written for each outcome
for (i in outcomes){
# Extract the FDR significant exposures
sensitivity_exposure <- fread(paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/output/MR_mre/meta_",i,"/meta_",i,"_FDR_sig_for_sensitivity.txt"))

# Get the trait names for the table later
gwas_names <- sensitivity_exposure %>% select(exposure,GWAS) %>% rename(GWAS = exposure, exposure = GWAS)

# Get the unique exposure values
exposure_pattern <- sensitivity_exposure$GWAS %>% unique()

# Construct the harmonised file names for pattern matching
harmonised_file_names <- paste0(exposure_pattern,"_",i,"_harmonised.txt")

# Set the path to the harmonised files directory
harmonised_filepath <- paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/harmonised/meta_",i)

# Get the full file paths for the harmonised files that match the pattern
harmonised_files <- list.files(path = harmonised_filepath, pattern = paste(harmonised_file_names, collapse = "|"), full.names = TRUE)


# Now you have the filtered exposure_files list
print(harmonised_files)




outcome_name <- i


for(j in harmonised_files){
 
  # read in GWAS file
	harmonised_file<-fread(j)

  print(harmonised_file)

  # Set the outcome name
	exposure_name <- str_remove(j, pattern = "_harmonised.txt")


	exposure_name <- str_remove(exposure_name, pattern = paste0("_",i))



mr_results <- mr(harmonised_file, method_list=c("mr_ivw_mre",
												"mr_penalised_weighted_median",
												"mr_weighted_median",
												"mr_egger_regression"))

print(mr_results)

mr_results <- mr_results %>% select(-id.outcome, -id.exposure)

harmonised_file <- harmonised_file %>% filter(mr_keep == "TRUE")

# Extract relevant fields for running contamination mixture and lasso MR in the MendelianRandomization package

conmix_input <- MendelianRandomization::mr_input(bx = harmonised_file$beta.exposure,
                                            bxse = harmonised_file$se.exposure, by = harmonised_file$beta.outcome, byse = harmonised_file$se.outcome,
                                            exposure = harmonised_file$exposure[1], outcome = harmonised_file$outcome[1], snps = harmonised_file$SNP,
                                            effect_allele = harmonised_file$effect_allele.exposure, other_allele = harmonised_file$other_allele.exposure,
                                            eaf = harmonised_file$eaf.exposure)

# Run conmix and lasso

conmix_output <- (MendelianRandomization::mr_conmix(conmix_input, alpha = 0.05))

conmix_results <- data.frame(exposure = harmonised_file$exposure[1], outcome = outcome_name, method = "Contamination Mixture",nsnp = length(conmix_output$ValidSNPs),b = (conmix_output$Estimate), pval = (conmix_output$Pvalue), se = (conmix_output@CIUpper[1] - conmix_output@CILower[1])/3.92)

lasso_output <- MendelianRandomization::mr_lasso(conmix_input, alpha = 0.05)



# lasso_output <- (MendelianRandomization::mr_lasso(conmix_input, alpha = 0.05, lambda = 0.8))

lasso_results <- data.frame(exposure = harmonised_file$exposure[1], outcome = outcome_name, method = "Lasso",nsnp = length(lasso_output$ValidSNPs),b = (lasso_output$Estimate), pval = (lasso_output$Pvalue), se = (lasso_output@CIUpper[1] - lasso_output@CILower[1])/3.92)


# Run the constrained maximum likelihood MR 

cML_MA_BIC_output <- MRcML::mr_cML(harmonised_file$beta.exposure,
                          harmonised_file$beta.outcome,
                          harmonised_file$se.exposure,
                          harmonised_file$se.outcome,
                          n = mean(harmonised_file$samplesize.exposure),
                          random_start = 100,
                          random_seed = 1)


cML_results <- data.frame(exposure = harmonised_file$exposure[1], outcome = outcome_name, method = "Constrained Maximum Likelihood",nsnp = (nrow(harmonised_file)) - (length(cML_MA_BIC_output$BIC_invalid)),b = (cML_MA_BIC_output$MA_BIC_theta), pval = (cML_MA_BIC_output$MA_BIC_p), se = (cML_MA_BIC_output$MA_BIC_se))



mr_results <- rbind(mr_results,conmix_results)

mr_results <- rbind(mr_results,lasso_results)

mr_results <- rbind(mr_results,cML_results)



print(mr_results)



mr_results$lower_OR <- exp(mr_results$b - 1.96*mr_results$se)

mr_results$OR <- exp(mr_results$b)

mr_results$upper_OR <- exp(mr_results$b + 1.96*mr_results$se)



mr_results$lower_b <- (mr_results$b - 1.96*mr_results$se)

mr_results$upper_b <- (mr_results$b + 1.96*mr_results$se)


num_pass <- sum(mr_results$pval < 0.05)

mr_results$sensitivity <- ifelse(num_pass >= 5 && length(unique(sign(mr_results$b))) == 1, "PASS", "FAIL")

# Set the exposure names to reflect
mr_results <- mr_results%>% inner_join(gwas_names) %>% mutate(exposure = ifelse(!is.na(GWAS), GWAS, exposure)) %>% select(-GWAS)



mre_meta_results = rbind(mre_meta_results,mr_results)


mr_het <- mr_heterogeneity(harmonised_file)

print(mr_het) 

# Set the exposure names to reflect
mr_het <- mr_het %>% inner_join(gwas_names) %>% mutate(exposure = ifelse(!is.na(GWAS), GWAS, exposure)) %>% select(-GWAS)


mr_pleio <- mr_pleiotropy_test(harmonised_file)

print(mr_pleio) 

# Set the exposure names to reflect
mr_pleio <- mr_pleio %>% inner_join(gwas_names) %>% mutate(exposure = ifelse(!is.na(GWAS), GWAS, exposure)) %>% select(-GWAS)



mre_meta_dep_het <- rbind(mre_meta_dep_het,mr_het)


mre_meta_dep_pleio <- rbind(mre_meta_dep_pleio,mr_pleio)



# Change the name of the harmonised file for plotting 
harmonised_file <- harmonised_file %>% inner_join(gwas_names) %>% mutate(exposure = ifelse(!is.na(GWAS), GWAS, exposure)) %>% select(-GWAS)

# Set name for plotting to avoid creating directories due to / etc

harmonised_file_exposure <- harmonised_file$exposure[1]

harmonised_file_exposure <- str_remove_all(harmonised_file_exposure,"/")

harmonised_file_exposure <- str_replace_all(harmonised_file_exposure,",","_")


harmonised_file_exposure <- str_replace_all(harmonised_file_exposure," ","_")




mr_leaveoneout <- mr_leaveoneout(harmonised_file, method = mr_ivw_mre)

mr_leaveoneout

mr_leave_plot <- mr_leaveoneout_plot(mr_leaveoneout)

ggsave(mr_leave_plot[[1]], filename = paste0(harmonised_file_exposure,"_",outcome_name,"_mr_leave_plot.png"), width = 7, height = 10)

all_leave_out_results <- rbind(all_leave_out_results,mr_leaveoneout)





mr_single_variant <- mr_singlesnp(harmonised_file, all_method = c("mr_ivw_mre","mr_egger_regression"))

mr_funnel <- mr_funnel_plot(mr_single_variant)

ggsave(mr_funnel[[1]], filename = paste0(harmonised_file_exposure,"_",outcome_name,"_mr_funnel_plot.png"), width = 7, height = 10)

all_single_results <- rbind(all_single_results,mr_single_variant)





mr_single <- mr_forest_plot(mr_single_variant)

ggsave(mr_single[[1]], filename = paste0(harmonised_file_exposure,"_",outcome_name,"_mr_forest_plot.png"), width = 7, height = 10)



}


}




# After everything has run write the outputs

write.table(mre_meta_results, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/output/MR_sensitivity/all_sensitivity/meta_sensitivity_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


write.table(mre_meta_dep_het, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/output/MR_sensitivity/all_sensitivity/meta_het_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


write.table(mre_meta_dep_pleio, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/output/MR_sensitivity/all_sensitivity/meta_pleio_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


write.table(all_leave_out_results, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/output/MR_sensitivity/all_sensitivity/meta_leave_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


write.table(all_single_results, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/output/MR_sensitivity/all_sensitivity/meta_single_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# Select sensitivity passing results for later
sensitivity_pass <- mre_meta_results %>% filter(sensitivity == "PASS")

write.table(sensitivity_pass, file = paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/2SMR/output/MR_sensitivity/all_sensitivity/meta_sensitivity_pass.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)






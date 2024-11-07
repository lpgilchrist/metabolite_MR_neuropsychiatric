library(data.table)
library(dplyr)
library(stringr)

setwd("/scratch/prj/proitsi/lachlan/PhD_project_2/PRS_metabs/outcome_phenos/")

# Load data
ukb_prs <- fread("/scratch/prj/proitsi/lachlan/PhD_project_2/PRS_metabs/outcome_phenos/metab_prs_phenos_combined_corrected_names.txt")
qc_exclusions_applied <- fread("/scratch/prj/ukbiobank/ukb82087/qc/ukb82087_post_qc_id_list_relatives_included.fam") %>% 
  rename(FID = V1) %>% select(FID)
related <- fread("/scratch/prj/proitsi/lachlan/PhD_project_2/PRS_metabs/outcome_phenos/rel_exclusions.txt") %>% select(V1) %>% rename(FID = V1)


# Exclude related individuals
ukb_prs <- ukb_prs %>% inner_join(qc_exclusions_applied, by = c("FID"))
ukb_prs <- ukb_prs %>% anti_join(related, by = c("FID"))

# Define outcome columns
outcome_cols <- names(ukb_prs)[25:32]

# Initialize results data table
all_results <- data.table()

# Loop through each outcome
for (i in outcome_cols) {
  
  if (i == "AD_occurance") {
    metab_cols <- names(ukb_prs)[2:4]
    ukb_test_data <- ukb_prs %>% filter(age >= 60)
  } else if (i == "ALS_occurance") {
    metab_cols <- names(ukb_prs)[5]
    ukb_test_data <- ukb_prs  
  } else if (i == "ANX_occurance") {
    metab_cols <- names(ukb_prs)[6:7]
    ukb_test_data <- ukb_prs
  } else if (i == "BP_occurance") {
    metab_cols <- names(ukb_prs)[8:10]
    ukb_test_data <- ukb_prs
  } else if (i == "DEP_occurance") {
    metab_cols <- names(ukb_prs)[11:13]
    ukb_test_data <- ukb_prs
  } else if (i == "MS_occurance") {
    metab_cols <- names(ukb_prs)[14:19]
    ukb_test_data <- ukb_prs 
  } else if (i == "PD_occurance") {
    metab_cols <- names(ukb_prs)[20:21]
    ukb_test_data <- ukb_prs %>% filter(age >= 50)
  } else if (i == "SCZ_occurance") {
    metab_cols <- names(ukb_prs)[22:24]
    ukb_test_data <- ukb_prs
  }


  # Create condition for all controls to have no occurance of any disorder
  condition_string <- paste0("!((", i, " == 0) & (", paste(outcome_cols[outcome_cols != i], " == 1", collapse = " | "), "))")

  condition_expr <- rlang::parse_expr(condition_string)

  # Apply filtering
  ukb_test_data <- ukb_test_data %>% filter(!!condition_expr)

  if (i == "DEP_occurance" | i == "ANX_occurance"){

  ukb_test_data <- ukb_test_data %>% filter(SCZ_occurance != 1, BP_occurance != 1)

  }

  
  # Loop through each metabolite
  for (j in metab_cols) {
    null_formula <- paste0(i, " ~ age + sex + as.factor(batch) + as.factor(assessment_centre) + 
                            pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")
    null_model <- glm(formula = null_formula, family = binomial, data = ukb_test_data)
    
    full_formula <- paste0(i, " ~ scale(", j,") + age + sex + as.factor(batch) + as.factor(assessment_centre) + 
                            pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")
    full_model <- glm(formula = full_formula, family = binomial, data = ukb_test_data)
    
    # Calculate Nagelkerke's R-squared
    r2 <- fmsb::NagelkerkeR2(full_model)$R2 - fmsb::NagelkerkeR2(null_model)$R2
    
    # Extract summary statistics
    summary_stats <- summary(full_model)
    N <- sum(!is.na(ukb_test_data[[i]]))
    N_cases <- sum(ukb_test_data[[i]] == 1)
    N_controls <- N-N_cases
    N_eff <- 4/(1/N_cases + 1/N_controls)

    filtered_nas <- ukb_test_data[!is.na(ukb_test_data[[i]]), ]
    age = mean(filtered_nas$age, na.rm = TRUE)
    percent_male <- sum(filtered_nas$sex == 1)/sum(!is.na(filtered_nas$sex))*100

    
    # Extract coefficients, standard errors, and p-values
    coef_estimate <- summary_stats$coefficients[paste0("scale(",j,")"), "Estimate"]
    se <- summary_stats$coefficients[paste0("scale(",j,")"), "Std. Error"]
    p_value <- summary_stats$coefficients[paste0("scale(",j,")"), "Pr(>|z|)"]
    
    # Create a data frame with results
    metab_prs_results <- data.frame(
      Metabolite = j,
      Outcome = i,
      N = N,
      N_cases = N_cases,
      N_controls = N_controls,
      N_eff = round(N_eff),
      percent_male = percent_male,
      mean_age = age,
      r2 = r2,
      estimate = coef_estimate,
      se = se,
      p_value = p_value
    )

    print(metab_prs_results)
    
    # Bind results to all_results data frame
    all_results <- rbind(all_results, metab_prs_results)
  }
}



p <- as.numeric(all_results$p_value)

fdr <- p.adjust(p, method = "fdr")

all_results$p_FDR <- fdr

all_results$OR <- exp(all_results$estimate)

all_results$OR_lower <- exp(all_results$estimate - 1.96*se)

all_results$OR_upper <- exp(all_results$estimate + 1.96*se)

all_results <- all_results %>% relocate(c(p_value,p_FDR), .after = OR_upper)

all_results <- all_results %>% relocate(c(OR_lower), .before = OR)

write.table(all_results, "/scratch/prj/proitsi/lachlan/PhD_project_2/PRS_metabs/outcome_phenos/metab_prs_results_glm.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)







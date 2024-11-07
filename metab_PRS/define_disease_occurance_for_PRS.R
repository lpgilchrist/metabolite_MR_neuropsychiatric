library(data.table)
library(dplyr)

setwd("/datasets/ukbiobank/ukb82087/phenotypes/")

withdrawal <- fread("/datasets/ukbiobank/ukb82087/raw/w82087_20240315.csv")


# Step 1: Read only the header to identify column names
age_sex_batch_centre <- fread("ukb50682.csv", select = c("eid", "31-0.0", "21022-0.0", "22000-0.0", "54-0.0", paste0("22009-0.", 1:10))) %>%
	   rename(
    sex = `31-0.0`,
    age = `21022-0.0`,
    batch = `22000-0.0`,
    assessment_centre = `54-0.0`,
    pc1 = `22009-0.1`,
    pc2 = `22009-0.2`,
    pc3 = `22009-0.3`,
    pc4 = `22009-0.4`,
    pc5 = `22009-0.5`,
    pc6 = `22009-0.6`,
    pc7 = `22009-0.7`,
    pc8 = `22009-0.8`,
    pc9 = `22009-0.9`,
    pc10 = `22009-0.10`
  )

age_sex_batch_centre <- age_sex_batch_centre[complete.cases(age_sex_batch_centre), ]

ukb <- fread("ukb675453.csv", select = c("eid","131022-0.0","130836-0.0","131042-0.0","42028-0.0","130892-0.0","130874-0.0","130906-0.0","130894-0.0")) %>%
         rename(PD_occurance = `131022-0.0`,
         		AD_occurance = `130836-0.0`,
         		MS_occurance = `131042-0.0`,
         		ALS_occurance = `42028-0.0`,
         		BP_occurance = `130892-0.0`,
         		SCZ_occurance = `130874-0.0`,
         		ANX_occurance = `130906-0.0`,
         		DEP_occurance = `130894-0.0`)


# Define exclusion dates
exclusion_dates <- c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07")

# List of columns to process
date_first_occurrence <- c("PD_occurance", "AD_occurance", "MS_occurance", "ALS_occurance", "BP_occurance", "SCZ_occurance", "ANX_occurance", "DEP_occurance")

for (col in date_first_occurrence) {
  # Set controls to 0
  ukb[[col]][is.na(ukb[[col]])] <- 0
  
  # Set values in exclusion_dates to NA
  ukb[[col]][ukb[[col]] %in% exclusion_dates] <- NA

  # Set all remaining non-zero and non-NA values to 1
  ukb[[col]][!is.na(ukb[[col]]) & ukb[[col]] != 0] <- 1
}

ukb_age_sex_batch_centre <- inner_join(age_sex_batch_centre,ukb)

ukb_age_sex_batch_centre <- ukb_age_sex_batch_centre %>% filter(!eid %in% withdrawal$V1)


write.table(ukb_age_sex_batch_centre, file=paste0("/scratch/prj/proitsi/lachlan/PhD_project_2/PRS_metabs/outcome_phenos/neuropsychiatric_phenos_for_PRS.txt"),sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





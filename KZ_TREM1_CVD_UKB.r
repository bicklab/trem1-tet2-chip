library(broom)
library(scales)
#library(stringr)
library(arrow)
library(survival)
library(tidyr)
library(dplyr)
library(haven)
#library(reshape2)
#library(tidyverse)
library(data.table)
library(lubridate)
library(ggplot2)

x <- load("UKB_data_all_pheno_20250228.Rdata")
x

dim(data_all_pheno_cox2)

#Formula
analyze_cox_model_refchip0 <- function(data, time_var, event_var, group_var, covariates) {
  df_filtered <- data %>%
    filter(chip == 0 | .data[[group_var]] == 1) %>%
    filter(.data[[time_var]] > 0)  
    
  if (nrow(df_filtered) != length(df_filtered[[time_var]]) || 
      nrow(df_filtered) != length(df_filtered[[event_var]])) {
    stop("Time and status are different lengths after filtering. Please check your data.")
  }
  
  if (!all(unique(df_filtered[[event_var]]) %in% c(0, 1))) {
    stop("Event variable must be binary (0 or 1).")
  }
  
  surv_obj <- Surv(time = df_filtered[[time_var]], event = df_filtered[[event_var]] == 1)
  
  covariate_formula <- paste(covariates, collapse = " + ")
  formula <- as.formula(paste("surv_obj ~", group_var, "+", covariate_formula))

  cox_fit <- coxph(formula, data = df_filtered)
  
  summary_cox <- summary(cox_fit)
  print(summary_cox)
  
  results <- as.data.frame(summary_cox$coefficients)
  colnames(results) <- c("Coefficient", "Exp(Coefficient)", "Standard Error", "z value", "p value")
  print(results)
  
  return(results)
}

batch_logistic_regression <- function(data, snp_list, outcome, covariates = NULL) {
  results <- list()
  
  base_formula <- paste0(outcome, " ~ ")
  
  for (snp in snp_list) {
    if (!is.null(covariates) && length(covariates) > 0) {
      formula <- as.formula(paste0(base_formula, snp, " + ", paste(covariates, collapse = " + ")))
    } else {
      formula <- as.formula(paste0(base_formula, snp))
    }
    

    print(formula)
    
    tryCatch({
      logistic_model <- glm(formula, data = data, family = binomial)
      summary_model <- summary(logistic_model)
      
      coef <- summary_model$coefficients[snp, "Estimate"]
      se <- summary_model$coefficients[snp, "Std. Error"]
      p_value <- summary_model$coefficients[snp, "Pr(>|z|)"]
      ci_lower <- exp(coef - 1.96 * se)
      ci_upper <- exp(coef + 1.96 * se)
      or <- exp(coef)
      
      results[[snp]] <- data.frame(
        SNP = snp,
        Beta = coef,
        OR = or,
        SE = se,
        CI_Lower = ci_lower,
        CI_Upper = ci_upper,
        P_Value = p_value
      )
    }, error = function(e) {
      warning(paste("Error in model for SNP:", snp, "- Skipping."))
    })
  }
  
  results_df <- do.call(rbind, results)
  return(results_df)
}

batch_logistic_by_outcome <- function(data, snp, outcome_list, covariates = NULL) {
  results <- list()
  
  for (outcome in outcome_list) {
    if (!is.null(covariates) && length(covariates) > 0) {
      formula <- as.formula(
        paste0(outcome, " ~ ", snp, " + ", paste(covariates, collapse = " + "))
      )
    } else {
      formula <- as.formula(paste0(outcome, " ~ ", snp))
    }
    
    print(formula)
    
    tryCatch({
      logistic_model <- glm(formula, data = data, family = binomial)
      summary_model <- summary(logistic_model)
      
      coef <- summary_model$coefficients[snp, "Estimate"]
      se <- summary_model$coefficients[snp, "Std. Error"]
      p_value <- summary_model$coefficients[snp, "Pr(>|z|)"]
      ci_lower <- exp(coef - 1.96 * se)
      ci_upper <- exp(coef + 1.96 * se)
      or <- exp(coef)
      
      results[[outcome]] <- data.frame(
        Outcome = outcome,
        SNP = snp,
        Beta = coef,
        OR = or,
        SE = se,
        CI_Lower = ci_lower,
        CI_Upper = ci_upper,
        P_Value = p_value
      )
    }, error = function(e) {
      warning(paste("Error in model for outcome:", outcome, "- Skipping."))
    })
  }
  
  results_df <- do.call(rbind, results)
  return(results_df)
}

table(data_all_pheno_cox2$chip)

y <- load("Chip_data_BUAT_0227.Rdata")
y

tet2_chip_ukb <- filter(chipcall_all, chipcall_all$cohort == "UKB" & chipcall_all$Gene.refGene == "TET2")
dim(tet2_chip_ukb)

data_all_pheno_cox2$tet2_chip <- NA
data_all_pheno_cox2$tet2_chip <- ifelse(data_all_pheno_cox2$ID_VUMC %in% tet2_chip_ukb$person_id, 1, 0)
table(data_all_pheno_cox2$tet2_chip)

table(data_all_pheno_cox2$CVD)

# TREM1

score_kz_rint <- readRDS("omicsscore_UKB_rint_KZ.rds")

trem1 <- score_kz_rint[,c("ID", "Triggering_receptor_expressed_on_myeloid_cells_1")]

data_all_pheno_cox2$ID_VUMC

data_all_pheno_cox2 <- data_all_pheno_cox2 %>% mutate(surv_ASCVD = pmin(surv_PVD,surv_AP,surv_hf,surv_AMI,surv_AHD,surv_death))
data_all_pheno_cox2 <- data_all_pheno_cox2 %>% mutate(surv_CAD = pmin(surv_AP,surv_AMI,surv_AHD,surv_death))

data_ana_trem1 <- data_all_pheno_cox2[,c("ID_VUMC","baseline_age", "age2", "genetic_sex", "smoking_0",
                                         "PC1","PC2","PC3","PCD4","PC5","BMI_0","baseline_HTN","baseline_DM",
                                         "systolicBP_0","diastolicBP_0","ldl_0","chip","tet2_chip","CVD","surv_PVD","surv_AP","surv_CVD",
                                         "surv_hf","surv_str","surv_AMI","surv_AHD","surv_ASCVD","surv_CAD","Coronary atherosclerosis [Atherosclerotic heart disease]",
                                        "Myocardial infarction [Heart attack]","Stroke","Heart failure","Angina pectoris","Peripheral vascular disease","death","surv_death")]

data_ana_trem1$ASCVD <- NA
data_ana_trem1$ASCVD <- ifelse(data_ana_trem1$`Coronary atherosclerosis [Atherosclerotic heart disease]`==1|
                        data_ana_trem1$`Myocardial infarction [Heart attack]`==1 |
                        data_ana_trem1$`Heart failure`==1|
                        data_ana_trem1$`Angina pectoris`==1|
                        data_ana_trem1$`Peripheral vascular disease`==1 |
                        data_ana_trem1$death==1,1, 0)
table(data_ana_trem1$ASCVD)

data_ana_trem1$CAD <- NA
data_ana_trem1$CAD <- ifelse(data_ana_trem1$`Coronary atherosclerosis [Atherosclerotic heart disease]`==1|
                        data_ana_trem1$`Myocardial infarction [Heart attack]`==1 |
                        data_ana_trem1$`Angina pectoris`==1|
                        data_ana_trem1$death==1,1, 0)
table(data_ana_trem1$CAD)

trem1$ID <- as.character(trem1$ID)
data_ana_trem1$ID_VUMC <- as.character(data_ana_trem1$ID_VUMC)

data_ana_pheno_trem1 <- merge(data_ana_trem1, trem1, by.x = "ID_VUMC", by.y = "ID", all = F)
dim(data_ana_pheno_trem1)

data_ana_pheno_trem1$CVD0 <- NA
data_ana_pheno_trem1$CVD0 <- ifelse(data_ana_pheno_trem1$`Coronary atherosclerosis [Atherosclerotic heart disease]`==1|
                        data_ana_pheno_trem1$`Myocardial infarction [Heart attack]`==1 |
                        data_ana_pheno_trem1$`Heart failure`==1|
                        data_ana_pheno_trem1$`Angina pectoris`==1|
                        data_ana_pheno_trem1$`Peripheral vascular disease`==1 |
                        data_ana_pheno_trem1$Stroke==1,1, 0)
table(data_ana_pheno_trem1$CVD0)

data_ana_pheno_trem1 <- data_ana_pheno_trem1 %>% mutate(surv_CVD0 = pmin(surv_PVD,surv_AP,surv_hf,surv_str,surv_AMI,surv_AHD))

saveRDS(data_ana_pheno_trem1, file = "data_ana_pheno_trem1.rds")

data_ana_pheno_trem1 <- readRDS("data_ana_pheno_trem1.rds")

dim(data_ana_pheno_trem1)

table(data_ana_pheno_trem1$chip)
table(data_ana_pheno_trem1$tet2_chip)

summary(data_ana_pheno_trem1$baseline_age)

data_nonchip_tet2chip <- filter(data_ana_pheno_trem1, data_ana_pheno_trem1$chip == 0 | data_ana_pheno_trem1$tet2_chip == 1)
dim(data_nonchip_tet2chip)

table(data_nonchip_tet2chip$tet2_chip)

analyze_cox_model <- function(data, time_var, event_var, group_var, covariates) {
  df_filtered <- data %>% filter(.data[[time_var]] > 0)
    
  if (nrow(df_filtered) != length(df_filtered[[time_var]]) || 
      nrow(df_filtered) != length(df_filtered[[event_var]])) {
    stop("Time and status are different lengths after filtering. Please check your data.")
  }
  
  if (!all(unique(df_filtered[[event_var]]) %in% c(0, 1))) {
    stop("Event variable must be binary (0 or 1).")
  }
  
  surv_obj <- Surv(time = df_filtered[[time_var]], event = df_filtered[[event_var]] == 1)
  
  covariate_formula <- paste(covariates, collapse = " + ")
  formula <- as.formula(paste("surv_obj ~", group_var, "+", covariate_formula))
  
  cox_fit <- coxph(formula, data = df_filtered)
  
  summary_cox <- summary(cox_fit)
  print(summary_cox)

  results <- as.data.frame(summary_cox$coefficients)
  colnames(results) <- c("Coefficient", "Exp(Coefficient)", "Standard Error", "z value", "p value")
  print(results)
  
  return(results)
}

results <- analyze_cox_model(
  data = data_nonchip_tet2chip, 
  time_var = "surv_CVD", 
  event_var = "CVD", 
  group_var = "tet2_chip", 
  covariates = c("baseline_age","age2", "genetic_sex", "smoking_0", "PC1", "PC2", "PC3", "PCD4", "PC5","ldl_0","baseline_DM","BMI_0","baseline_HTN")
)

results

results <- analyze_cox_model(
  data = data_nonchip_tet2chip, 
  time_var = "surv_CVD0", 
  event_var = "CVD0", 
  group_var = "tet2_chip", 
  covariates = c("baseline_age","age2", "genetic_sex", "smoking_0", "PC1", "PC2", "PC3", "PCD4", "PC5","ldl_0","baseline_DM","BMI_0","baseline_HTN")
)

logistic_mca <- batch_logistic_regression(
  data = data_ana_pheno_trem1,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome = "CVD",
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0", "PC1","PC2","PC3","PCD4","PC5","BMI_0","baseline_DM","systolicBP_0","diastolicBP_0","ldl_0")
)
logistic_mca

results <- analyze_cox_model(
  data = data_ana_pheno_trem1, 
  time_var = "surv_CVD", 
  event_var = "CVD", 
  group_var = "Triggering_receptor_expressed_on_myeloid_cells_1", 
  covariates = c("baseline_age","age2", "genetic_sex", "smoking_0", "PC1", "PC2", "PC3", "PCD4", "PC5","ldl_0","baseline_DM","BMI_0","baseline_HTN")
)

results

data_ana_pheno_trem1_chip <- filter(data_ana_pheno_trem1, data_ana_pheno_trem1$chip == 1)

logistic_mca <- batch_logistic_regression(
  data = data_ana_pheno_trem1_chip,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome = "CVD",
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0", "PC1","PC2","PC3","PCD4","PC5","BMI_0","baseline_DM","systolicBP_0","diastolicBP_0","ldl_0")
)
logistic_mca

results <- analyze_cox_model(
  data = data_ana_pheno_trem1_chip, 
  time_var = "surv_CVD", 
  event_var = "CVD", 
  group_var = "Triggering_receptor_expressed_on_myeloid_cells_1", 
  covariates = c("baseline_age","age2", "genetic_sex", "smoking_0", "PC1", "PC2", "PC3", "PCD4", "PC5","ldl_0","baseline_DM","BMI_0","baseline_HTN")
)
results

data_ana_pheno_trem1_tet2 <- filter(data_ana_pheno_trem1, data_ana_pheno_trem1$tet2_chip == 1)

dim(data_ana_pheno_trem1_tet2)

logistic_mca <- batch_logistic_regression(
  data = data_ana_pheno_trem1_tet2,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome = "`Coronary atherosclerosis [Atherosclerotic heart disease]`",
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0", "PC1","PC2","PC3","PCD4","PC5","BMI_0","baseline_DM","baseline_HTN","ldl_0")
)
logistic_mca

results <- analyze_cox_model(
  data = data_ana_pheno_trem1_tet2, 
  time_var = "surv_CVD0", 
  event_var = "CVD0", 
  group_var = "Triggering_receptor_expressed_on_myeloid_cells_1", 
  covariates = c("baseline_age","age2", "genetic_sex", "smoking_0", "PC1", "PC2", "PC3", "PCD4", "PC5","ldl_0","baseline_DM","BMI_0","baseline_HTN")
)
results

data_ana_pheno_trem1$trem1_333 <- 2
data_ana_pheno_trem1$trem1_333[data_ana_pheno_trem1$Triggering_receptor_expressed_on_myeloid_cells_1 >= quantile(data_ana_pheno_trem1$Triggering_receptor_expressed_on_myeloid_cells_1,2/3)] <- 3
data_ana_pheno_trem1$trem1_333[data_ana_pheno_trem1$Triggering_receptor_expressed_on_myeloid_cells_1 < quantile(data_ana_pheno_trem1$Triggering_receptor_expressed_on_myeloid_cells_1,1/3)] <- 1
data_ana_pheno_trem1$trem1_333<-as.factor(data_ana_pheno_trem1$trem1_333)
table(data_ana_pheno_trem1$trem1_333)

table(data_ana_pheno_trem1$trem1_333, data_ana_pheno_trem1$tet2_chip)

data_ana_pheno_trem1$trem1_262 <- 2
data_ana_pheno_trem1$trem1_262[data_ana_pheno_trem1$Triggering_receptor_expressed_on_myeloid_cells_1 >= quantile(data_ana_pheno_trem1$Triggering_receptor_expressed_on_myeloid_cells_1,0.8)] <- 3
data_ana_pheno_trem1$trem1_262[data_ana_pheno_trem1$Triggering_receptor_expressed_on_myeloid_cells_1 < quantile(data_ana_pheno_trem1$Triggering_receptor_expressed_on_myeloid_cells_1,0.2)] <- 1
data_ana_pheno_trem1$trem1_262 <- as.factor(data_ana_pheno_trem1$trem1_262)
table(data_ana_pheno_trem1$trem1_262)

data_ana_pheno_trem1_low <- filter(data_ana_pheno_trem1, data_ana_pheno_trem1$trem1_333 == 1)
data_ana_pheno_trem1_inter <- filter(data_ana_pheno_trem1, data_ana_pheno_trem1$trem1_333 == 2)
data_ana_pheno_trem1_high <- filter(data_ana_pheno_trem1, data_ana_pheno_trem1$trem1_333 == 3)

results <- analyze_cox_model(
  data = data_ana_pheno_trem1_low, 
  time_var = "surv_CVD0", 
  event_var = "CVD0", 
  group_var = "tet2_chip", 
  covariates = c("baseline_age","age2", "genetic_sex", "smoking_0", "PC1", "PC2", "PC3", "PCD4", "PC5","ldl_0","baseline_DM","BMI_0","baseline_HTN")
)

results <- analyze_cox_model(
  data = data_ana_pheno_trem1_inter, 
  time_var = "surv_CVD0", 
  event_var = "CVD0", 
  group_var = "tet2_chip", 
  covariates = c("baseline_age","age2", "genetic_sex", "smoking_0", "PC1", "PC2", "PC3", "PCD4", "PC5","ldl_0","baseline_DM","BMI_0","baseline_HTN")
)

results <- analyze_cox_model(
  data = data_ana_pheno_trem1_high, 
  time_var = "surv_CVD0", 
  event_var = "CVD0", 
  group_var = "tet2_chip", 
  covariates = c("baseline_age","age2", "genetic_sex", "smoking_0", "PC1", "PC2", "PC3", "PCD4", "PC5","ldl_0","baseline_DM","baseline_HTN")
)

data_all_pheno_filter <- filter(data_ana_pheno_trem1, data_ana_pheno_trem1$surv_CVD0 > 0)
table(data_all_pheno_filter$CVD0, data_all_pheno_filter$tet2_chip,data_all_pheno_filter$trem1_33)

data_ana_pheno_trem1_low <- filter(data_ana_pheno_trem1, data_ana_pheno_trem1$trem1_262 == 1)
data_ana_pheno_trem1_inter <- filter(data_ana_pheno_trem1, data_ana_pheno_trem1$trem1_262 == 2)
data_ana_pheno_trem1_high <- filter(data_ana_pheno_trem1, data_ana_pheno_trem1$trem1_262 == 3)

data_filter <- filter(data_all_pheno_filter,data_all_pheno_filter$surv_CVD0 >0)

cox_fit <- coxph(Surv(surv_CVD0, CVD0)~ tet2_chip*Triggering_receptor_expressed_on_myeloid_cells_1 + baseline_age+age2+genetic_sex+smoking_0+PC1+PC2+PC3+PCD4+PC5+ baseline_DM+ BMI_0+baseline_HTN, data = data_filter)
summary(cox_fit)

tet2_af <- tet2_chip_ukb[,c("person_id","AF")]

tet2_af$AF <- as.numeric(tet2_af$AF)
summary(tet2_af$AF)

tet2_af$tet2_af333 <- 2
tet2_af$tet2_af333[tet2_af$AF >= quantile(tet2_af$AF,2/3)] <- 3
tet2_af$tet2_af333[tet2_af$AF < quantile(tet2_af$AF,1/3)] <- 1
tet2_af$tet2_af333 <- as.factor(tet2_af$tet2_af333)
table(tet2_af$tet2_af333)

tet2_af$person_id <- as.character(tet2_af$person_id)
data_filter <- merge(data_filter,tet2_af, by.x = "ID_VUMC" ,by.y = "person_id",all.x = T)

data_all_pheno_tet2 <- filter(data_filter,data_filter$tet2_chip == 1 | data_filter$chip == 0)

dim(data_all_pheno_tet2)

data_all_pheno_tet2$tet2_af333 <- ifelse(is.na(data_all_pheno_tet2$tet2_af333),0, data_all_pheno_tet2$tet2_af333)
table(data_all_pheno_tet2$tet2_af333)

data_all_pheno_tet2$tet2_af333 <- as.numeric(data_all_pheno_tet2$tet2_af333)

cox_fit <- coxph(Surv(surv_CVD0, CVD0)~ tet2_af333*Triggering_receptor_expressed_on_myeloid_cells_1 + baseline_age+age2+genetic_sex+smoking_0+PC1+PC2+PC3+PCD4+PC5+ baseline_DM+ BMI_0+baseline_HTN, data = data_all_pheno_tet2)
summary(cox_fit)

x <- load("UKB_data_tet2_specific_0320.Rdata")
x

data_all_pheno_cox2$ID_VUMC <- as.character(data_all_pheno_cox2$ID_VUMC)
data_trem1_sep_tet2 <- merge(data_all_pheno_cox2, trem1, by.x = "ID_VUMC", by.y = "ID", all = F)
dim(data_trem1_sep_tet2)

data_trem1_sep_tet2$CVD0 <- NA
data_trem1_sep_tet2$CVD0 <- ifelse(data_trem1_sep_tet2$`Coronary atherosclerosis [Atherosclerotic heart disease]`==1|
                        data_trem1_sep_tet2$`Myocardial infarction [Heart attack]`==1 |
                        data_trem1_sep_tet2$`Heart failure`==1|
                        data_trem1_sep_tet2$`Angina pectoris`==1|
                        data_trem1_sep_tet2$`Peripheral vascular disease`==1 |
                        data_trem1_sep_tet2$Stroke==1,1, 0)
table(data_trem1_sep_tet2$CVD0)

data_trem1_sep_tet2 <- data_trem1_sep_tet2 %>% mutate(surv_CVD0 = pmin(surv_PVD,surv_AP,surv_hf,surv_str,surv_AMI,surv_AHD))

data_all_pheno_filter <- filter(data_trem1_sep_tet2, data_trem1_sep_tet2$surv_CVD0 > 0)

cox_fit <- coxph(Surv(surv_CVD0, CVD0)~ tet2_miss6*Triggering_receptor_expressed_on_myeloid_cells_1 + baseline_age+age2+genetic_sex+smoking_0+PC1+PC2+PC3+PCD4+PC5+ baseline_DM+ BMI_0+baseline_HTN, data = data_all_pheno_filter)
summary(cox_fit)

data_KM_trem1 <- data_ana_pheno_trem1_tet2[,c("ID_VUMC", "tet2_chip", "trem1_333","CVD0","surv_CVD0")]

save(data_KM_trem1, file= "data_KM_trem1.Rdata")

score_yp <- fread("ukb_52k_proteome_dec12_yp.tsv")

head(score_yp)
score_yp$eid <- as.character(score_yp$eid)

data_trem1_tet2 <- merge(data_ana_pheno_trem1_tet2, score_yp, by.x = "ID_VUMC", by.y = "eid", all = F)

dim(data_trem1_tet2)
table(data_trem1_tet2$CVD)

pro_list <- colnames(score_yp)[3:1465]

logistic_mca <- batch_logistic_regression(
  data = data_trem1_tet2,
  snp_list = pro_list,
  outcome = "CVD",
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0", "PC1","PC2","PC3","PCD4","PC5","BMI_0","baseline_DM","systolicBP_0","diastolicBP_0","ldl_0")
)

logistic_mca <- logistic_mca %>% arrange(P_Value)
logistic_mca

ukb_measure <- read.csv("ukb_measurement_data.csv")

head(ukb_measure)

names(data_all_pheno_cox2)[1] <- "person_id"

ukb_measure_all <- ukb_measure %>%
  left_join(data_all_pheno_cox2 %>% select(person_id, min_date), by = "person_id")

ukb_measure_all$measurement_datetime <- as.Date(ukb_measure_all$measurement_datetime)
head(ukb_measure_all)

dim(ukb_measure_all)

summary_stats <- ukb_measure_all %>%
  arrange(person_id, measurement_datetime) %>%
  group_by(person_id) %>%
  summarise(
    hgb_baseline = first(hgb),
    hgb_latest   = last(hgb),
    hgb_mean     = mean(hgb, na.rm = TRUE),
    hgb_dna      = hgb[which.min(abs(difftime(measurement_datetime, min_date, units = "days")))],

    mcv_baseline = first(mcv),
    mcv_latest   = last(mcv),
    mcv_mean     = mean(mcv, na.rm = TRUE),
    mcv_dna      = mcv[which.min(abs(difftime(measurement_datetime, min_date, units = "days")))],

    plt_baseline = first(plt),
    plt_latest   = last(plt),
    plt_mean     = mean(plt, na.rm = TRUE),
    plt_dna      = plt[which.min(abs(difftime(measurement_datetime, min_date, units = "days")))],

    rdw_baseline = first(rdw),
    rdw_latest   = last(rdw),
    rdw_mean     = mean(rdw, na.rm = TRUE),
    rdw_dna      = rdw[which.min(abs(difftime(measurement_datetime, min_date, units = "days")))],

    wbc_baseline = first(wbc),
    wbc_latest   = last(wbc),
    wbc_mean     = mean(wbc, na.rm = TRUE),
    wbc_dna      = wbc[which.min(abs(difftime(measurement_datetime, min_date, units = "days")))],

    .groups = "drop"
  )

dim(summary_stats)
head(summary_stats)

write.table(summary_stats, file = "Lab_data_ukb_KZ.txt")

dim(data_ana_pheno_trem1)
dim(data_nonchip_tet2chip)

table(data_nonchip_tet2chip$tet2_chip, useNA = "always")

summary_stats$person_id <- as.character(summary_stats$person_id)
data_all_tet2_heme <- merge(data_nonchip_tet2chip,summary_stats,by.x = "ID_VUMC", by.y = "person_id", all = F)

dim(data_all_tet2_heme)

blood_vars <- c("wbc_baseline", "wbc_latest", "wbc_dna", "wbc_mean",
                "hgb_baseline", "hgb_latest", "hgb_dna", "hgb_mean",
                "mcv_baseline", "mcv_latest", "mcv_dna", "mcv_mean",
                "plt_baseline", "plt_latest", "plt_dna", "plt_mean",
                "rdw_baseline", "rdw_latest", "rdw_dna", "rdw_mean")

analyze_tet2_assoc <- function(var, data) {
  formula <- as.formula(paste("tet2_chip ~", var, "+baseline_age+age2+genetic_sex+smoking_0+PC1+PC2+PC3+PCD4+PC5"))
  model <- glm(formula, data = data, family = binomial)

  tidy_result <- tidy(model) %>%
    filter(term == var) %>%
    mutate(
      OR = exp(estimate),
      CI_lower = exp(estimate - 1.96 * std.error),
      CI_upper = exp(estimate + 1.96 * std.error)
    ) %>%
    select(term, OR, CI_lower, CI_upper, p.value) %>%
    rename(Variable = term, P = p.value)

  return(tidy_result)
}

results <- lapply(blood_vars, analyze_tet2_assoc, data = data_all_tet2_heme)
results_df <- bind_rows(results)

results_df

tet2_heme <- filter(data_all_tet2_heme, data_all_tet2_heme$tet2_chip == 1)
dim(tet2_heme)

analyze_trem1_assoc <- function(var, data) {
  formula <- as.formula(paste("Triggering_receptor_expressed_on_myeloid_cells_1 ~", var, 
                              "+ baseline_age + age2 + genetic_sex + smoking_0 + PC1 + PC2 + PC3 + PCD4 + PC5"))
  
  model <- lm(formula, data = data)

  tidy_result <- tidy(model) %>%
    filter(term == var) %>%
    mutate(
      CI_lower = estimate - 1.96 * std.error,
      CI_upper = estimate + 1.96 * std.error
    ) %>%
    select(term, estimate, CI_lower, CI_upper, p.value) %>%
    rename(Variable = term, Beta = estimate, P = p.value)

  return(tidy_result)
}

results_tet2 <- lapply(blood_vars, analyze_trem1_assoc, data = tet2_heme)
results_tet2 <- bind_rows(results_tet2)

results_tet2

library(dplyr)
library(broom)

analyze_interaction <- function(var, data) {
  formula_str <- paste(var, "~ Triggering_receptor_expressed_on_myeloid_cells_1 * tet2_chip + baseline_age + age2 + genetic_sex + smoking_0 + PC1 + PC2 + PC3 + PCD4 + PC5")
  model <- lm(as.formula(formula_str), data = data)

  tidy_result <- tidy(model) %>%
    filter(term == "Triggering_receptor_expressed_on_myeloid_cells_1:tet2_chip") %>%
    mutate(
      CI_lower = estimate - 1.96 * std.error,
      CI_upper = estimate + 1.96 * std.error
    ) %>%
    select(term, estimate, CI_lower, CI_upper, p.value) %>%
    rename(
      Variable = term,
      Beta = estimate,
      P = p.value
    ) %>%
    mutate(Blood_Measure = var) %>%
    select(Blood_Measure, everything())

  return(tidy_result)
}

interaction_results <- lapply(blood_vars, analyze_interaction, data = data_all_tet2_heme)
interaction_summary <- bind_rows(interaction_results)

interaction_summary$P_Bonferroni <- p.adjust(interaction_summary$P, method = "fdr")
interaction_summary

library(metafor)

meta_hr_metamodel <- function(hr1, lower1, upper1, hr2, lower2, upper2,
                              method = "FE"  
) {
  yi  <- log(c(hr1, hr2))
  sei <- (log(c(upper1, upper2)) - log(c(lower1, lower2))) / (2 * 1.96)
  
  res <- rma(yi = yi, sei = sei, method = method)
  
  data.frame(
    HR       = as.numeric(exp(res$beta)),
    CI_lower = as.numeric(exp(res$ci.lb)),
    CI_upper = as.numeric(exp(res$ci.ub)),
    p_value  = as.numeric(res$pval)
  )
}

meta_hr_metamodel(
  hr1 = 0.85, lower1 = 0.47, upper1 = 1.45,
  hr2 = 1.05, lower2 = 0.86, upper2 = 1.28,
  method = "REML"
)

cytopenia_ukb <- read.csv("cytopenia_ukb_cohort.csv")

colnames(cytopenia_ukb)
dim(cytopenia_ukb)

table(cytopenia_ukb$gene, useNA = "always")

cytopenia_ukb

cytopenia_ukb2 <- cytopenia_ukb[,c("person_id","persistent_cytopenia","persistent_cytopenia_datetime","last_cbc_datetime",
                                   "persistent_anemia", "persistent_leukopenia", "persistent_thrombocytopenia",
                                   "index_datetime","gender", "race", "ethnicity","ever_smoker", "age")]

anyDuplicated(cytopenia_ukb2)

cytopenia_ukb3 <- distinct(cytopenia_ukb2)

dim(cytopenia_ukb3)

anyDuplicated(cytopenia_ukb3$person_id)

data_ccus <- merge(cytopenia_ukb3, trem1, by.x = "person_id",by.y = "ID", all = F)

data_ccus <- data_ccus %>%
  mutate(
    persistent_cytopenia_datetime = as.Date(ymd_hms(persistent_cytopenia_datetime)),
    last_cbc_datetime = as.Date(ymd_hms(last_cbc_datetime)),
    index_datetime = as.Date(ymd_hms(index_datetime))
  )

data_cova <- data_all_pheno_cox2[,c("ID_VUMC","baseline_age", "age2", "genetic_sex", "smoking_0",
                                         "PC1","PC2","PC3","PCD4","PC5","BMI_0","baseline_HTN","baseline_DM",
                                         "systolicBP_0","diastolicBP_0","ldl_0","chip","tet2_chip")]

data_ccus2 <- merge(data_ccus, data_cova, by.x = "person_id", by.y = "ID_VUMC", all = F)

dim(data_ccus)
dim(data_ccus2)

library(survival)
library(dplyr)
library(rlang)
library(broom)

run_cox_model_each_var <- function(data,
                                   event_col,
                                   event_time_col,
                                   baseline_time_col,
                                   last_record_time_col,
                                   covariates,
                                   variables_of_interest) {
  
  event_sym <- sym(event_col)
  event_time_sym <- sym(event_time_col)
  baseline_sym <- sym(baseline_time_col)
  last_sym <- sym(last_record_time_col)
  
  df <- data %>%
    mutate(
      followup_time = as.numeric(pmin(!!event_time_sym, !!last_sym, na.rm = TRUE) - !!baseline_sym),
      status = !!event_sym
    ) %>%
    filter(!is.na(followup_time), followup_time > 0)
  
  total_event <- sum(df$status == 1, na.rm = TRUE)

  surv_obj <- Surv(time = df$followup_time, event = df$status)

  results_list <- list()

  for (var in variables_of_interest) {
    var_sym <- sym(var)

    N <- df %>% filter(!!var_sym == 1) %>% nrow()
    N_event <- df %>% filter(!!var_sym == 1, status == 1) %>% nrow()
    
    event_rate <- if (total_event > 0) N_event / total_event else NA_real_

    formula_str <- paste("surv_obj ~", paste(c(covariates, var), collapse = " + "))
    cox_formula <- as.formula(formula_str)

    model <- coxph(cox_formula, data = df)

    result <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == var) %>%
      mutate(
        variable = var,
        N = N,
        N_event = N_event,
        event_rate = event_rate,
        HR = estimate,
        CI_lower = conf.low,
        CI_upper = conf.high,
        p_value = p.value
      ) %>%
      select(variable, N, N_event, event_rate, HR, CI_lower, CI_upper, p_value)
    
    results_list[[var]] <- result
  }

  final_result <- bind_rows(results_list)

  return(final_result)
}

ls(data_ccus2)

cox_result <- run_cox_model_each_var(
  data = data_ccus2,
  event_col = "persistent_cytopenia",
  event_time_col = "persistent_cytopenia_datetime",
  baseline_time_col = "index_datetime",
  last_record_time_col = "last_cbc_datetime",
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0","PC1","PC2","PC3","PCD4","PC5"),
  variables_of_interest = "Triggering_receptor_expressed_on_myeloid_cells_1"
)

cox_result

data_ccus_chip <- filter(data_ccus2, data_ccus2$chip == 1)

cox_result <- run_cox_model_each_var(
  data = data_ccus_chip,
  event_col = "persistent_cytopenia",
  event_time_col = "persistent_cytopenia_datetime",
  baseline_time_col = "index_datetime",
  last_record_time_col = "last_cbc_datetime",
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0","PC1","PC2","PC3","PCD4","PC5"),
  variables_of_interest = "Triggering_receptor_expressed_on_myeloid_cells_1"
)

cox_result

data_ccus_tet2 <- filter(data_ccus2, data_ccus2$tet2_chip == 1)

cox_result <- run_cox_model_each_var(
  data = data_ccus_tet2,
  event_col = "persistent_cytopenia",
  event_time_col = "persistent_cytopenia_datetime",
  baseline_time_col = "index_datetime",
  last_record_time_col = "last_cbc_datetime",
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0","PC1","PC2","PC3","PCD4","PC5"),
  variables_of_interest = "Triggering_receptor_expressed_on_myeloid_cells_1"
)

cox_result

tet2_chip_ukb_large <- filter(tet2_chip_ukb, tet2_chip_ukb$AF >=0.1)

data_ccus_tet2_large <- filter(data_ccus2, data_ccus2$person_id %in% tet2_chip_ukb_large$person_id)
dim(data_ccus_tet2_large)

cox_result <- run_cox_model_each_var(
  data = data_ccus_tet2_large,
  event_col = "persistent_cytopenia",
  event_time_col = "persistent_cytopenia_datetime",
  baseline_time_col = "index_datetime",
  last_record_time_col = "last_cbc_datetime",
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0","PC1","PC2","PC3","PCD4","PC5"),
  variables_of_interest = "Triggering_receptor_expressed_on_myeloid_cells_1"
)

cox_result

outcome_list <- c("persistent_cytopenia", "persistent_anemia", "persistent_leukopenia", "persistent_thrombocytopenia")

result_df <- batch_logistic_by_outcome(
  data = data_ccus2,
  snp = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome_list = outcome_list,
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0", "PC1","PC2","PC3","PCD4","PC5")
)
result_df

result_df <- batch_logistic_by_outcome(
  data = data_ccus_chip,
  snp = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome_list = outcome_list,
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0", "PC1","PC2","PC3","PCD4","PC5")
)
result_df

result_df <- batch_logistic_by_outcome(
  data = data_ccus_tet2,
  snp = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome_list = outcome_list,
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0", "PC1","PC2","PC3","PCD4","PC5","BMI_0")
)
result_df

result_df <- batch_logistic_by_outcome(
  data = data_ccus_tet2_large,
  snp = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome_list = outcome_list,
  covariates = c("baseline_age", "age2", "genetic_sex", "smoking_0", "PC1","PC2","PC3","PCD4","PC5", "BMI_0")
)
result_df



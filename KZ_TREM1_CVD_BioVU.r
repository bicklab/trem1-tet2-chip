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

data_final3 <- readRDS("pheno_data_biovu_250k_KZ.rds")

score_kz_rint <- readRDS("omicsscore_biovu_rint_KZ.rds")

chip_call_biovu_250k <- fread("agd_chip_calls_250k__2025_05_06.tsv")
head(chip_call_biovu_250k)

tet2 <- filter(chip_call_biovu_250k, chip_call_biovu_250k$Gene.refGene == "TET2")
d3a <- filter(chip_call_biovu_250k, chip_call_biovu_250k$Gene.refGene == "DNMT3A")

trem1 <- score_kz_rint[,c("ID", "Triggering_receptor_expressed_on_myeloid_cells_1")]

data_ana_trem1 <- data_final3[,c("person_id","GRID","CVD","surv_CVD","age","age2","gender","race","t2d","hypertension","ldl_median",
                              "ASCVD","surv_ASCVD","CAD","surv_CAD","surv_stroke","stroke","heartfailure","surv_hf",
                              "PVD","surv_PVD","MI","surv_MI","CV_431.1","CV_431.11","CVD1","surv_CVD1","baseline_HTN","baseline_T2D")]

data_ana_trem1$chip <- ifelse(data_ana_trem1$GRID %in% chip_call_biovu_250k$GRID, 1, 0)
data_ana_trem1$tet2 <- ifelse(data_ana_trem1$GRID %in% tet2$GRID, 1, 0)
data_ana_trem1$d3a <- ifelse(data_ana_trem1$GRID %in% d3a$GRID, 1, 0)

#EUR1

ancestry_250k_biovu <- fread("biovu_prs_scores/EUR_BCF_sum_scores_0.001.txt")
dim(ancestry_250k_biovu)

#eur_250k <- filter(ancestry_250k_biovu, ancestry_250k_biovu$supervised_ancestry_cluster == "EUR")

#Ancestry Yash
data_ana_trem1_eur <- filter(data_ana_trem1, data_ana_trem1$GRID %in% eur_250k$GRID)
dim(data_ana_trem1_eur)

#Ancestry HannahP
data_ana_trem1_eur <- filter(data_ana_trem1, data_ana_trem1$GRID %in% ancestry_250k_biovu$FID)
dim(data_ana_trem1_eur)

data_ana_trem1_eur <- merge(data_ana_trem1_eur, trem1, by.x = "GRID", by.y = "ID", all = F)

dim(data_ana_trem1_eur)
summary(data_ana_trem1_eur$Triggering_receptor_expressed_on_myeloid_cells_1)

data_ana_trem1_eur_4090 <- filter(data_ana_trem1_eur, data_ana_trem1_eur$age >=40 & data_ana_trem1_eur$age <=90)
dim(data_ana_trem1_eur_4090)

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

batch_cox_regression <- function(data, snp_list, surv_time, outcome, covariates = NULL) {
  results <- list()

  surv_obj <- Surv(time = data[[surv_time]], event = data[[outcome]])

  for (snp in snp_list) {
    formula_str <- paste("surv_obj ~", snp)
    if (!is.null(covariates) && length(covariates) > 0) {
      formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
    }
    formula <- as.formula(formula_str)

    tryCatch({
      data_for_model <- data
      data_for_model$surv_obj <- surv_obj

      cox_model <- coxph(formula, data = data_for_model)
      cox_summary <- summary(cox_model)

      coef <- cox_summary$coefficients[snp, "coef"]
      se <- cox_summary$coefficients[snp, "se(coef)"]
      p_value <- cox_summary$coefficients[snp, "Pr(>|z|)"]
      hr <- exp(coef)
      ci_lower <- exp(coef - 1.96 * se)
      ci_upper <- exp(coef + 1.96 * se)

      results[[snp]] <- data.frame(
        SNP = snp,
        Beta = coef,
        HR = hr,
        SE = se,
        CI_Lower = ci_lower,
        CI_Upper = ci_upper,
        P_Value = p_value
      )
    }, error = function(e) {
      warning(paste("Error in model for SNP:", snp, "-", e$message))
    })
  }

  results_df <- do.call(rbind, results)
  return(results_df)
}

data_filter <- filter(data_ana_trem1_eur_4090, data_ana_trem1_eur_4090$surv_CVD1 >0)

cox_results <- batch_cox_regression(
  data = data_filter,
  snp_list = "tet2",
  surv_time = "surv_CVD1",
  outcome = "CVD1",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

data_filter <- filter(data_ana_trem1_eur_4090, data_ana_trem1_eur_4090$surv_CVD1 >0)

cox_results <- batch_cox_regression(
  data = data_filter,
  snp_list = "chip",
  surv_time = "surv_CVD1",
  outcome = "CVD1",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

logistic_mca <- batch_logistic_regression(
  data = data_ana_trem1_eur_4090,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
logistic_mca

data_filter <- filter(data_ana_trem1_eur_4090, data_ana_trem1_eur_4090$surv_CVD >0)

cox_results <- batch_cox_regression(
  data = data_filter,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_CVD",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

data_ana_trem1_eur$tet2 <- NA
data_ana_trem1_eur$tet2 <- ifelse(data_ana_trem1_eur$GRID %in% tet2$GRID, 1, 0 )
table(data_ana_trem1_eur$tet2)

data_tet2_trem1_eur <- filter(data_ana_trem1_eur, data_ana_trem1_eur$tet2 == 1)
dim(data_tet2_trem1_eur)

data_tet2_trem1_eur_4090 <- filter(data_tet2_trem1_eur, data_tet2_trem1_eur$age >=40 & data_tet2_trem1_eur$age <=90)
dim(data_tet2_trem1_eur_4090)

logistic_mca <- batch_logistic_regression(
  data = data_tet2_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
logistic_mca

data_filter <- filter(data_tet2_trem1_eur, data_tet2_trem1_eur$surv_CVD >0)

cox_results <- batch_cox_regression(
  data = data_tet2_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_CVD",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

CVDlist <- c("CVD","CVD1","ASCVD","CAD","MI", "heartfailure","stroke","PVD","CV_431.1","CV_431.11")

results_list <- list()

for (cvd in CVDlist) {
  cat("Analyzing:", cvd, "\n")
  
  formula <- as.formula(paste0(
    "`", cvd, "` ~ Triggering_receptor_expressed_on_myeloid_cells_1 + age + age2 + gender + ",
    "t2d + baseline_HTN + ldl_median"
  ))
  
  model <- glm(formula, data = data_tet2_trem1_eur, family = binomial)
  
  trem1_row <- summary(model)$coefficients["Triggering_receptor_expressed_on_myeloid_cells_1", ]
  
  results_list[[cvd]] <- data.frame(
    CVD = cvd,
    beta = trem1_row["Estimate"],
    std_error = trem1_row["Std. Error"],
    p_value = trem1_row["Pr(>|z|)"]
  )
}

cvd_results <- do.call(rbind, results_list)
cvd_results

data_filter <- filter(data_tet2_trem1_eur, data_tet2_trem1_eur$surv_CVD >0)
cox_results <- batch_cox_regression(
  data = data_tet2_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_CVD",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

data_filter <- filter(data_tet2_trem1_eur, data_tet2_trem1_eur$surv_CAD >0)
cox_results <- batch_cox_regression(
  data = data_tet2_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_CAD",
  outcome = "CAD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

data_filter <- filter(data_tet2_trem1_eur, data_tet2_trem1_eur$surv_MI >0)
cox_results <- batch_cox_regression(
  data = data_tet2_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_MI",
  outcome = "MI",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

data_filter <- filter(data_tet2_trem1_eur, data_tet2_trem1_eur$surv_stroke >0)
cox_results <- batch_cox_regression(
  data = data_tet2_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_stroke",
  outcome = "CV_431.11",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

data_filter <- filter(data_tet2_trem1_eur, data_tet2_trem1_eur$surv_hf >0)
cox_results <- batch_cox_regression(
  data = data_tet2_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_hf",
  outcome = "heartfailure",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

data_filter <- filter(data_tet2_trem1_eur, data_tet2_trem1_eur$surv_PVD >0)
cox_results <- batch_cox_regression(
  data = data_tet2_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_PVD",
  outcome = "PVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results



data_ana_trem1_eur$chip <- NA
data_ana_trem1_eur$chip <- ifelse(data_ana_trem1_eur$GRID %in% chip_call_biovu_250k$GRID, 1, 0 )
table(data_ana_trem1_eur$chip)

109766+4522

4522/(4522+109766)

dim(data_ana_trem1_eur)
ls(data_ana_trem1_eur)

data_chip_trem1_eur <- filter(data_ana_trem1_eur, data_ana_trem1_eur$chip == 1)
dim(data_chip_trem1_eur)

logistic_mca <- batch_logistic_regression(
  data = data_chip_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
logistic_mca

data_filter <- filter(data_chip_trem1_eur, data_chip_trem1_eur$surv_CVD >0)

cox_results <- batch_cox_regression(
  data = data_filter,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_CVD",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

data_d3a_trem1_eur <- filter(data_ana_trem1_eur, data_ana_trem1_eur$d3a == 1)
dim(data_d3a_trem1_eur)

logistic_mca <- batch_logistic_regression(
  data = data_d3a_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
logistic_mca

data_filter <- filter(data_d3a_trem1_eur, data_d3a_trem1_eur$surv_CVD >0)

cox_results <- batch_cox_regression(
  data = data_filter,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_CVD",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

data_nonchip_trem1_eur <- filter(data_ana_trem1_eur, data_ana_trem1_eur$chip == 0)
dim(data_nonchip_trem1_eur)

logistic_mca <- batch_logistic_regression(
  data = data_nonchip_trem1_eur,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
logistic_mca

data_filter <- filter(data_nonchip_trem1_eur, data_nonchip_trem1_eur$surv_CVD >0)

cox_results <- batch_cox_regression(
  data = data_filter,
  snp_list = "Triggering_receptor_expressed_on_myeloid_cells_1",
  surv_time = "surv_CVD",
  outcome = "CVD",
  covariates = c("age", "age2", "gender","baseline_HTN","t2d","ldl_median")
)
cox_results

data_nonchip_tet2chip <- filter(data_ana_trem1_eur, data_ana_trem1_eur$chip == 0 | data_ana_trem1_eur$tet2 == 1)
dim(data_nonchip_tet2chip)

cox_fit <- coxph(Surv(surv_CVD1, CVD1)~ tet2*Triggering_receptor_expressed_on_myeloid_cells_1 + age + age2 + gender + t2d+ baseline_HTN + ldl_median, data = data_nonchip_tet2chip)
summary(cox_fit)

tet2_af <- tet2[,c("GRID","AF")]

tet2_af$AF <- as.numeric(tet2_af$AF)
summary(tet2_af$AF)

tet2_af$tet2_af333 <- 2
tet2_af$tet2_af333[tet2_af$AF >= quantile(tet2_af$AF,2/3)] <- 3
tet2_af$tet2_af333[tet2_af$AF < quantile(tet2_af$AF,1/3)] <- 1
tet2_af$tet2_af333 <- as.factor(tet2_af$tet2_af333)
table(tet2_af$tet2_af333)

tet2_af$GRID <- as.character(tet2_af$GRID)
data_ana_trem1_eur <- merge(data_ana_trem1_eur,tet2_af, by = "GRID",all.x = T)

data_ana_trem1_eur$tet2_af333 <- ifelse(is.na(data_ana_trem1_eur$tet2_af333),0, data_ana_trem1_eur$tet2_af333)
table(data_ana_trem1_eur$tet2_af333)

data_ana_trem1_eur$tet2_af333 <- as.numeric(data_ana_trem1_eur$tet2_af333)

data_nonchip_tet2chip <- filter(data_ana_trem1_eur, data_ana_trem1_eur$chip == 0 | data_ana_trem1_eur$tet2 == 1)
dim(data_nonchip_tet2chip)

cox_fit <- coxph(Surv(surv_CVD1, CVD1)~ tet2_af333*Triggering_receptor_expressed_on_myeloid_cells_1 + age + age2+ race + gender + t2d+ hypertension + ldl_median, data = data_nonchip_tet2chip)
summary(cox_fit)

data_ana_trem1_eur$trem1_333 <- 2
data_ana_trem1_eur$trem1_333[data_ana_trem1_eur$Triggering_receptor_expressed_on_myeloid_cells_1 >= quantile(data_ana_trem1_eur$Triggering_receptor_expressed_on_myeloid_cells_1,2/3)] <- 3
data_ana_trem1_eur$trem1_333[data_ana_trem1_eur$Triggering_receptor_expressed_on_myeloid_cells_1 < quantile(data_ana_trem1_eur$Triggering_receptor_expressed_on_myeloid_cells_1,1/3)] <- 1
data_ana_trem1_eur$trem1_333<-as.factor(data_ana_trem1_eur$trem1_333)
table(data_ana_trem1_eur$trem1_333)

data_nonchip_tet2chip <- filter(data_ana_trem1_eur, data_ana_trem1_eur$chip == 0 | data_ana_trem1_eur$tet2 == 1)
dim(data_nonchip_tet2chip)

data_filter <- filter(data_nonchip_tet2chip, data_nonchip_tet2chip$age >40 & data_nonchip_tet2chip$age <90)
data_filter <- filter(data_filter, data_filter$surv_CVD1 >0)

cox_results <- batch_cox_regression(
  data = data_filter,
  snp_list = "tet2",
  surv_time = "surv_CVD1",
  outcome = "CVD1",
  covariates = c("age", "age2", "gender","hypertension","t2d","ldl_median")
)
cox_results

data_filter_low <- filter(data_filter, data_filter$trem1_333 == 1)

cox_results <- batch_cox_regression(
  data = data_filter_low,
  snp_list = "tet2",
  surv_time = "surv_CVD1",
  outcome = "CVD1",
  covariates = c("age", "age2", "gender","hypertension","t2d","ldl_median")
)
cox_results

data_filter_inter <- filter(data_filter, data_filter$trem1_333 == 2)

cox_results <- batch_cox_regression(
  data = data_filter_inter,
  snp_list = "tet2",
  surv_time = "surv_CVD1",
  outcome = "CVD1",
  covariates = c("age", "age2", "gender","hypertension","t2d","ldl_median")
)
cox_results

data_filter_high <- filter(data_filter, data_filter$trem1_333 == 3)

cox_results <- batch_cox_regression(
  data = data_filter_high,
  snp_list = "tet2",
  surv_time = "surv_CVD1",
  outcome = "CVD1",
  covariates = c("age", "age2", "gender","hypertension","t2d","ldl_median")
)
cox_results

data_km_biovu_trem1 <- data_filter[,c("GRID","trem1_333","tet2","CVD","CVD1","surv_CVD","surv_CVD1")]

save(data_km_biovu_trem1, file = "data_km_biovu_trem1.Rdata")

tet2_chip_biovu <- filter(chip_call_biovu_250k, chip_call_biovu_250k$Gene.refGene == "TET2")

library(stringr)

tet2_chip_biovu <- tet2_chip_biovu %>%
  mutate(NonsynOI_numeric = str_extract(aachange, "\\d+"))

tet2_chip_biovu$NonsynOI_numeric <- as.numeric(tet2_chip_biovu$NonsynOI_numeric)

tet2_chip_biovu_fssg <- filter(tet2_chip_biovu, tet2_chip_biovu$ExonicFunc.refGene == "frameshift deletion" |
                            tet2_chip_biovu$ExonicFunc.refGene == "frameshift insertion" | tet2_chip_biovu$ExonicFunc.refGene == "stopgain")

summary(tet2_chip_biovu$NonsynOI_numeric)

tet2_chip_biovu_fssg$fssg_early <- NA
tet2_chip_biovu_fssg$fssg_early <- ifelse(tet2_chip_biovu_fssg$NonsynOI_numeric < 1129, 1, 0)
tet2_chip_biovu_fssg$fssg_late <- NA
tet2_chip_biovu_fssg$fssg_late <- ifelse(tet2_chip_biovu_fssg$NonsynOI_numeric >= 1129, 1, 0)

tet2_chip_biovu_fssg_early <- filter(tet2_chip_biovu_fssg, tet2_chip_biovu_fssg$fssg_early == 1)
tet2_chip_biovu_fssg_late <- filter(tet2_chip_biovu_fssg, tet2_chip_biovu_fssg$fssg_late == 1)

tet2_chip_biovu_miss <- filter(tet2_chip_biovu, tet2_chip_biovu$ExonicFunc.refGene == "nonsynonymous SNV")

tet2_miss1 <- filter(tet2_chip_biovu_miss, tet2_chip_biovu_miss$NonsynOI_numeric <= 1128)
tet2_miss2 <- filter(tet2_chip_biovu_miss, 1129 <= tet2_chip_biovu_miss$NonsynOI_numeric & tet2_chip_biovu_miss$NonsynOI_numeric <= 1312)
tet2_miss3 <- filter(tet2_chip_biovu_miss, 1313 <= tet2_chip_biovu_miss$NonsynOI_numeric & tet2_chip_biovu_miss$NonsynOI_numeric <= 1380)
tet2_miss4 <- filter(tet2_chip_biovu_miss, 1381 <= tet2_chip_biovu_miss$NonsynOI_numeric & tet2_chip_biovu_miss$NonsynOI_numeric <= 1750)
tet2_miss5 <- filter(tet2_chip_biovu_miss, 1751 <= tet2_chip_biovu_miss$NonsynOI_numeric & tet2_chip_biovu_miss$NonsynOI_numeric <= 1850)
tet2_miss6 <- filter(tet2_chip_biovu_miss, 1851 <= tet2_chip_biovu_miss$NonsynOI_numeric & tet2_chip_biovu_miss$NonsynOI_numeric <= 2002)
dim(tet2_miss1)
dim(tet2_miss2)
dim(tet2_miss3)
dim(tet2_miss4)
dim(tet2_miss5)
dim(tet2_miss6)

data_ana_trem1_eur$tet2 <- NA
data_ana_trem1_eur$tet2 <- ifelse(data_ana_trem1_eur$GRID %in% tet2$GRID, 1, 0)

data_ana_trem1_eur$tet2_fssg <- ifelse(data_ana_trem1_eur$GRID %in% tet2_chip_biovu_fssg$GRID, 1, 0)

data_ana_trem1_eur$tet2_fssg_early <- ifelse(data_ana_trem1_eur$GRID %in% tet2_chip_biovu_fssg_early$GRID, 1, 0)
data_ana_trem1_eur$tet2_fssg_late <- ifelse(data_ana_trem1_eur$GRID %in% tet2_chip_biovu_fssg_late$GRID, 1, 0)

data_ana_trem1_eur$tet2_miss <- ifelse(data_ana_trem1_eur$GRID %in% tet2_chip_biovu_miss$GRID, 1, 0)

data_ana_trem1_eur$tet2_miss2 <- ifelse(data_ana_trem1_eur$GRID %in% tet2_miss2$GRID, 1, 0)
data_ana_trem1_eur$tet2_miss3 <- ifelse(data_ana_trem1_eur$GRID %in% tet2_miss3$GRID, 1, 0)
data_ana_trem1_eur$tet2_miss4 <- ifelse(data_ana_trem1_eur$GRID %in% tet2_miss4$GRID, 1, 0)
data_ana_trem1_eur$tet2_miss6 <- ifelse(data_ana_trem1_eur$GRID %in% tet2_miss6$GRID, 1, 0)

data_nonchip_tet2chip <- filter(data_ana_trem1_eur, data_ana_trem1_eur$chip == 0 | data_ana_trem1_eur$tet2 == 1)
dim(data_nonchip_tet2chip)

data_filter <- filter(data_nonchip_tet2chip,data_nonchip_tet2chip$surv_CVD >0)

cox_fit <- coxph(Surv(surv_CVD1, CVD1)~ tet2*Triggering_receptor_expressed_on_myeloid_cells_1 + age + age2 + gender + t2d+ baseline_HTN + ldl_median, data = data_nonchip_tet2chip)
summary(cox_fit)

cox_fit <- coxph(Surv(surv_CVD1, CVD1)~ tet2_fssg*Triggering_receptor_expressed_on_myeloid_cells_1 + age + age2 + gender + t2d+ baseline_HTN + ldl_median, data = data_nonchip_tet2chip)
summary(cox_fit)

cox_fit <- coxph(Surv(surv_CVD1, CVD1)~ tet2_fssg*Triggering_receptor_expressed_on_myeloid_cells_1 + age + age2 + gender + t2d+ baseline_HTN + ldl_median, data = data_filter)
summary(cox_fit)

cox_fit <- coxph(Surv(surv_CVD1, CVD1)~ tet2_fssg_early*Triggering_receptor_expressed_on_myeloid_cells_1 + age + age2 + gender + t2d+ baseline_HTN + ldl_median, data = data_nonchip_tet2chip)
summary(cox_fit)

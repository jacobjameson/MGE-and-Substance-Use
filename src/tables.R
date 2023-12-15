#-------------------------------------------------------------------------
# AUTHOR:             Jacob Jameson
# PURPOSE:            Generate main table results
#-------------------------------------------------------------------------
#
# Create dataset ----------------------------------------------------------

source('src/Construct Analytical Dataset.R')

# load packages ----------------------------------------------------------
library(lmtest)
library(sandwich)
library(broom)
library(tidyverse)
library(sjstats)
library(MASS)
library(stargazer)
library(gtsummary)

#----------------------------------------------------------------------------------
##########################################################################
# TABLE 1 RESULTS
#
# Table 1: Demographic characteristics of self-identified males enrolled 
# in the National Longitudinal Study of Adolescent to Adult Health with 
# complete school network survey data and GE scores, stratified by
# gender expression (GE) status 
##########################################################################

t1.vars <- c("race", 'pseudo.gpa', 'sespc_al', 'nhood1_d', 'edu', 'insurance',
             "w1.cigarettes", "w4.cigarettes.bin.30",
             "w1.marijuana", "w4.marijuana.bin.year", "w4.marijuana.bin.30",
             "w1.drunk", "w4.drunk.bin.year", "w4.drunk.bin.30",
             "w1.recreational", "w4.fav.bin.year", "w4.fav.bin.30",
             "w4.prescription")

t1.unweighted <- analytical_dataset %>% 
  filter(in_sample == 1) %>%
  select(one_of(t1.vars), increasing) %>%
  tbl_summary(by = increasing,
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c("{mean} ({sd})"),
              missing_text = "(Missing)") %>% 
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2),
        test = list(all_continuous() ~ "aov",
                    all_categorical() ~ "chisq.test")) %>%
  add_overall() %>% modify_header(label ~ "**Variable**") %>%
  modify_caption("**Patient Characteristics (unweighted)**") %>%
  bold_labels()

gt::gtsave(as_gt(t1.unweighted), "outputs/tables/Table 1 Unweighted.pdf")

weighted <- survey::svydesign(id=~psuscid, strata=~region,
                            weights=~gswgt4_2, 
                            data=analytical_dataset, nest=TRUE)

t1.weighted <- subset(weighted, in_sample == 1) %>%
  tbl_svysummary(
    by = increasing, 
    type = all_continuous() ~ "continuous2",
    statistic = 
      all_continuous() ~ c("{mean} ({sd})"),
    missing = "always",
    include = c(one_of(t1.vars), increasing)) %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2),
        test = list(all_continuous() ~ "svy.kruskal.test",
                    all_categorical() ~ "svy.chisq.test")) %>%
  add_overall() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Table 1. Patient Characteristics (weighted)**") %>%
  bold_labels() 

gt::gtsave(as_gt(t1.weighted), "outputs/tables/Table 1 Weighted.pdf")

##########################################################################
# TABLE 2 REGRESSION RESULTS
#
# Table 2. Coefficients from logistic regressions models assessing 
# association of gender expression (GE) with adult substance use 
# behaviors over prior 30-days
##########################################################################

# Define the models
models <- list(
  w4.cigarettes.bin.30 = "w1.cigarettes",
  w4.marijuana.bin.30 = "w1.marijuana",
  w4.drunk.bin.30 = "w1.drunk",
  w4.fav.bin.30 = "w1.recreational",
  w4.prescription = NULL
)

# Define the list of predictors
predictors <- c("w1.GE_male_std", "w4.GE_male_std", "delta_w1_w4_GE")
base_predictors <- "race + pseudo.gpa + sespc_al + nhood1_d"

# Start the sink to capture output
sink("outputs/tables/Table 2.txt")

# Loop through each predictor
for (predictor in predictors) {
  
  cat(paste("\nResults for predictor:", predictor, "\n"))
  
  model_list <- list()
  robust_list <- list()
  
  # Loop through each model and run regression for current predictor
  for (var in names(models)) {
    
    # Construct the formula
    unique_predictor <- models[[var]]
    if (is.null(unique_predictor)) {
      formula_str <- paste(var, "~", base_predictors, "+", predictor)
    } else {
      formula_str <- paste(var, "~", base_predictors, "+", predictor, "+", unique_predictor)
    }
    formula_obj <- as.formula(formula_str)
    
    # Run regression
    mod <- glm(formula = formula_obj, data = subset(analytical_dataset, in_sample == 1),
               weights = weights, family = "quasibinomial")
    
    # Store the model and robust SEs
    model_list[[var]] <- mod
    robust_list[[var]] <- coeftest(x = mod, vcov = vcovCL(mod, type = "HC0", cluster = ~ cluster))
  }
  
  stargazer(robust_list, type = "text", report = "vc*p", header = FALSE, 
            title = paste("Table 2: Results for predictor:", predictor), 
            column.labels = names(models))  
}

sink()
##########################################################################
##########################################################################
# APPENDIX 4 REGRESSION RESULTS
#
# Appendix 4: Associations between adolescent school and social network
# variables and changes in gender expression (GE) between adolescence 
# and adulthood
##########################################################################

# Define the list of predictors
predictors <- c('w1.GE_male_std_school', 'katz_centrality_R', '
                nominations', 'num_bff_noms', 'mbff_reciprocity', 'fbff_reciprocity')

# Initialize a list to store the models and robust standard errors
model_list <- list()
robust_list <- list()

# Run regression for each predictor individually
for (predictor in predictors) {
  formula_str <- paste("delta_w1_w4_GE ~", predictor, "+ race + pseudo.gpa + sespc_al + nhood1_d")
  mod <- glm(formula_str, data = subset(analytical_dataset, in_sample == 1), weights = weights)
  robust_list[[predictor]] <- coeftest(x = mod, vcov = vcovCL(mod, type = "HC0", cluster = ~ cluster))
}

# Run regression with all predictors together
formula_all <- as.formula(paste("delta_w1_w4_GE ~", paste(paste(predictors, collapse = " + "), "+ race + pseudo.gpa + sespc_al + nhood1_d")))
mod_all <- glm(formula_all, data = subset(analytical_dataset, in_sample == 1), weights = weights)
robust_list[["All together"]] <- coeftest(x = mod_all, vcov = vcovCL(mod_all, type = "HC0", cluster = ~ cluster))

# Save the results to a .txt file
sink("outputs/tables/Appendix 4.txt")

stargazer(robust_list, type = "text", report = "vc*p", 
          header = FALSE, title = "Appendix 4")

sink()
##########################################################################
##########################################################################
# APPENDIX 5 REGRESSION RESULTS
#
# Appendix 5. Coefficients from logistic regression models assessing 
# association of gender expression (GE) with adult substance use behaviors 
# over prior 12-months
##########################################################################

# Define the models
models <- list(
  w4.marijuana.bin.year = "w1.marijuana",
  w4.drunk.bin.year = "w1.drunk",
  w4.fav.bin.year = "w1.recreational"
)

# Define the list of predictors
predictors <- c("w1.GE_male_std", "w4.GE_male_std", "delta_w1_w4_GE")
base_predictors <- "race + pseudo.gpa + sespc_al + nhood1_d"

# Start the sink to capture output
sink("outputs/tables/Appendix 5.txt")

# Loop through each predictor
for (predictor in predictors) {
  
  cat(paste("\nResults for predictor:", predictor, "\n"))
  
  model_list <- list()
  robust_list <- list()
  
  # Loop through each model and run regression for current predictor
  for (var in names(models)) {
    
    # Construct the formula
    unique_predictor <- models[[var]]
    if (is.null(unique_predictor)) {
      formula_str <- paste(var, "~", base_predictors, "+", predictor)
    } else {
      formula_str <- paste(var, "~", base_predictors, "+", predictor, "+", unique_predictor)
    }
    formula_obj <- as.formula(formula_str)
    
    # Run regression
    mod <- glm(formula = formula_obj, data = subset(analytical_dataset, in_sample == 1),
               weights = weights, family = "quasibinomial")
    
    # Store the model and robust SEs
    model_list[[var]] <- mod
    robust_list[[var]] <- coeftest(x = mod, vcov = vcovCL(mod, type = "HC0", cluster = ~ cluster))
  }
  
  stargazer(robust_list, type = "text", report = "vc*p", header = FALSE, 
            title = paste("Appendix 5: Results for predictor:", predictor), 
            column.labels = names(models))  
}

sink()

##########################################################################
# APPENDIX 6 REGRESSION RESULTS
#
# Appendix 6. Coefficients from logistic regression models 
# assessing impact of gender expression (GE) over time on adult 
# lifetime substance use
##########################################################################

# Define the models
models <- list(
  w4.cigarettes = "w1.cigarettes",
  w4.marijuana = "w1.marijuana",
  w4.drunk = "w1.drunk",
  w4.recreational = "w1.recreational",
  w4.prescription = NULL
)

# Define the list of predictors
predictors <- c("w1.GE_male_std", "w4.GE_male_std", "delta_w1_w4_GE")
base_predictors <- "race + pseudo.gpa + sespc_al + nhood1_d"

# Start the sink to capture output
sink("outputs/tables/Appendix 6.txt")

# Loop through each predictor
for (predictor in predictors) {
  
  cat(paste("\nResults for predictor:", predictor, "\n"))
  
  model_list <- list()
  robust_list <- list()
  
  # Loop through each model and run regression for current predictor
  for (var in names(models)) {
    
    # Construct the formula
    unique_predictor <- models[[var]]
    if (is.null(unique_predictor)) {
      formula_str <- paste(var, "~", base_predictors, "+", predictor)
    } else {
      formula_str <- paste(var, "~", base_predictors, "+", predictor, "+", unique_predictor)
    }
    formula_obj <- as.formula(formula_str)
    
    # Run regression
    mod <- glm(formula = formula_obj, data = subset(analytical_dataset, in_sample == 1),
               weights = weights, family = "quasibinomial")
    
    # Store the model and robust SEs
    model_list[[var]] <- mod
    robust_list[[var]] <- coeftest(x = mod, vcov = vcovCL(mod, type = "HC0", cluster = ~ cluster))
  }
  
  stargazer(robust_list, type = "text", report = "vc*p", header = FALSE, 
            title = paste("Appendix 6: Results for predictor:", predictor), 
            column.labels = names(models))  
}

sink()


##########################################################################
# APPENDIX 7 REGRESSION RESULTS
#
# Appendix 7. Incidence rate ratios from negative binomial regression 
# models assessing impact of gender expression (GE) over time on adult 
# substance use over prior 30-days
##########################################################################

# Define the models
models <- list(
  w4.cigarettes.bin.30 = "w1.cigarettes",
  w4.marijuana.bin.30 = "w1.marijuana",
  w4.drunk.bin.30 = "w1.drunk",
  w4.fav.bin.30 = "w1.recreational",
  w4.prescription = NULL
)

# Define the list of predictors
predictors <- c("w1.GE_male_std", "w4.GE_male_std", "delta_w1_w4_GE")
base_predictors <- "race + pseudo.gpa + sespc_al + nhood1_d"

# Start the sink to capture output
sink("outputs/tables/Appendix 7.txt")

# Loop through each predictor
for (predictor in predictors) {
  
  cat(paste("\nResults for predictor:", predictor, "\n"))
  
  model_list <- list()
  robust_list <- list()
  
  # Loop through each model and run regression for current predictor
  for (var in names(models)) {
    
    # Construct the formula
    unique_predictor <- models[[var]]
    if (is.null(unique_predictor)) {
      formula_str <- paste(var, "~", base_predictors, "+", predictor)
    } else {
      formula_str <- paste(var, "~", base_predictors, "+", predictor, "+", unique_predictor)
    }
    formula_obj <- as.formula(formula_str)
    
    # Run regression
    mod <- glm.nb(formula = formula_obj, data = subset(analytical_dataset, in_sample == 1),
                  weights = weights)
    
    # Store the model and robust SEs
    model_list[[var]] <- mod
    robust <- coeftest(x = mod, vcov = vcovCL(mod, type = "HC0", cluster = ~ cluster))
    
    # Exponentiate to get the IRR
    robust[, "Estimate"] <- exp(robust[, "Estimate"])
    robust[, "Std. Error"] <- robust[, "Estimate"] * robust[, "Std. Error"]  # Delta method approximation
    
    robust_list[[var]] <- robust
  }
  
  stargazer(robust_list, type = "text", report = "vc*p", header = FALSE, 
            title = paste("Appendix 7: Results for predictor:", predictor), 
            column.labels = names(models))  
}

sink()



##########################################################################
# APPENDIX 8 REGRESSION RESULTS
#
# Appendix 8. Incidence rate ratios from negative binomial regression 
# models assessing impact of gender expression (GE) over time on adult 
# substance use over prior 12-months
##########################################################################

# Define the models
models <- list(
  w4.marijuana.bin.year = "w1.marijuana",
  w4.drunk.bin.year = "w1.drunk",
  w4.fav.bin.year = "w1.recreational"
)

# Define the list of predictors
predictors <- c("w1.GE_male_std", "w4.GE_male_std", "delta_w1_w4_GE")
base_predictors <- "race + pseudo.gpa + sespc_al + nhood1_d"

# Start the sink to capture output
sink("outputs/tables/Appendix 8.txt")

# Loop through each predictor
for (predictor in predictors) {
  
  cat(paste("\nResults for predictor:", predictor, "\n"))
  
  model_list <- list()
  robust_list <- list()
  
  # Loop through each model and run regression for current predictor
  for (var in names(models)) {
    
    # Construct the formula
    unique_predictor <- models[[var]]
    if (is.null(unique_predictor)) {
      formula_str <- paste(var, "~", base_predictors, "+", predictor)
    } else {
      formula_str <- paste(var, "~", base_predictors, "+", predictor, "+", unique_predictor)
    }
    formula_obj <- as.formula(formula_str)
    
    # Run regression
    mod <- glm.nb(formula = formula_obj, data = subset(analytical_dataset, in_sample == 1),
                  weights = weights)
    
    # Store the model and robust SEs
    model_list[[var]] <- mod
    robust <- coeftest(x = mod, vcov = vcovCL(mod, type = "HC0", cluster = ~ cluster))
    
    # Exponentiate to get the IRR
    robust[, "Estimate"] <- exp(robust[, "Estimate"])
    robust[, "Std. Error"] <- robust[, "Estimate"] * robust[, "Std. Error"]  # Delta method approximation
    
    robust_list[[var]] <- robust
  }
  
  stargazer(robust_list, type = "text", report = "vc*p", header = FALSE, 
            title = paste("Appendix 8: Results for predictor:", predictor), 
            column.labels = names(models))  
}

sink()
##########################################################################
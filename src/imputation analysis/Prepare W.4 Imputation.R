#-------------------------------------------------------------------------
# AUTHOR:             Jacob Jameson
# PURPOSE:            Construct and Prepare Wave 4 (Median Imputation Analysis)
#-------------------------------------------------------------------------
#
# load packages ----------------------------------------------------------
library(tidyverse)
library(haven)
library(scales)

# Data paths ------------------------------------------------------------
data_path <- '~/Add Health Data'
inhome_path <-  paste0(data_path, '/Wave IV In Home Interview Data/wave4')
weights_path <-  paste0(data_path, '/Wave IV Grand Sample Weights/weights4')

# Load the wave 4 data and wave 4 weights -------------------------------
allwave.4 <- read_xpt(paste0(inhome_path, '/wave4.xpt'))
homeweights.4 <- read_xpt(paste0(weights_path, '/weights4.xpt'))

# Merge wave 4 data with the weights ------------------------------------
wave.4 <- merge(allwave.4, homeweights.4, by='AID')

# Rename variables lowercase --------------------------------------------
names(wave.4) <- tolower(names(wave.4))

# Remove data no longer using -------------------------------------------
rm(allwave.4, homeweights.4)

#-------------------------------------------------------------------------
# Create variables that will be used in the analysis
#
#   - w4.GE_male: GE score created based on Flemming et al.
#   - w4.GE_male_std: GE score by males standardized
#-------------------------------------------------------------------------


# Variables indicated by Flemming et.all for GE --------------------------

wave.4 <- wave.4 %>%
  mutate(h4to25 = ifelse(h4to25 > 1 | h4to25 < 0, NA, h4to25),
         h4cj1 = ifelse(h4cj1 > 1 | h4cj1 < 0, NA, h4cj1),
         h4da16 = ifelse(h4da16 > 3 | h4da16 < 1, NA, h4da16),
         h4pe5 = ifelse(h4pe5 > 5 | h4pe5 < 1, NA, h4pe5),
         h4pe9 = ifelse(h4pe9 > 5 | h4pe9 < 1, NA, h4pe9),
         h4pe2 = ifelse(h4pe2 > 5 | h4pe2 < 1, NA, h4pe2),
         h4da6 = ifelse(h4da6 > 7 | h4da6 < 0, NA, h4da6),
         h4da23 = ifelse(h4da23 > 105 | h4da23 < 0, NA, h4da23),
         h4da8 = ifelse(h4da8 > 7 | h4da8 < 0, NA, h4da8),
         h4pe4 = ifelse(h4pe4 > 5 | h4pe4 < 1, NA, h4pe4),
         h4re10 = ifelse(h4re10 > 7 | h4re10 < 0, NA, h4re10),
         h4da17 = ifelse(h4da17 > 99 | h4da17 < 0, NA, h4da17),
         h4mi1 = ifelse(h4mi1 > 1 | h4mi1 < 0, NA, h4mi1),
         h4da4 = ifelse(h4da4 > 7 | h4da4 < 0, NA, h4da4),
         h4mh23 = ifelse(h4mh23 > 3 | h4mh23 < 0, NA, h4mh23),
         h4pe6 = ifelse(h4pe6 > 5 | h4pe6 < 1, NA, h4pe6),
         h4mh7 = ifelse(h4mh7 > 6 | h4mh7 < 1, NA, h4mh7),
         h4pe10 = ifelse(h4pe10 > 5 | h4pe10 < 1, NA, h4pe10),
         h4pe35 = ifelse(h4pe35 > 5 | h4pe35 < 1, NA, h4pe35),
         h4da11 = ifelse(h4da11 > 1 | h4da11 < 0, NA, h4da11),
         h4pe22 = ifelse(h4pe22 > 5 | h4pe22 < 1, NA, h4pe22),
         h4pe26 = ifelse(h4pe26 > 5 | h4pe26 < 1, NA, h4pe26))


# Variables indicated by Flemming et.all for GE

ge.vars <- c('h4to25', 'h4cj1', 'h4da16', 'h4pe5', 'h4pe9', 'h4pe2', 'h4da6',
             'h4da23', 'h4da8', 'h4pe4', 'h4re10', 'h4da17', 'h4mi1', 'h4da4',
             'h4mh23', 'h4pe6', 'h4mh7', 'h4pe10',  'h4pe35', 'h4da11', 
             'h4pe22', 'h4pe26')


####################################################################

# create a data frame with the ge.vars
ge.df <- wave.4[, c('bio_sex4', ge.vars)]
colSums(is.na(filter(ge.df, bio_sex4 == 1)))

# count the number of missing values by row (ID)
missing <- rowSums(is.na(ge.df)) > 0
num_missing <- sum(missing)

wave.4 <- wave.4 %>% 
  group_by(bio_sex4) %>% 
  mutate_at(ge.vars, ~ifelse(is.na(.), median(., na.rm = TRUE), .))

# Logistic regression to predict male using GE variables ------------
wave.4$w4_male <- factor(ifelse(wave.4$bio_sex4 == 1, 1, 0))
wave.4$w4_female <- factor(ifelse(wave.4$bio_sex4 == 2, 1, 0))

predict_male <- glm(w4_male ~ factor(h4to25) + factor(h4cj1) + 
                      factor(h4da16) + factor(h4pe5) + 
                      factor(h4pe9) + factor(h4pe2) + factor(h4da6) + 
                      h4da23 + factor(h4da8) + factor(h4pe4) + 
                      factor(h4re10) + h4da17 + factor(h4mi1) + 
                      factor(h4da4) + factor(h4mh23) + factor(h4pe6) + 
                      factor(h4mh7) + factor(h4pe10) + factor(h4pe35) + 
                      factor(h4da11) + factor(h4pe22) + factor(h4pe26),
                    data = wave.4, family = "binomial")

# Get a prediction score for male. This will be our GE score ------------
wave.4$w4.GE_male <- predict(predict_male, wave.4, type="response")

wave.4 <- wave.4 %>%
  group_by(w4_male) %>%
  mutate(w4.GE_male_std = scale(w4.GE_male)  %>% as.vector) %>%
  ungroup()

#-------------------------------------------------------------------------

wave.4 <- wave.4 %>%
  mutate(edu = case_when(
    h4ed2 %in% c("1", "2") ~ "Some HS or Less",
    h4ed2 == "3" ~ "HS Diploma/GED",
    h4ed2 %in% c("4", "5", "6") ~ "Some College or Tech/Assoc Degree",
    h4ed2 %in% c("7", "8", "9", "10", "11", "12", "13") ~ "College Degree or More",
    TRUE ~ NA_character_
  ))

# convert to factor with specified labels
wave.4$edu <- factor(wave.4$edu, levels = c("Some HS or Less", "HS Diploma/GED",
                                            "Some College or Tech/Assoc Degree",
                                            "College Degree or More"))


wave.4 <- wave.4 %>% 
  mutate(insurance = case_when(
    h4hs1 %in% c(1, 11, 96, 98) ~ "Uninsured or DK",
    h4hs1 %in% c(2, 3, 4, 5, 6, 7, 8) ~ "Private or Employer-Based",
    h4hs1 %in% c(9) ~ "Medicaid/Medicare",
    h4hs1 %in% c(10, 7) ~ "Other Government"
  ))

wave.4$insurance <- factor(wave.4$insurance,
                           levels = c("Private or Employer-Based", 
                                      "Medicaid/Medicare", 
                                      "Other Government",
                                      "Uninsured or DK"))

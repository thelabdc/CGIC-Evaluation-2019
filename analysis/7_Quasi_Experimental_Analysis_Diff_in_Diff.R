###############################################################################
# Author: FYO
# CGIC evaluations: Implementing DiD design
# 
# Last Updated: 
###############################################################################

# Packages needed for analysis
library(tidyverse)
library(vtable)
library(lubridate)
library(digest)
library(purrr)
library(readr)
library(reshape)
library(lmtest)
library(stargazer)
library(R.utils)

################################ Loading data #################################

# To load, open this file AFTER loading the .Rproj file, which sets the 
# working directory to /CGIC_Evaluation_2019/analysis/

shotspotter  <- read.delim(file = "../data/did_clean_shotspotter.csv", sep = ",", stringsAsFactors = FALSE)
calls  <- read.delim(file = "../data/did_calls_df.csv", sep = ",", stringsAsFactors = FALSE)
dcr_violent <- read.delim(file = "../data/did_violent_crime_df.csv", sep = ",", stringsAsFactors = FALSE)
dcr_gun <- read.delim(file = "../data/did_gun_crime_df.csv", sep = ",", stringsAsFactors = FALSE)
gun_broad_arrests <- read.delim(file = "../data/did_broad_arrests_df.csv", sep = ",", stringsAsFactors = FALSE)
gun_weapons_violation_arrests <- read.delim(file = "../data/did_violent_gun_arrests_df.csv", sep = ",", stringsAsFactors = FALSE)
gun_violent_arrests <- read.delim(file = "../data/did_violent_gun_arrests_df.csv", sep = ",", stringsAsFactors = FALSE)


# making variable names lower case
names(shotspotter) <- tolower(names(shotspotter))
names(calls) <- tolower(names(calls))
names(dcr_violent) <- tolower(names(dcr_violent))
names(dcr_gun) <- tolower(names(dcr_gun))
names(gun_broad_arrests) <- tolower(names(gun_broad_arrests))
names(gun_weapons_violation_arrests) <- tolower(names(gun_weapons_violation_arrests))
names(gun_violent_arrests) <- tolower(names(gun_violent_arrests))


######## Additional processing steps for regressions and rolling means ########

data_prep <- function(outcome_raw, outcome_count_name) {
  
  # Description: This fxn takes raw data on outcomes and generates daily counts by PSA
  # Args: outcome_raw = file containing raw data for given outcome
  #       outcome_count_name = what to name the processed outcome variable
  # Returns: df containing the daily counts 
  
  daily_outcome <- outcome_raw %>% 
    filter(event_psa > 0 & !is.na(event_psa)) %>% # drops obs. w/ PSA == -1 or missing PSA
    mutate(event_date_time = ymd_hms(event_time)) %>%
    mutate(event_date = date(event_date_time)) %>%
    group_by(event_date, event_psa) %>%
    summarise(outcome_count = n()) %>%
    mutate(Treated = substr(event_psa, 1,1) == 7) %>%
    mutate(year = year(event_date)) %>%
    mutate(month = month(event_date)) %>%
    mutate(month_year =  ymd(paste(year, month, "01", sep = "-")))
  
  # workaround for naming count variable based on passed in fxn argument
  daily_outcome[,outcome_count_name] <- daily_outcome$outcome_count
  daily_outcome <- daily_outcome %>% select(-outcome_count)
  
  return(daily_outcome)
}

daily_shotspotter <- data_prep(shotspotter, "shotspotter_count")
daily_calls <- data_prep(calls, "calls_count")
daily_dcr_violent <- data_prep(dcr_violent, "dcr_violent_count")
daily_dcr_gun <- data_prep(dcr_gun, "dcr_gun_count")
daily_arrest_broad <- data_prep(gun_broad_arrests, "arrest_broad_count")
daily_arrest_weapon_violation <- data_prep(gun_weapons_violation_arrests, "arrest_weapon_violation_count")
daily_arrest_violent <- data_prep(gun_violent_arrests, "arrest_violent_count")


# Saving results as .csv for later use (e.g. calculating rolling means)
write.csv(daily_shotspotter, file = "../data/daily_shotspotter.csv", row.names = FALSE)
write.csv(daily_calls, file = "../data/daily_calls.csv", row.names = FALSE)
write.csv(daily_dcr_violent, file = "../data/daily_dcr_violent.csv", row.names = FALSE)
write.csv(daily_dcr_gun, file = "../data/daily_dcr_gun.csv", row.names = FALSE)
write.csv(daily_arrest_broad, file = "../data/daily_arrest_broad.csv", row.names = FALSE)
write.csv(daily_arrest_weapon_violation, file = "../data/daily_arrest_violation.csv", row.names = FALSE)
write.csv(daily_arrest_violent, file = "../data/daily_arrest_violent.csv", row.names = FALSE)


######################### Processing for regressions ##########################

### Prepping data frame for adding zeros on days with no events
start_date = ymd('2016-11-01')
end_date = ymd('2019-04-30')

# Enumerating list of PSAs
psa_list <- daily_calls %>% group_by(event_psa) %>% tally() %>% filter(event_psa > 0) %>% select(event_psa)

# Enumerating list of dates to be included
date_list <- seq(start_date, end_date, by = "days")

# Creating df with date/psa combinations to merge with other data
df <- data.frame(event_date = rep(date_list, nrow(psa_list)),
                 event_psa = rep(psa_list$event_psa, each = length(date_list)))


# Function that merges with exhaustive list of PSA-date combos, replaces NAs with 0
add_zeros <- function(outcome_df) {
  
  # Description: This function adds zeros to the outcome data
  # Args: outcome_df = file containing processed daily data for given outcome
  # Returns: df containing the daily counts with zeros for use in regressions
  
  # # Joining with full date list
   outcome_df <- outcome_df %>%
     right_join(df, by = c("event_date", "event_psa")) %>%
     select(-year, -month, -month_year)
  
  # Replacing with zeros
  outcome_df <- outcome_df %>%
    mutate(days_since_treatment = event_date - ymd("2017-11-01")) %>%
    mutate(quarters_since_treatment = floor(days_since_treatment / 90)) %>%
    mutate(Treated = (substr(event_psa, 1,1) == 7))
  
  # Replacing NAs with with zeros
  outcome_df[is.na(outcome_df)] <- 0
  
  return(outcome_df)
}

daily_shotspotter <- add_zeros(daily_shotspotter) %>%
  filter(event_psa < 200 | event_psa >= 300) # Removing 2D psas since shotspotter is not deployed in that district
daily_calls <- add_zeros(daily_calls)
daily_dcr_violent <- add_zeros(daily_dcr_violent)
daily_dcr_gun <- add_zeros(daily_dcr_gun)
daily_arrest_broad <- add_zeros(daily_arrest_broad)
daily_arrest_weapon_violation <- add_zeros(daily_arrest_weapon_violation)
daily_arrest_violent <- add_zeros(daily_arrest_violent)

# Saving results as .Rdata objects for use in regression function
R.utils::saveObject(daily_shotspotter, file = "../data/daily_shotspotter_reg.RData")
R.utils::saveObject(daily_calls, file = "../data/daily_calls_reg.RData")
R.utils::saveObject(daily_dcr_violent, file = "../data/daily_dcr_violent_reg.RData")
R.utils::saveObject(daily_dcr_gun, file = "../data/daily_dcr_gun_reg.RData")
R.utils::saveObject(daily_arrest_broad, file = "../data/daily_arrest_broad_reg.RData")
R.utils::saveObject(daily_arrest_weapon_violation, file = "../data/daily_arrest_violation_reg.RData")
R.utils::saveObject(daily_arrest_violent, file = "../data/daily_arrest_violent_reg.RData")

########################### Estimating Regressions ############################

# Function for generating confidence intervals. This is called by the estimate_models function
generate_ci <- function(results, outname) {
  # Description: This takes regression estimates and generates 95% confidence intervals
  # Args: results = object containing regression results
  #       outname = file name for output table
  # Returns: object containing confidence interval results
  
  
  # Estimate confidence intervals
  ci <- as.data.frame(confint(results))
  
  ci$variable_names <- row.names(ci)
  ci$estimate <- (ci$`2.5 %` + ci$`97.5 %`)/2
  
  # Export results as .csv
  output_name <- paste('../output/', outname, ".csv", sep = "")
  write.csv(ci, file = output_name, row.names = TRUE)
  
  # Return value
  return(ci)
}


# Function for running regression analysis
estimate_models <- function(data_path = "../data/",
                            data_file = "daily_shotspotter_reg.RData",
                            outcome = "shotspotter_count", 
                            ci_output_name = "shotspotter_ci",
                            model = " ~ factor(event_psa) + factor(quarters_since_treatment) + factor(Treated) + factor(Treated):factor(quarters_since_treatment)") {
  
  
  # Description: Estimates regression model, generates confidence interval for estimates
  #Args: data_path = file path for input data
  #      data_file = name of file containing input data.
  #                   Must be .RData object saved using "saveObject" command from R.utils package
  #       outcome = name of dependant variable
  #       model = model string to be called by lm().
  #Returns: object containing coefficient estimates and clustered standard errors
  
  
  # Importing data for regression as temp file (must use loadObject from R.utils to rename while loading)
  temp <<- loadObject(paste(data_path, data_file, sep = ""))
  
  # Adding outcome of interest to model string
  model_string <- paste(outcome, model, sep = "")
  
  # Estimaing the model
  model_results.raw <- lm(model_string, data = temp)
  
  # Estimating the clustered variance/covariance matrix
  model_results.vcov <- sandwich::vcovCL(model_results.raw,cluster = ~ event_psa)
  #model_results.vcov_test <- cluster.vcov(model_results.raw, data_df$event_psa) These are equivalent, the first is faster
  
  # Returning estimates with robust standard errors
  model_results.clustered <- coeftest(model_results.raw, vcov = model_results.vcov)
  
  # Generating confidence intervals as well.
  model_results.ci <- generate_ci(model_results.clustered, ci_output_name)
  
  # Returning regression results for stargazer
  return(model_results.clustered)
 
}


# Write confidence interval results to csv in the "outputs" folder
shotspotter_model.results <- estimate_models(data_file = "daily_shotspotter_reg.Rdata", 
                                             outcome = "shotspotter_count", 
                                             ci_output_name = "shotspotter_ci")
calls_model.results <- estimate_models(data_file = "daily_calls_reg.RData",
                                       outcome = "calls_count", 
                                       ci_output_name = "calls_ci")
dcr_violent_model.results <- estimate_models(data_file = "daily_dcr_violent_reg.RData", 
                                             outcome = "dcr_violent_count", 
                                             ci_output_name = "dcr_violent_ci")
dcr_gun_model.results <- estimate_models(data_file = "daily_dcr_gun_reg.RData", 
                                         outcome = "dcr_gun_count", 
                                         ci_output_name = "dcr_gun_ci")
arrest_broad_model.results <- estimate_models(data_file = "daily_arrest_broad_reg.RData", 
                                              outcome = "arrest_broad_count", 
                                              ci_output_name = "arrest_broad_ci")
arrest_violation_model.results <- estimate_models(data_file = "daily_arrest_violation_reg.RData", 
                                                  outcome = "arrest_weapon_violation_count", 
                                                  ci_output_name = "arrest_violation_ci")
arrest_violent_model.results <- estimate_models(data_file = "daily_arrest_violent_reg.RData", 
                                                outcome = "arrest_violent_count", 
                                                ci_output_name = "arrest_violent_ci")
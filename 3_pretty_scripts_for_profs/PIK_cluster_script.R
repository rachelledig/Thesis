## Rachel Ledig
#Thesis - Bayesian Hierarchial Model of Irrigaton Expansion

## Script for PIK Cluster

####################################################################################
####################################################################################
####################################################################################

#Notes
#Fix file paths
#Force install packages?


####################################################################################
####################################################################################
####################################################################################




## Packages
library(dplyr)
library(tidyr)
library(purrr)
library(brms)

## Plotting time series per country

#load unstandardized data
aei <- read.csv("/Volumes/RachelExternal/Thesis/Data_and_Plots/aei_raw.csv")

#nest data
by_country <-
  aei %>%
  group_by(Countryname) %>%
  nest()

#run first model based on time series
id_fits_time_zib <-
  brm(data = by_country$data[[1]],
      formula = irrcrop ~ 1 + yearcount,
      family = zero_one_inflated_beta(),
      prior = 
      control = list(adapt_delta = 0.95, 
                     max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/Country_fits_pred_intercept_years.Rds")

# map to all other countries
yearcount_fits <- 
  by_country %>%
  mutate(model = map(data, ~update(id_fits_time_zib, newdata = ., seed = 2)))

#save nested df for later analysis on local machine
saveRDS(yearcount_fits, "/Volumes/RachelExternal/Thesis/Data_and_Plots/Trial2.Rds")

####################################################################################
####################################################################################
####################################################################################

# Simple Model 1- yearcount, country, ID

####################################################################################

#Unconditional Means Model ID

uncond_means_id <-
  brm(data = aei_std,
      formula = irrcrop ~ 1 + (1|id),
      family = zero_one_inflated_beta(),
      prior = 
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_uncond_means_id.Rds")

####################################################################################

#Unconditional Growth Model ID

uncond_growth_id <-
  brm(data = aei_std,
      formula = irrcrop ~ 1 + yearcount + (1 + yearcount|id),
      family = zero_one_inflated_beta(),
      prior = 
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_uncond_growth_id.Rds")

####################################################################################

#Full Simple model ID


full_simple_id <-
  brm(data = aei_std,
      formula = irrcrop ~ 0 + Intercept + yearcount + ISO + ISO:yearcount (1 + yearcount|id),
      family = zero_one_inflated_beta(),
      prior = 
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_full_id.Rds")


####################################################################################
####################################################################################
####################################################################################

# Simple Model - yearcount, regions, Country

####################################################################################

#Unconditional Means Model ISO

uncond_means_ISO <-
  brm(data = aei_std,
      formula = irrcrop ~ 1 + (1|ISO),
      family = zero_one_inflated_beta(),
      prior = 
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_uncond_means_ISO.Rds")

####################################################################################

#Unconditional Growth Model ISO

uncond_growth_ISO <-
  brm(data = aei_std,
      formula = irrcrop ~ 1 + yearcount + (1 + yearcount|ISO),
      family = zero_one_inflated_beta(),
      prior = 
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_uncond_growth_ISO.Rds")

####################################################################################

#Full Simple model ISO


full_simple_ISO <-
  brm(data = aei_std,
      formula = irrcrop ~ 0 + Intercept +  + six_regions + six_regions:yearcount (1 + yearcount|ISO),
      family = zero_one_inflated_beta(),
      prior = 
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_full_ISO.Rds")


####################################################################################
####################################################################################
####################################################################################

# Complex model 

full_simple <-
  brm(data = aei_std,
      formula = irrcrop ~ 0 + Intercept +  + six_regions + six_regions:yearcount (1 + yearcount|ISO),
      family = zero_one_inflated_beta(),
      prior = 
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_full.Rds")

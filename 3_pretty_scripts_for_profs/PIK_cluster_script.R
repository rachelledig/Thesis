## Rachel Ledig
#Thesis - Bayesian Hierarchical Model of Irrigaton Expansion

## Script for PIK Cluster

####################################################################################
####################################################################################
####################################################################################

#Notes
#Fix file paths
#Force install packages?
#saving outputs?
#cores? 


####################################################################################
####################################################################################
####################################################################################




## Packages
library(dplyr)
library(brms)

## Individual cell change model

#load unstandardized data
aei_std <- read.csv("/Volumes/RachelExternal/Thesis/Data_and_Plots/aei_std.csv")

#nest data
by_country <-
  aei_std %>%
  group_by(Countryname) %>%
  nest()

#run first model based on time series
id_fits_time_zib <-
  brm(data = by_country$data[[1]],
      formula = bf(irrcrop ~ 1 + yearcount,
                   coi ~ 1 + yearcount,
                   coi ~ 1 + yearcount),
      family = zero_one_inflated_beta(),
      prior = c(prior(normal(0, 100), class = b)), #use default priors except for slopes
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


####################################################################################
####################################################################################
####################################################################################

# Simple Model 1- yearcount, country, ID

####################################################################################

#Unconditional Means Model ID
#aims to look at general effects across cells without time as a component

#default priors are fine here

uncond_means_id <-
  brm(data = aei_std,
      formula = bf(irrcrop ~ 1 + (1|id),
                   coi ~ 1 + (1|id),
                   zoi ~ 1 + (1|id)),
      family = zero_one_inflated_beta(),
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_uncond_means_id.Rds")

####################################################################################

#Unconditional Growth Model ID
#aims to look at general effects across cells and over time prior to adding other preds

uncond_growth_id <-
  brm(data = aei_std,
      formula = bf(irrcrop ~ 1 + yearcount + (1 + yearcount|id),
                   coi ~ 1 + yearcount + (1 + yearcount|id),
                   zoi ~ 1 + yearcount + (1 + yearcount|id)),
      family = zero_one_inflated_beta(),
      prior = c(prior(normal(0, 100), class = b)),
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
      formula = bf(irrcrop ~ 0 + Intercept + yearcount + ISO + ISO:yearcount (1 + yearcount|id),
                   coi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount (1 + yearcount|id),
                   zoi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount (1 + yearcount|id)),
      family = zero_one_inflated_beta(),
      prior = c(prior(normal(0, 100), class = b)),
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
      formula = bf(irrcrop ~ 1 + (1|ISO),
                   coi ~ 1 + (1|ISO),
                   zoi ~ 1 + (1|ISO)),
      family = zero_one_inflated_beta(),
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
      formula = bf(irrcrop ~ 1 + yearcount + (1 + yearcount|ISO),
                   zoi ~ 1 + yearcount + (1 + yearcount|ISO),
                   coi ~ 1 + yearcount + (1 + yearcount|ISO)),
      family = zero_one_inflated_beta(),
      prior = c(prior(normal(0, 100), class = b)),
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
      formula = bf(irrcrop ~ 0 + Intercept +  + six_regions + yearcount + six_regions:yearcount (1 + yearcount|ISO)
                   coi ~ 0 + Intercept +  + six_regions + yearcount + six_regions:yearcount (1 + yearcount|ISO),
                   zoi ~ 0 + Intercept +  + six_regions + yearcount + six_regions:yearcount (1 + yearcount|ISO)),
      family = zero_one_inflated_beta(),
      prior = c(prior(normal(0, 100), class = b)),
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

#create a list of data frames to run with brm_multiple

num_groups <- 10
folds <- 
  aei_std %>% 
  drop_na(years, rugged, precip, gdppc, popdens, dist, six_regions, pet, DD.regime) %>% 
  group_by((row_number()-1) %/% (n()/num_groups)) %>%
  nest %>% 
  pull(data)


####################################################################################
####################################################################################
####################################################################################

# Complex model 
# without interactions, at lower levels they are easier to sort out but here its too complex

full_simple <-
  brm_multiple(data = folds,
               formula = bf(irrcrop ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id),
                            coi ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id),
                            zoi ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id)),
               family = zero_one_inflated_beta(),
               prior = c(prior(normal(0, 100), class = b)),
               control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
               cores = 2, 
               iter = 1000, 
               chains = 4,
               seed = 2,
               file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/complex_full.Rds")

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
library(tidyr)
library(brms)


# changes <- aei %>% 
#   subset(yearcount == 0 | yearcount == 45) %>% 
#   group_by(id) %>% 
#   mutate(change = irrfrac - lag(irrfrac, default = irrfrac[1])) %>% 
#   subset(change > 0)
# 
# changes <- changes[changes$irrfrac==changes$change, ] #all of these started from 0
# 
# ids <- changes[sample(1:nrow(changes), 50), ]
# ids <- as.list(ids$id)
# 
# aei_small <- aei %>% 
#   filter(., id %in% ids) %>% 
#   subset(., subset = yearcount %in% c("0",  "5", "10", "15", "20", "25", "30", "35", "40", "45"))
# 
# # #high performing cells
# ids <- c("35802", "47090", "46708", "46832", "35802", "46710", "46707", "46961", "46833", "46962")

# aei_small_std <-
#   aei_small %>% 
#   mutate(across(c(17:19, 21), centered)) %>% 
#   mutate(across(c(23), normalized)) %>% 
#   mutate(across(c(2,3,6,26), as.factor)) %>% 
#   mutate(yearcount = yearcount +1)


## Individual cell change model

#load unstandardized data
#aei_std <- read.csv("/home/ledig/thesis/aei_std.csv")

# #nest data
# by_country <-
#   aei_std %>%
#   group_by(Countryname) %>%
#   nest()

#run first model based on time series
id_fits_time_zib <-
  brm(data = by_country$data[[1]],
      formula = bf(irrfrac ~ 1 + yearcount,
                   zi ~ 1 + yearcount),
      family = zero_inflated_beta(),
      prior = c(prior(normal(0, 100), class = b)), #use default priors except for slopes
      control = list(adapt_delta = 0.95, 
                     max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/Country_fits_pred_intercept_years_small.Rds")

# map to all other countries
yearcount_fits <- 
  by_country %>%
  mutate(model = map(data, ~update(id_fits_time_zib, newdata = ., seed = 2)))

#save nested df for later analysis on local machine
saveRDS(yearcount_fits, "/Volumes/RachelExternal/Thesis/Data_and_Plots/Trial2_small.Rds")

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
  brm(data = aei_small_std,
      formula = bf(irrfrac ~ 1 + (1|id),
                   zi ~ 1 + (1|id)),
      family = zero_inflated_beta(),
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_uncond_means_id_small.Rds")

####################################################################################

#Unconditional Growth Model ID
#aims to look at general effects across cells and over time prior to adding other preds

uncond_growth_id <-
  brm(data = aei_small_std,
      formula = bf(irrfrac ~ 1 + yearcount + (1 + yearcount|id),
                   zi ~ 1 + yearcount + (1 + yearcount|id)),
      family = zero_inflated_beta(),
      prior = c(prior(normal(0, 100), class = b)),
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      inits = "0", #had trouble initializing beta_lpdf. second shape parameter was 0. Solution found here https://discourse.mc-stan.org/t/rejecting-initial-value/7152/4
      iter = 1000, 
      chains = 4,
      seed = 348,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_uncond_growth_id_small.Rds")

####################################################################################

#Full Simple model ID


full_simple_id <-
  brm(data = aei_small_std,
      formula = bf(irrfrac ~ 0 + Intercept + yearcount + ISO + ISO:yearcount + (1 + yearcount|id), 
                   zi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount + (1 + yearcount|id)),
      family = zero_inflated_beta(),
      prior = c(prior(normal(0, 100), class = b)),
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      inits = "0",
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_full_id_small2.Rds")


####################################################################################
####################################################################################
####################################################################################

# Simple Model - yearcount, regions, Country

####################################################################################

#Unconditional Means Model ISO
job::job({
uncond_means_ISO <-
  brm(data = aei_small_std,
      formula = bf(irrfrac ~ 1 + (1|ISO),
                   zi ~ 1 + (1|ISO)),
      family = zero_inflated_beta(),
      control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_uncond_means_ISO_small.Rds")
})


####################################################################################

#Unconditional Growth Model ISO
job::job({
uncond_growth_ISO <-
  brm(data = aei_small_std,
      formula = bf(irrfrac ~ 1 + yearcount + (1 + yearcount|ISO),
                   zi ~ 1 + yearcount + (1 + yearcount|ISO)),
      family = zero_inflated_beta(),
      prior = c(prior(normal(0, 100), class = b)),
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      inits = "0",
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_uncond_growth_ISO_small.Rds")
})
####################################################################################

#Full Simple model ISO


full_simple_ISO <-
  brm(data = aei_small_std,
      formula = bf(irrfrac ~ 0 + Intercept +  + six_regions + yearcount + six_regions:yearcount + (1 + yearcount|ISO)
                   zi ~ 0 + Intercept +  + six_regions + yearcount + six_regions:yearcount + (1 + yearcount|ISO)),
      family = zero_inflated_beta(),
      prior = c(prior(normal(0, 100), class = b)),
        control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
      cores = 2, 
      iter = 1000, 
      chains = 4,
      seed = 2,
      file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/simple_full_ISO_small.Rds")


####################################################################################
####################################################################################
####################################################################################

#create a list of data frames to run with brm_multiple

num_groups <- 10
folds <- 
  aei_small_std %>% 
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
               formula = bf(irrfrac ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id),
                            zi ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id)),
               family = zero_inflated_beta(),
               prior = c(prior(normal(0, 100), class = b)),
               control = list(adapt_delta = 0.95, 
                       max_treedepth=11),
               cores = 2, 
               iter = 1000, 
               chains = 4,
               seed = 2,
               file = "/Volumes/RachelExternal/Thesis/Thesis/3_pretty_scripts_for_profs/Fits/complex_full_small.Rds")

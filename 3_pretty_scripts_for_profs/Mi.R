## Packages
library(dplyr)
library(rstan)
library(brms)


years <- seq(5, 45, by = 10)

aei_raw_test <- aei[aei$yearcount %in% years, ]


normalized<-function(y) {
  x<-y[!is.na(y)]
  x<-(x - min(x)) / (max(x) - min(x))
  y[!is.na(y)]<-x
  return(y)
}

centered <-function(y) {
  x<-y[!is.na(y)]
  x<-(x - mean(x)) / (sd(x)*2)
  y[!is.na(y)]<-x
  return(y)
}

aei_std_test <-
  aei_raw_test %>% 
  mutate(across(c(17:19, 21), centered)) %>% 
  mutate(across(c(23), normalized)) %>% 
  mutate(across(c(2,3,6,26), as.factor))


summary(aei_std_test)


#Unconditional Means Model ID####################################################################################
#aims to look at general effects across cells without time as a component

uncond_means_data_id <- 
  make_standata(data = aei_std_test,
                formula = bf(irrcrop ~ 1 + (1|id),
                             coi ~ 1 + (1|id),
                             zoi ~ 1 + (1|id)),
                family = zero_one_inflated_beta())

uncond_means_code_id <-
  make_stancode(data = aei_std_test,
                formula = bf(irrcrop ~ 1 + (1|id),
                             coi ~ 1 + (1|id),
                             zoi ~ 1 + (1|id)),
                family = zero_one_inflated_beta(),
                sample_prior = "only")

uncond_means_stan_id <-
  stan(model_code = uncond_means_code_id, 
       data = uncond_means_data_id, 
       iter = 2000, 
       chains = 4,
       cores = 4,
       seed = 17,
       control = list(adapt_delta = 0.9, 
                      max_treedepth=11))

saveRDS(uncond_means_stan_id, "/p/projects/open/ledig/stan_uncond_means_id.rds")

#Unconditional Growth Model ID####################################################################################
#aims to look at general effects across cells and over time prior to adding other preds

uncond_growth_data_id <- 
  make_standata(data = aei_std,
                formula = bf(irrcrop ~ 1 + yearcount + (1 + yearcount|id),
                             coi ~ 1 + yearcount + (1 + yearcount|id),
                             zoi ~ 1 + yearcount + (1 + yearcount|id)),
                family = zero_one_inflated_beta())

uncond_growth_code_id <-
  make_stancode(data = aei_std,
                formula = bf(irrcrop ~ 1 + yearcount + (1 + yearcount|id),
                             coi ~ 1 + yearcount + (1 + yearcount|id),
                             zoi ~ 1 + yearcount + (1 + yearcount|id)),
                family = zero_one_inflated_beta(),
                prior = c(prior(normal(0, 100), class = b)),
                sample_prior = "only")


uncond_growth_stan_id <-
  stan(model_code = uncond_growth_code_id,
       data = uncond_growth_data_id,
       iter = 2000, 
       chains = 4,
       cores = 4,
       seed = 17,
       control = list(adapt_delta = 0.9, 
                      max_treedepth=11))

saveRDS(uncond_growth_stan_id, "/p/projects/open/ledig/stan_uncond_growth_id.rds")

#Full Simple model ID####################################################################################

full_simple_data_id <- 
  make_standata(data = aei_std,
                formula = formula = bf(irrcrop ~ 0 + Intercept + yearcount + ISO + ISO:yearcount (1 + yearcount|id),
                                       coi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount (1 + yearcount|id),
                                       zoi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount (1 + yearcount|id)),
                family = zero_one_inflated_beta())

full_simple_code_id <-
  make_stancode(data = aei_std,
                formula = formula = bf(irrcrop ~ 0 + Intercept + yearcount + ISO + ISO:yearcount (1 + yearcount|id),
                                       coi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount (1 + yearcount|id),
                                       zoi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount (1 + yearcount|id)),
                family = zero_one_inflated_beta(),
                prior = c(prior(normal(0, 100), class = b)),
                sample_prior = "only")


full_simple_stan_id <-
  stan(model_code = full_simple_code_id,
       data = full_simple_data_id,
       iter = 2000, 
       chains = 4,
       cores = 4,
       seed = 17,
       control = list(adapt_delta = 0.9, 
                      max_treedepth=11))

saveRDS(full_simple_stan_id, "/p/projects/open/ledig/stan_uncond_growth_id.rds")



#Unconditional Means Model ISO####################################################################################



n

saveRDS(uncond_means_stan_id, "/p/projects/open/ledig/stan_uncond_means_ISO.rds")

#Unconditional Growth Model ISO####################################################################################



uncond_growth_data_ISO <- 
  make_standata(data = aei_std,
                formula = bf(irrcrop ~ 1 + yearcount + (1 + yearcount|ISO),
                             coi ~ 1 + yearcount + (1 + yearcount|ISO),
                             zoi ~ 1 + yearcount + (1 + yearcount|ISO)),
                family = zero_one_inflated_beta())

uncond_growth_code_ISO <-
  make_stancode(data = aei_std,
                formula = bf(irrcrop ~ 1 + yearcount + (1 + yearcount|ISO),
                             coi ~ 1 + yearcount + (1 + yearcount|ISO),
                             zoi ~ 1 + yearcount + (1 + yearcount|ISO)),
                family = zero_one_inflated_beta(),
                prior = c(prior(normal(0, 100), class = b)),
                sample_prior = "only")


uncond_growth_stan_ISO <-
  stan(model_code = uncond_growth_code_ISO,
       data = uncond_growth_data_ISO,
       iter = 2000, 
       chains = 4,
       cores = 4,
       seed = 17,
       control = list(adapt_delta = 0.9, 
                      max_treedepth=11))

saveRDS(uncond_growth_stan_id, "/p/projects/open/ledig/stan_uncond_growth_ISO.rds")

#Full Simple model ISO####################################################################################


full_simple_data_ISO <- 
  make_standata(data = aei_std,
                formula = bf(irrcrop ~ 0 + Intercept + yearcount + six_regions + six_regions:yearcount (1 + yearcount|ISO),
                             coi ~ 0 + Intercept + yearcount + six_regions + six_regions:yearcount (1 + yearcount|ISO),
                             zoi ~ 0 + Intercept + yearcount + six_regions + six_regions:yearcount (1 + yearcount|ISO)),
                family = zero_one_inflated_beta())

full_simple_code_ISO <-
  make_stancode(data = aei_std,
                formula = bf(irrcrop ~ 0 + Intercept + yearcount + six_regions + six_regions:yearcount (1 + yearcount|ISO),
                             coi ~ 0 + Intercept + yearcount + six_regions + six_regions:yearcount (1 + yearcount|ISO),
                             zoi ~ 0 + Intercept + yearcount + six_regions + six_regions:yearcount (1 + yearcount|ISO)),
                family = zero_one_inflated_beta(),
                prior = c(prior(normal(0, 100), class = b)),
                sample_prior = "only")


full_simple_stan_id <-
  stan(model_code = full_simple_code_ISO,
       data = full_simple_data_ISO,
       iter = 2000, 
       chains = 4,
       cores = 4,
       seed = 17,
       control = list(adapt_delta = 0.9, 
                      max_treedepth=11))

saveRDS(full_simple_stan_id, "/p/projects/open/ledig/stan_full_simple_ISO.rds")



# Complex model####################################################################################

# without interactions, at lower levels they are easier to sort out but here its too complex

full_complex_data <- 
  make_standata(data = aei_std,
                formula = bf(irrcrop ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id),
                             coi ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id),
                             zoi ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id)),
                family = zero_one_inflated_beta())

full_complex_code <-
  make_stancode(data = aei_std,
                formula = bf(irrcrop ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id),
                             coi ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id),
                             zoi ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id)),
                family = zero_one_inflated_beta(),
                prior = c(prior(normal(0, 100), class = b)),
                sample_prior = "only")


full_complex_stan <-
  stan(model_code = full_complex_code,
       data = full_complex_data,
       iter = 2000, 
       chains = 4,
       cores = 4,
       seed = 17,
       control = list(adapt_delta = 0.9, 
                      max_treedepth=11))

saveRDS(full_simple_stan_id, "/p/projects/open/ledig/stan_full_complex.rds")




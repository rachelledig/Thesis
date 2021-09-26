## Packages
library(dplyr)
library(rstan)
library(brms)


changes <- aei %>% 
  subset(yearcount == 0 | yearcount == 45) %>% 
  group_by(id) %>% 
  mutate(change = irrfrac - lag(irrfrac, default = irrfrac[1])) %>% 
  subset(change > 0)

changes <- changes[changes$irrfrac==changes$change, ] #all of these started from 0

ids <- changes[sample(1:nrow(changes), 50), ]
ids <- as.list(ids$id)

aei_small <- aei %>% 
  filter(., id %in% ids) %>% 
  subset(., subset = yearcount %in% c("0",  "5", "10", "15", "20", "25", "30", "35", "40", "45"))

# #high performing cells
# ids <- c("35802", "47090", "46708", "46832", "35802", "46710", "46707", "46961", "46833", "46962")

aei_small_std <-
  aei_small %>% 
  mutate(across(c(17:19, 21), centered)) %>% 
  mutate(across(c(23), normalized)) %>% 
  mutate(across(c(2,3,6,26), as.factor)) %>% 
  mutate(yearcount = yearcount +1)


summary(aei_small_std)

#Unconditional Means Model ID####################################################################################
#aims to look at general effects across cells without time as a component

uncond_means_data_id <- 
  make_standata(data = aei_small_std,
                formula = bf(irrcrop ~ 1 + (1|id),
                             coi ~ 1 + (1|id),
                             zoi ~ 1 + (1|id)),
                family = zero_one_inflated_beta())

uncond_means_code_id <-
  make_stancode(data = aei_small_std,
                formula = bf(irrcrop ~ 1 + (1|id),
                             coi ~ 1 + (1|id),
                             zoi ~ 1 + (1|id)),
                family = zero_one_inflated_beta(),
                sample_prior = "only")

job::job({
  uncond_means_stan_id <-
  stan(model_code = uncond_means_code_id, 
       data = uncond_means_data_id, 
       iter = 2000, 
       chains = 4,
       cores = 4,
       seed = 17,
       control = list(adapt_delta = 0.9, 
                      max_treedepth=11))
})

saveRDS(uncond_means_stan_id, "/p/projects/open/ledig/thesis/fits/stan_uncond_means_id.rds")

#Unconditional Growth Model ID####################################################################################
#aims to look at general effects across cells and over time prior to adding other preds

uncond_growth_data_id <- 
  make_standata(data = aei_small_std,
                formula = bf(irrcrop ~ 1 + yearcount + (1 + yearcount|id),
                             coi ~ 1 + yearcount + (1 + yearcount|id),
                             zoi ~ 1 + yearcount + (1 + yearcount|id)),
                family = zero_one_inflated_beta())

uncond_growth_code_id <-
  make_stancode(data = aei_small_std,
                formula = bf(irrcrop ~ 1 + yearcount + (1 + yearcount|id),
                             coi ~ 1 + yearcount + (1 + yearcount|id),
                             zoi ~ 1 + yearcount + (1 + yearcount|id)),
                family = zero_one_inflated_beta(),
                prior = c(prior(normal(0, 100), class = b)))

    

saveRDS(uncond_growth_stan_id, "/p/projects/open/ledig/thesis/fits/stan_uncond_growth_id.rds")

#Full Simple model ID####################################################################################

full_simple_data_id <- 
  make_standata(data = aei_small_std,
                formula = bf(irrcrop ~ 0 + Intercept + yearcount + ISO + ISO:yearcount + (1 + yearcount|id),
                                       coi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount + (1 + yearcount|id),
                                       zoi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount + (1 + yearcount|id)),
                family = zero_one_inflated_beta())

full_simple_code_id <-
  make_stancode(data = aei_small_std,
                formula =  bf(irrcrop ~ 0 + Intercept + yearcount + ISO + ISO:yearcount + (1 + yearcount|id),
                                       coi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount + (1 + yearcount|id),
                                       zoi ~ 0 + Intercept + yearcount + ISO + ISO:yearcount + (1 + yearcount|id)),
                family = zero_one_inflated_beta(),
                prior = c(prior(normal(0, 100), class = b)))

job::job({
full_simple_stan_id <-
  stan(model_code = full_simple_code_id,
       data = full_simple_data_id,
       iter = 2000, 
       chains = 4,
       cores = 4,
       seed = 17,
       control = list(adapt_delta = 0.9, 
                      max_treedepth=11))
})

saveRDS(full_simple_stan_id, "/p/projects/open/ledig/thesis/fits/stan_full_simple_id.rds")



#Unconditional Means Model ISO####################################################################################



uncond_means_data_ISO <- 
  make_standata(data = aei_small_std,
                formula = bf(irrcrop ~ 1 + (1|ISO),
                             coi ~ 1 + (1|ISO),
                             zoi ~ 1 + (1|ISO)),
                family = zero_one_inflated_beta())

uncond_means_code_ISO <-
  make_stancode(data = aei_small_std,
                formula = bf(irrcrop ~ 1 + (1|ISO),
                             coi ~ 1 + (1|ISO),
                             zoi ~ 1 + (1|ISO)),
                family = zero_one_inflated_beta(),
                sample_prior = "only")

uncond_means_stan_ISO <-
  stan(model_code = uncond_means_code_ISO, 
       data = uncond_means_data_ISO, 
       iter = 2000, 
       chains = 4,
       cores = 4,
       seed = 17,
       control = list(adapt_delta = 0.9, 
                      max_treedepth=11))

saveRDS(uncond_means_stan_id, "/p/projects/open/ledig/thesis/fits/stan_uncond_means_ISO.rds")

#Unconditional Growth Model ISO####################################################################################



uncond_growth_data_ISO <- 
  make_standata(data = aei_small_std,
                formula = bf(irrcrop ~ 1 + yearcount + (1 + yearcount|ISO),
                             coi ~ 1 + yearcount + (1 + yearcount|ISO),
                             zoi ~ 1 + yearcount + (1 + yearcount|ISO)),
                family = zero_one_inflated_beta())

uncond_growth_code_ISO <-
  make_stancode(data = aei_small_std,
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

saveRDS(uncond_growth_stan_id, "/p/projects/open/ledig/thesis/fits/stan_uncond_growth_ISO.rds")

#Full Simple model ISO####################################################################################


full_simple_data_ISO <- 
  make_standata(data = aei_small_std,
                formula = bf(irrcrop ~ 0 + Intercept + yearcount + six_regions + six_regions:yearcount (1 + yearcount|ISO),
                             coi ~ 0 + Intercept + yearcount + six_regions + six_regions:yearcount (1 + yearcount|ISO),
                             zoi ~ 0 + Intercept + yearcount + six_regions + six_regions:yearcount (1 + yearcount|ISO)),
                family = zero_one_inflated_beta())

full_simple_code_ISO <-
  make_stancode(data = aei_small_std,
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

saveRDS(full_simple_stan_id, "/p/projects/open/ledig/thesis/fits/stan_full_simple_ISO.rds")



# Complex model####################################################################################

# without interactions, at lower levels they are easier to sort out but here its too complex

full_complex_data <- 
  make_standata(data = aei_small_std,
                formula = bf(irrcrop ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id),
                             coi ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id),
                             zoi ~ 0 + Intercept + yearcount + DD.regime + dist + popdens + pet + gdppc + rugged + (1 + yearcount + popdens + DD.regime + gdppc|id)),
                family = zero_one_inflated_beta())

full_complex_code <-
  make_stancode(data = aei_small_std,
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

saveRDS(full_complex_stan, "/p/projects/open/ledig/thesis/fits/stan_full_complex.rds")




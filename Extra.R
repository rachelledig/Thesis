## Packages
library(dplyr)
library(rstan)

## Individual cell change model

#load unstandardized data
aei_std <- read.csv("/Volumes/RachelExternal/Thesis/Data_and_Plots/aei_std.csv")

#nest data
by_country <-
  aei_std %>%
  group_by(Countryname) %>%
  nest()

#run first model based on time series
// generated with brms 2.15.0
code1 <- "functions {
  /* zero-one-inflated beta log-PDF of a single response 
  * Args: 
    *   y: response value 
  *   mu: mean parameter of the beta part
  *   phi: precision parameter of the beta part
  *   zoi: zero-one-inflation probability
  *   coi: conditional one-inflation probability
  * Returns:  
    *   a scalar to be added to the log posterior 
  */ 
    real zero_one_inflated_beta_lpdf(real y, real mu, real phi,
                                     real zoi, real coi) {
      row_vector[2] shape = [mu * phi, (1 - mu) * phi]; 
      if (y == 0) { 
        return bernoulli_lpmf(1 | zoi) + bernoulli_lpmf(0 | coi); 
      } else if (y == 1) {
        return bernoulli_lpmf(1 | zoi) + bernoulli_lpmf(1 | coi);
      } else { 
        return bernoulli_lpmf(0 | zoi) + beta_lpdf(y | shape[1], shape[2]);
      } 
    }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_zoi;  // number of population-level effects
  matrix[N, K_zoi] X_zoi;  // population-level design matrix
  int<lower=1> K_coi;  // number of population-level effects
  matrix[N, K_coi] X_coi;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  int Kc_zoi = K_zoi - 1;
  matrix[N, Kc_zoi] Xc_zoi;  // centered version of X_zoi without an intercept
  vector[Kc_zoi] means_X_zoi;  // column means of X_zoi before centering
  int Kc_coi = K_coi - 1;
  matrix[N, Kc_coi] Xc_coi;  // centered version of X_coi without an intercept
  vector[Kc_coi] means_X_coi;  // column means of X_coi before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 2:K_zoi) {
    means_X_zoi[i - 1] = mean(X_zoi[, i]);
    Xc_zoi[, i - 1] = X_zoi[, i] - means_X_zoi[i - 1];
  }
  for (i in 2:K_coi) {
    means_X_coi[i - 1] = mean(X_coi[, i]);
    Xc_coi[, i - 1] = X_coi[, i] - means_X_coi[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> phi;  // precision parameter
  vector[Kc_zoi] b_zoi;  // population-level effects
  real Intercept_zoi;  // temporary intercept for centered predictors
  vector[Kc_coi] b_coi;  // population-level effects
  real Intercept_coi;  // temporary intercept for centered predictors
}
transformed parameters {
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    // initialize linear predictor term
    vector[N] zoi = Intercept_zoi + Xc_zoi * b_zoi;
    // initialize linear predictor term
    vector[N] coi = Intercept_coi + Xc_coi * b_coi;
    for (n in 1:N) {
      // apply the inverse link function
      zoi[n] = inv_logit(zoi[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      coi[n] = inv_logit(coi[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      mu[n] = inv_logit(mu[n]);
    }
    for (n in 1:N) {
      target += zero_one_inflated_beta_lpdf(Y[n] | mu[n], phi, zoi[n], coi[n]);
    }
  }
  // priors including constants
  target += normal_lpdf(b | 0, 100);
  target += student_t_lpdf(Intercept | 3, 0, 2.5);
  target += gamma_lpdf(phi | 0.01, 0.01);
  target += logistic_lpdf(Intercept_zoi | 0, 1);
  target += logistic_lpdf(Intercept_coi | 0, 1);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_zoi_Intercept = Intercept_zoi - dot_product(means_X_zoi, b_zoi);
  // actual population-level intercept
  real b_coi_Intercept = Intercept_coi - dot_product(means_X_coi, b_coi);
}"

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

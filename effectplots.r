# "In Bayesian analyses, the key to your inference is the parameter of interest’s
# posterior distribution. It fulfills every property of a probability distribution
# and quantifies how probable it is for the population parameter to lie in certain
# regions. On the one hand, you can characterize the posterior by its mode. This is
# the parameter value that, given the data and its prior probability, is most 
# probable in the population. Alternatively, you can use the posterior’s mean or
# median. Using the same distribution, you can construct a 95% credibility interval,
# the counterpart to the confidence interval in frequentist statistics. Other than
# the confidence interval, the Bayesian counterpart directly quantifies the 
# probability that the population value lies within certain limits. There is a 95% 
# probability that the parameter value of interest lies within the boundaries of the
# 95% credibility interval. Unlike the confidence interval, this is not merely a 
# simulation quantity, but a concise and intuitive probability statement. For more 
# on how to interpret Bayesian analysis, check Van de Schoot et al. 2014."

library(tidyverse)
library(brms)
library(tidybayes)
library(schoenberg)
library(ape)

#
# This sets a global theme for all my plots. 
theme_set(theme_bw() +
            theme(
              plot.background = element_blank()
              ,panel.grid.major = element_blank()
              ,panel.grid.minor = element_blank()
              ,panel.background = element_blank()
              ,axis.text.x  = element_text(angle=90, vjust=0.5, size=8)
            ))

# Points at my working folder.
setwd("~/Documents/BEMP/bemp_to_r/litterfall/")

# Data cleaning for modeling. 
flood_model <- read_rds("./data/interim/litterfall_annual_data_for_modeling.RDS")
head(flood_model)

# Flood model

# Flipping the depth to groundwater to be negative numbers. 
flood_model <-flood_model %>% 
  mutate(mean_depth_gw_neg = -1 * mean_annual_depth_gw )
head(flood_model)

rowSums(is.na(flood_model))

summary(flood_model)

clean_bemp_data <- flood_model[complete.cases(flood_model), ]
head(clean_bemp_data)

unique(clean_bemp_data$site)
summary(clean_bemp_data)

# Check for weird spatial magic. Shouldn't be an issue. 
dist = with(clean_bemp_data, dist(cbind(lat, lon), diag = T, upper = T))
inv_dist = as.matrix(1/dist)
str(inv_dist)
inv_dist[which(!is.finite(inv_dist))] <- 0
inv_dist
ape::Moran.I(clean_bemp_data$cw, weight = inv_dist, scaled=T, na.rm=TRUE)

# Assumptions of the linear models from Data Analysis Using Regression and 
# Multilevel/Hierarchical Models, 2002, Gelman and Hill. 

# 1. Validity. Most importantly, the data you are analyzing should map to the
# research question you are trying to answer. 

# 2. Additivity and linearity. The most important mathematical assumption of 
# the regression model is that its deterministic component is a linear function
# of the separate predictors.
 
# 3. Independence of errors.
 
# 4. Equal variance of errors.
 
# 5. Normality of errors.

# Sets up the linear models to be combined into a piecewise structural equation model.
cw_mod <- bf(cw ~ fire_impact + flood_impact + cleared_impact + mean_depth_gw_neg +
               (1|site) + s(year, bs='gp', m=1) + (1|age))
will_mod <- bf(will ~ fire_impact + flood_impact + cleared_impact + mean_depth_gw_neg + 
                 (1|site) + s(year, bs='gp', m=1))
nmol_mod <- bf(nmol ~ fire_impact + flood_impact + cleared_impact + mean_depth_gw_neg + 
                 (1|site) + s(year, bs='gp', m=1))
ro_mod <- bf(ro ~ fire_impact + flood_impact + cleared_impact + mean_depth_gw_neg + 
               (1|site) + s(year, bs='gp', m=1))
sc_mod <- bf(sc ~ fire_impact + flood_impact + cleared_impact  + mean_depth_gw_neg + 
               (1|site) + s(year, bs='gp', m=1))
woody_mod <- bf(woody ~ fire_impact + flood_impact + cleared_impact + mean_depth_gw_neg + 
                  (1|site) + s(year, bs='gp', m=1))

# See what brms suggest for default priors (weakly informed)
get_prior(cw_mod + will_mod + sc_mod + ro_mod + woody_mod + nmol_mod, data = clean_bemp_data)

# Priors used in this analysis
prior2 <- c(set_prior("normal(-1, 1)", class = "b", coef= "fire_impact", resp="cw"),
            set_prior("normal(0, 1)", class = "b", coef= "flood_impact", resp="cw"),
            set_prior("normal(0, 1)", class = "b", coef= "cleared_impact", resp="cw"),
            set_prior("normal(0.5, 0.5)", class = "b", coef= "mean_depth_gw_neg",
                      resp="cw"),
            set_prior("normal(-1, 1)", class = "b", coef= "fire_impact", resp="will"),
            set_prior("normal(0.5, 0.5)", class = "b", coef= "flood_impact", resp="will"),
            set_prior("normal(0, 1)", class = "b", coef= "cleared_impact", resp="will"),
            set_prior("normal(0.5, 0.5)", class = "b", coef= "mean_depth_gw_neg",
                      resp="will"),
            set_prior("normal(-1, 1)", class = "b", coef= "fire_impact", resp="sc"),
            set_prior("normal(0, 1)", class = "b", coef= "flood_impact", resp="sc"),
            set_prior("normal(-1, 0.5)", class = "b", coef= "cleared_impact", resp="sc"),
            set_prior("normal(0, 1)", class = "b", coef= "mean_depth_gw_neg",
                      resp="sc"),
            set_prior("normal(-1, 1)", class = "b", coef= "fire_impact", resp="ro"),
            set_prior("normal(0, 1)", class = "b", coef= "flood_impact", resp="ro"),
            set_prior("normal(-1, 0.5)", class = "b", coef= "cleared_impact", resp="ro"),
            set_prior("normal(0, 1)", class = "b", coef= "mean_depth_gw_neg",
                      resp="ro"),
            set_prior("normal(0, 1)", class = "b", coef= "fire_impact", resp="nmol"),
            set_prior("normal(0, 1)", class = "b", coef= "flood_impact", resp="nmol"),
            set_prior("normal(0, 1)", class = "b", coef= "cleared_impact", resp="nmol"),
            set_prior("normal(0.5, 0.5)", class = "b", coef= "mean_depth_gw_neg",
                      resp="nmol"),
            set_prior("normal(1, 0.5)", class = "b", coef= "fire_impact", resp="woody"),
            set_prior("normal(-1, 0.5)", class = "b", coef= "flood_impact", resp="woody"),
            set_prior("normal(0, 1)", class = "b", coef= "cleared_impact", resp="woody"),
            set_prior("normal(-1, 0.5)", class = "b", coef= "mean_depth_gw_neg",
                      resp="woody"))

k_fit_brms <- brm(cw_mod + will_mod + sc_mod + ro_mod + woody_mod + nmol_mod +
                    set_rescor(TRUE), 
                  data=clean_bemp_data, prior = prior2, family = "gaussian",
                  cores=6, chains = 4, iter = 2500, 
                  control = list(adapt_delta=0.99, max_treedepth = 12))

### Model summary at the 50% uncertainity interval
# Group summary and population level estimates
summary(k_fit_brms, prob = 0.5)

# Smooth Terms from the summary)
# Estimate is how wiggly the smoothing line is. Higher values, more wiggle.
# If your 50% uncertainity interval does not include 0 then you can take that
# to mean there is evidence that a smoothing term is recommended. 

ms_k_fit_grms <- marginal_smooths(k_fit_brms)
plot(ms_k_fit_grms)

# Full group level estimates
ranef(k_fit_brms, groups="site", probs = 0.5)

### Model checking plots. Posterior predict. 
pp_check(k_fit_brms, resp = "cw", nsamples = 100)
pp_check(k_fit_brms, resp = "will", nsamples = 100)
pp_check(k_fit_brms, resp = "nmol", nsamples = 100)
pp_check(k_fit_brms, resp = "sc", nsamples = 100)
pp_check(k_fit_brms, resp = "ro", nsamples = 100)
pp_check(k_fit_brms, resp = "woody", nsamples = 100)

###
pp_check(k_fit_brms, resp = "cw", type = "ecdf_overlay")


# Where is the variance in the τ parameters? (group level)
posterior_samples(k_fit_brms) %>% 
  select(starts_with("sd_site__")) %>% 
  gather(key, tau) %>% 
  mutate(key = str_remove(key, "sd_site__") %>% str_remove(., "__Intercept")) %>% 

  ggplot(aes(x = tau, fill = key)) +
  geom_density(color = "transparent", alpha = 2/3) +
  scale_fill_viridis_d(NULL, end = .85) +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab(expression(tau)) +
  theme(panel.grid = element_blank())

posterior_samples(k_fit_brms) %>% 
  select(starts_with("sd_year")) %>% 
  gather(key, tau) %>% 
  mutate(key = str_remove(key, "sd_") %>% str_remove(., "__Intercept")) %>% 
  
  ggplot(aes(x = tau, fill = key)) +
  geom_density(color = "transparent", alpha = 2/3) +
  scale_fill_viridis_d(NULL, end = .85) +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab(expression(tau)) +
  theme(panel.grid = element_blank())

posterior_samples(k_fit_brms) %>% 
  select(starts_with("sd_age")) %>% 
  gather(key, tau) %>% 
  mutate(key = str_remove(key, "sd_") %>% str_remove(., "__Intercept")) %>% 
  
  ggplot(aes(x = tau, fill = key)) +
  geom_density(color = "transparent", alpha = 2/3) +
  scale_fill_viridis_d(NULL, end = .85) +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab(expression(tau)) +
  theme(panel.grid = element_blank())

### Marginal effects plots

# Cottonwood
p_depth_cw <- marginal_effects(k_fit_brms, "mean_depth_gw_neg", resp = "cw")
plot(p_depth_cw)[[1]] + ylim(0, 350)

# Willow
p_depth_will <- marginal_effects(k_fit_brms, "cleared_impact", resp = "will")
plot(p_depth_will)[[1]] + ylim(-5, 25)

# Saltcedar
p_depth_sc <- marginal_effects(k_fit_brms, "mean_depth_gw_neg", resp = "sc")
plot(p_depth_sc)[[1]] + ylim(-10, 30)

# Russian olive
p_depth_ro <- marginal_effects(k_fit_brms, "mean_depth_gw_neg", resp = "ro")
plot(p_depth_ro)[[1]] + ylim(-10, 30)

# Woody
p_depth_woody <- marginal_effects(k_fit_brms, "mean_depth_gw_neg", resp = "woody")
plot(p_depth_woody)[[1]] + ylim(-10, 75)

### New data give the old data and model. Does it predict back reasonable values?
# Solid line is the mean and the dashed is the 50% uncertainity interval.
# Cottonwood
k_fit_brms %>%
  spread_draws(b_cw_Intercept, r_site__cw[site,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = b_cw_Intercept + r_site__cw) %>%
  ungroup() %>%
  mutate(site = str_replace_all(site, "[.]", " ")) %>% 
  
  # plot
  ggplot(aes(x = mu, y = reorder(site, mu))) +
  geom_vline(xintercept = fixef(k_fit_brms)[1, 1], color = "#839496", size = 1) +
  geom_vline(xintercept = fixef(k_fit_brms)[1, 3:4], color = "#839496", linetype = 2) +
  geom_halfeyeh(.width = .5, size = 2/3, fill = "#859900") +
  labs(x = expression("Cottonwood litterfall (g/m^2)"),
       y = "BEMP sites ordered by mean predicted litterfall") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        text = element_text(family = "Ubuntu")) 

### Solid line is the mean and the dashed is the 50% uncertainity interval.
# Willow
k_fit_brms %>%
  spread_draws(b_will_Intercept, r_site__will[site,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = b_will_Intercept + r_site__will) %>%
  ungroup() %>%
  mutate(site = str_replace_all(site, "[.]", " ")) %>% 
  
  # plot
  ggplot(aes(x = mu, y = reorder(site, mu))) +
  geom_vline(xintercept = fixef(k_fit_brms)[2, 1], color = "#839496", size = 1) +
  geom_vline(xintercept = fixef(k_fit_brms)[2, 3:4], color = "#839496", linetype = 2) +
  geom_halfeyeh(.width = .5, size = 2/3, fill = "#859900") +
  labs(x = expression("Willow litterfall (g/m^2)"),
       y = NULL) +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        text = element_text(family = "Ubuntu"))

### Solid line is the mean and the dashed is the 50% uncertainity interval.
# Saltcedar
k_fit_brms %>%
  spread_draws(b_sc_Intercept, r_site__sc[site,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = b_sc_Intercept + r_site__sc) %>%
  ungroup() %>%
  mutate(site = str_replace_all(site, "[.]", " ")) %>% 
  
  # plot
  ggplot(aes(x = mu, y = reorder(site, mu))) +
  geom_vline(xintercept = fixef(k_fit_brms)[3, 1], color = "#839496", size = 1) +
  geom_vline(xintercept = fixef(k_fit_brms)[3, 3:4], color = "#839496", linetype = 2) +
  geom_halfeyeh(.width = .5, size = 2/3, fill = "#6c71c4") +
  labs(x = expression(italic("Saltcedar litterfall (g/m^2)")),
       y = NULL) +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        text = element_text(family = "Ubuntu")) 

### Solid line is the mean and the dashed is the 50% uncertainity interval.
# Russian olive
k_fit_brms %>%
  spread_draws(b_ro_Intercept, r_site__ro[site,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = b_ro_Intercept + r_site__ro) %>%
  ungroup() %>%
  mutate(site = str_replace_all(site, "[.]", " ")) %>% 
  
  # plot
  ggplot(aes(x = mu, y = reorder(site, mu))) +
  geom_vline(xintercept = fixef(k_fit_brms)[4, 1], color = "#839496", size = 1) +
  geom_vline(xintercept = fixef(k_fit_brms)[4, 3:4], color = "#839496", linetype = 2) +
  geom_halfeyeh(.width = .5, size = 2/3, fill = "#6c71c4") +
  labs(x = expression("Russian olive litterfall (g/m^2)"),
       y = NULL) +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        text = element_text(family = "Ubuntu")) 

# Model forward a 50cm drop on average assuming var. is constant. 

# Model forward a 50cm drop on average assuming var. is increasing.
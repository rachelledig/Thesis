

library(tidyverse)

library(tidybayes)
library(brms)

#
# This sets a global theme for all my plots. 
theme_set(theme_minimal() +
            theme(
              plot.background = element_blank()
              ,panel.grid.major = element_blank()
              ,panel.background = element_blank()
              ,axis.text.x  = element_text(angle=90, vjust=0.5, size=8)
            ))


### Model summary at the 50% uncertainity interval
# Group summary and population level estimates
summary(composite_reg, prob = 0.5)

#fixed effects
fixef(composite_reg, groups="eight_regions", probs = 0.5)

#random effects
ranef(composite_reg, groups="eight_regions", probs = 0.5)

#combined
coef(composite_reg)$eight_regions

### Model checking plots. Posterior predict. 
pp_check(composite_reg, ndraws = 100)

###
pp_check(composite_ISO, type = "ecdf_overlay", ndraws = 100)


### Marginal effects plots mu
conds <- make_conditions(data, vars = "eight_regions")

# yearcount:precip
p_yearcount_precip <- conditional_effects(composite_reg, "yearcount:precip", dpar = "mu", conditions = conds, re_formula = NULL)
plot(p_yearcount_precip)[[1]] 

# yearcount:rugged
p_yearcount_rugged <- conditional_effects(composite_reg, "yearcount:rugged", dpar = "mu", conditions = conds, re_formula = NULL)
plot(p_yearcount_rugged)[[1]]

# yearcount:gdppc
p_yearcount_gdppc <- conditional_effects(composite_reg, "yearcount:gdppc", dpar = "mu", conditions = conds, re_formula = NULL)
plot(p_yearcount_gdppc)[[1]] 

# yearcount:popdens
p_yearcount_popdens <- conditional_effects(composite_reg, "yearcount:popdens", dpar = "mu", conditions = conds, re_formula = NULL)
plot(p_yearcount_popdens)[[1]] 

# yearcount:DD.regime
p_yearcount_DD_regime <- conditional_effects(composite_reg, "yearcount:DD.regime", dpar = "mu", conditions = conds, re_formula = NULL)
plot(p_yearcount_DD_regime)[[1]] 



### Marginal effects plots zi
conds <- make_conditions(data, vars = c("ISO", "eight_regions")) 

# rugged
p_zi_rugged <- conditional_effects(`165aa_5yronlyISO_notimeserieszi`, "rugged", dpar = "zi", conditions = conds, re_formula = NULL)
plot(p_zi_rugged)[[1]] 

# precip
p_zi_precip <- conditional_effects(`165aa_5yronlyISO_notimeserieszi`, "precip", dpar = "zi", conditions = conds, re_formula = NULL)
plot(p_zi_precip)[[1]]

# dist
p_zi_dist <- conditional_effects(`165aa_5yronlyISO_notimeserieszi`, "dist", dpar = "zi", conditions = conds, re_formula = NULL)
plot(p_zi_dist)[[1]] 

# gdppc
p_zi_gdppc <- conditional_effects(`165aa_5yronlyISO_notimeserieszi`, "gdppc", dpar = "zi", conditions = conds, re_formula = NULL)
plot(p_zi_gdppc)[[1]] 


### New data give the old data and model. Does it predict back reasonable values?
# Solid line is the mean and the dashed is the 50% uncertainity interval.


#posterior draws to construct mu
composite_reg %>%
  spread_draws(b_Intercept, r_eight_regions[eight_regions,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = b_Intercept + r_eight_regions) %>%
  ungroup() %>%
  mutate(mu = exp(mu)) %>% 
  
  # plot
  ggplot(aes(x = mu, y = reorder(eight_regions, mu))) +
  geom_vline(xintercept = exp(fixef(composite_reg)[1, 1]), color = "#839496", size = 1) +
  geom_vline(xintercept = exp(fixef(composite_reg)[1, 3:4]), color = "#839496", linetype = 2) +
  geom_halfeyeh(.width = .5, size = 2/3, fill = "#859900") +
  xlim(0,0.10)


#posterior draws to construct zi
`165aa_5yronlyISO_notimeserieszi` %>%
  spread_draws(b_zi_Intercept, r_ISO__zi[ISO,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(zi = b_zi_Intercept + r_ISO__zi) %>%
  ungroup() %>%
  mutate(zi = exp(zi)) %>% 
  
  # plot
  ggplot(aes(x = zi, y = reorder(ISO, zi))) +
  geom_vline(xintercept = exp(fixef(`165aa_5yronlyISO_notimeserieszi`)[13, 1]), color = "#839496", size = 1) +
  geom_vline(xintercept = exp(fixef(`165aa_5yronlyISO_notimeserieszi`)[13, 3:4]), color = "#839496", linetype = 2) +
  stat_halfeye(.width = .5, size = 2/3, fill = "#859900") +
  xlim(0,1)

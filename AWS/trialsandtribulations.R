centered <-function(y) {
  x<-y[!is.na(y)]
  x<-(x - mean(x)) / (sd(x)*2)
  y[!is.na(y)]<-x
  return(y)
}

aei_nocrop_samp165aa_5yr <- aei_nocrop_samp165aa_5yr[,-1]
aei_nocrop_samp165aa_5yr <-
  aei_nocrop_samp165aa_5yr %>% 
  mutate(across(c(13:20), centered)) %>% 
  mutate(across(c(1:3,6,7,21,22), as.factor)) 


#This works - this is predicted - plots varying intercepts of an intercept only model
bound <- 
  aei_nocrop_samp165aa_5yr %>%
  data_grid(ISO) %>%
  add_epred_draws(`165aa_5yronlyISO_uncondmeans`) %>%
  group_by(ISO) %>%
  mutate(meanperclass = mean(.epred))%>%
  ungroup() 

aei_nocrop_samp165aa_5yr %>% 
  select(ISO, eight_regions) %>% 
  unique() %>% 
  left_join(bound, ., by = "ISO") %>%
  ggplot(aes(x = .epred, y = reorder(as.factor(ISO), meanperclass), color = eight_regions)) +
  stat_pointinterval(.width = c(.66, .95))



posterior_samples(`165aa_5yrfixed_effects_nodemo`) %>%
  select(starts_with(c("b_", "Intercept", "sd_", "phi"))) %>% 
  mutate(across(.cols = starts_with("phi"), exp)) %>% #hmmmmmm
  mutate(across(.cols = starts_with(c("b_", "Intercept")), plogis)) %>% 
  posterior_summary() %>% 
  round(digits = 3)



re_model_only <- crossing(yearcount = seq(min(aei_nocrop_samp165aa_5yr$yearcount), 
                                        max(aei_nocrop_samp165aa_5yr$yearcount), length.out=100),
                          ISO = unique(aei_nocrop_samp165aa_5yr$ISO)) %>%
  add_fitted_draws(`165aa_5yronlyISO_uncondgrowth`)

re_model_summary <- re_model_only %>%
  group_by(ISO, yearcount) %>%
  summarize(.value = mean(.value))

ggplot(re_model_only,
       aes(x = yearcount, y = .value)) +
  facet_wrap(~ISO) +
  stat_interval() 



aei_nocrop_samp165aa_5yr %>%
  subset(eight_regions == "south_america") %>% 
  group_by(ISO) %>%
  data_grid(yearcount = seq_range(yearcount, n = 45)) %>%
 add_epred_draws(`165aa_5yronlyISO_uncondgrowth`) %>%
  ggplot(aes(x = yearcount, y = irrfrac, color = ISO)) +
  stat_lineribbon(aes(y = .epred)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")



#####   THIS ONE WORKS - this is posterior samples
`165aa_5yronlyISO_uncondmeans` %>%
  +     spread_draws(b_Intercept, r_ISO[ISO,]) %>%
  +     median_qi(ISO_mean = b_Intercept + r_ISO, .width = c(.95, .66)) %>%
  +     ggplot(aes(y = ISO, x = ISO_mean, xmin = .lower, xmax = .upper)) +
  +     geom_pointinterval() 


### this one doesnt
`165aa_5yronlyISO_uncondmeans` %>%
  +     spread_rvars(`b_zi_Intercept`, r_ISO__zi[ISO,]) %>%
  +     median_qi(ISO_mean = b_zi_Intercept + r_ISO__zi, .width = c(.95, .66)) %>%
  +     ggplot(aes(y = ISO, x = ISO_mean, xmin = .lower, xmax = .upper)) +
  +     geom_pointinterval() 


bound2 <- 
  aei_nocrop_samp165aa_5yr %>%
  data_grid(ISO, yearcount) %>%
  add_epred_draws(`165aa_5yronlyISO_uncondgrowth`, dpar = "zi") %>%
  group_by(ISO, yearcount) %>%
  summarise(meanperclass = mean(zi))

aei_nocrop_samp165aa_5yr %>% 
  select(ISO, eight_regions) %>% 
  unique() %>% 
  left_join(bound2, ., by = "ISO") %>%
  ggplot(aes(x = yearcount, y = meanperclass, color = ISO, group = ISO)) +
  geom_line() +
  theme(legend.position = "none") +
  facet_wrap(~eight_regions) 


bound2 <- 
  aei_nocrop_samp165aa_5yr %>%
  add_epred_draws(`165aa_5yronlyISO_varying_effects_nodemo`, dpar = "mu") %>%
  group_by(ISO, yearcount) %>%
  summarise(meanperclass = mean(.epred))

aei_nocrop_samp165aa_5yr %>% 
  select(ISO, eight_regions) %>% 
  unique() %>% 
  left_join(bound2, ., by = "ISO") %>%
  ggplot(aes(x = yearcount, y = meanperclass, color = ISO, group = ISO)) +
  geom_line() +
  theme(legend.position = "none") +
  facet_wrap(~eight_regions) 


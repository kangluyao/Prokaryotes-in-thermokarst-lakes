# SEM
## load the package
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(nlme)

climate_vars <- env_div_agg[ , c('MAP')]
DOM_vars <- env_div_agg[ , c('FrI')]
physio_chemical_vars <- env_div_agg[ , c('Salinity')]
ALT <- env_div_agg[ , c('ALT')]

library(piecewiseSEM)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI + Salinity +ALT, data = env_div_agg),
  lm(ALT ~ MAP, data = env_div_agg),
  lm(FrI ~ MAP + ALT, data = env_div_agg),
  lm(Salinity ~ MAP + ALT + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI + Salinity +ALT, data = env_div_agg),
  lm(ALT ~ MAP, data = env_div_agg),
  lm(FrI ~ MAP + ALT, data = env_div_agg),
  lm(Salinity ~ MAP + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI + Salinity +ALT, data = env_div_agg),
  lm(FrI ~ MAP + ALT, data = env_div_agg),
  lm(Salinity ~ MAP + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI + Salinity, data = env_div_agg),
  lm(FrI ~ MAP + ALT, data = env_div_agg),
  lm(Salinity ~ MAP + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI, data = env_div_agg),
  lm(FrI ~ MAP + ALT, data = env_div_agg),
  lm(Salinity ~ MAP + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI, data = env_div_agg),
  lm(FrI ~ MAP, data = env_div_agg),
  lm(Salinity ~ MAP + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI, data = env_div_agg),
  lm(FrI ~ MAP, data = env_div_agg),
  lm(Salinity ~ MAP, data = env_div_agg))
summary(sem_mod)
fisherC(sem_mod)
rsquared(sem_mod)

coeffs_sem <- coefs(sem_mod, standardize = "scale")[,c(1,2,8)]
coeffs_sem
## Direct effect of climate variables on total community composition.
clim_dir_effect <- coeffs_sem[1,3]
## Ind effect of climate variables on total community composition.
clim_ind_effect <- coeffs_sem[3,3]*coeffs_sem[2,3]

total_clim_effect <- clim_dir_effect + clim_ind_effect

## Direct effect of ALT on total community composition.
ALT_dir_effect <- 0
## Ind effect of ALT on total community composition.
ALT_ind_effect <- 0
total_ALT_effect <- ALT_dir_effect + ALT_ind_effect

## Direct effect of DOM variables on total community composition.
DOM_dir_effect <- coeffs_sem[2,3]
## Ind effect of DOM variables on total community composition.
DOM_ind_effect <- 0
total_DOM_effect <- DOM_dir_effect + DOM_ind_effect

## Direct effect of physio_chemical variables on total community composition.
physio_chemical_dir_effect <- 0
## Ind effect of climate variables on total community composition.
physio_chemical_ind_effect <- 0
total_physio_chemical_effect <- physio_chemical_dir_effect + physio_chemical_ind_effect
## bar plot
sem_effect <- data.frame(variable = c(rep('Climate elements',3),
                                      rep('ALT', 3),
                                      rep('DOM properties', 3),
                                      rep('physio-chemical factors', 3)),
                         type_of_effects = c(rep(c('Total effects', 'Direct effects', 
                                                   'Indirect effects'), 4)),
                         value = c(total_clim_effect, clim_dir_effect, clim_ind_effect,
                                   total_ALT_effect, ALT_dir_effect, ALT_ind_effect,
                                   total_DOM_effect, DOM_dir_effect, DOM_ind_effect,
                                   total_physio_chemical_effect, physio_chemical_dir_effect,
                                   physio_chemical_ind_effect))
sem_effect$type_of_effects <- factor(sem_effect$type_of_effects, ordered = T,
                                     levels = c('Total effects', 'Direct effects', 'Indirect effects'))
sem_effect$variable <- factor(sem_effect$variable, ordered = T,
                              levels = c('Climate elements', 'DOM properties', 'physio-chemical factors', 'ALT'))
sem_effect_plot <-  ggplot(data = sem_effect, 
                           aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat="identity", width = 0.6) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-0.8, 0.5)) +
  labs(x = 'Variables', y = 'Standardized effects from SEM') +
  facet_wrap(~ type_of_effects) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(color='black',size=14),
        strip.text = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

sem_effect_plot
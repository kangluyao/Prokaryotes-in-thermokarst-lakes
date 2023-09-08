# Loading packages
library(ggplot2)
library(tidyverse)
library(adespatial)
#####################################################
## Determine the LCBD of total community
beta_meta_div <- beta.div(t(as.matrix(otu_table(meta_physeq))), 
                          method = "hellinger", sqrt.D = FALSE, samp = TRUE, 
                          nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)

env_div <- data.frame(LCBD = beta_meta_div$LCBD, sample_data(meta_physeq))

# random forest analysis for Northern Hemisphere
## Regression:
library(randomForest)
library(rfPermute)
library(A3)
env_div_rf <- env_div %>%
  dplyr::select(c('latitude', 'longitude', 'LCBD',  
                  'MAT', 'MAP', 'DOC', 'SUVA254', 
                  'a320', 'pH'))
set.seed(123)
out.rf <- a3(LCBD ~ . + 0, data = (env_div_rf)[3:9], randomForest, 
             model.args = list(ntree = 999, num.rep = 999), p.acc = 0.001) # It took about 2 hours on a 12-core laptop 
out.rf

meta.lcbd.rf <- rfPermute(LCBD ~ ., data = (env_div_rf)[3:9], ntree = 999, num.rep = 999,
                          importance = TRUE, na.action = na.omit)
## Show "importance" of variables: higher value mean more important:
mean(meta.lcbd.rf$rf$rsq)
sd(meta.lcbd.rf$rf$rsq)
meta.lcbd.rf$pval

impor.dat <-data.frame(variables = rownames(importance(meta.lcbd.rf)), 
                       IncMSE = importance(meta.lcbd.rf)[,1])
library(dplyr)
impor_meta_plot <- impor.dat %>% 
  arrange(IncMSE)  %>%
  mutate(variables = factor(variables, levels=variables)) %>%
  ggplot(aes(x = variables, y = IncMSE)) +
  geom_bar(stat = "identity", fill = rev(c("#da5724", '#3b86bc', '#915c83', '#eb8f70', '#c6dfed', '#f4c2c2')),
           colour = NA, width = 0.50) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 63)) +
  ylab('Increases in MSE (%)') + xlab(NULL) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", size=0.5, fill=NA),
        panel.grid = element_blank(), 
        axis.title = element_text(color='black',size = 14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color = 'black'),
        axis.text.y = element_text(colour='black',size = 12),
        axis.text.x = element_text(colour='black', size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.6, 0.75)) +
  coord_flip()
impor_meta_plot
# test the relationship between LCBD and MAP across the Northern Hemisphere
meta_lcbd_map <- ggplot(env_div, aes(MAP, LCBD)) +
  geom_point(aes(color = Sitegroup), size = 1, alpha = 0.8)+
  scale_color_manual(values = c("#a58aff", '#c49a00', '#c6dfed', 
                                '#00c094', '#3b86bc', '#eb8f70', '#fb61d7')) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2),
              size = 1, se = T, colour = 'black') +
  xlab('MAP (mm)') + ylab('LCBD') +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", size = 0.5, fill = NA),
        panel.grid = element_blank(), 
        axis.title = element_text(color = 'black',size = 14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color = 'black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.text.x = element_text(colour = 'black', size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.6, 0.75))
# check the spatial autocorrelation
## determine the distance matrix
env_div_rf <- na.omit(env_div_rf)
# sum(apply(scale(env_div_rf), 2, is.nan))
# apply(env_div_rf, 2, var) == 0
library(geodist)
x <- tibble::tibble (x = env_div_rf$longitude,
                     y = env_div_rf$latitude)
distance_matrix <- geodist (x)/100000

#names of the response variable and the predictors
dependent.variable.name <- "LCBD"
predictor.variable.names <- colnames(env_div_rf)[4:9]

#coordinates of the cases
xy <- env_div_rf[, c('longitude', 'latitude')]

#distance matrix
distance.matrix <- distance_matrix

#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 10, 50, 100, 150, 180)

#random seed for reproducibility
random.seed <- 1
#Fitting a non-spatial Random Forest model with rf()
model.spatial <- spatialRF::rf_spatial(
  data = env_div_rf,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  verbose = FALSE
)
#shows the Moran’s I of the residuals of the spatial model
Moran_residual_plot_NH <- spatialRF::plot_moran(model.spatial, verbose = FALSE)

#############################################
# random forest analysis for Tibetan Plateau
env_div_rf_tp <- env_div %>%
  filter(Region == 'Tibetan Plateau') %>%
  dplyr::select(c('latitude', 'longitude','LCBD',  'MAT', 
                  'MAP', 'DOC', 'SUVA254', 'a320', 'pH'))
set.seed(123)
out.rf.tp <- a3(LCBD ~ . + 0, data = (env_div_rf_tp)[3:9], randomForest, 
                model.args = list(ntree = 999, num.rep = 999), p.acc = 0.001)
out.rf.tp
tp.lcbd.rf <- rfPermute(LCBD ~ ., data = (env_div_rf_tp)[3:9], ntree = 999, num.rep = 999,
                     importance = TRUE, na.action = na.omit)
## Show "importance" of variables: higher value mean more important:
mean(tp.lcbd.rf$rf$rsq)
sd(tp.lcbd.rf$rf$rsq)
tp.lcbd.rf$pval

impor.dat <- data.frame(variables = rownames(importance(tp.lcbd.rf)), 
                       IncMSE = importance(tp.lcbd.rf)[,1])

library(dplyr)
impor_tp_plot <- impor.dat %>% 
  arrange(IncMSE)  %>%
  mutate(variables=factor(variables, levels = variables)) %>%
  ggplot(aes(x = variables, y = IncMSE)) +
  geom_bar(stat = "identity", fill = rev(c("#da5724", '#3b86bc', '#c6dfed', '#eb8f70', '#f4c2c2', '#915c83')),
           colour = NA, width = 0.50) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
  ylab('Increases in MSE (%)') + xlab(NULL) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", size=0.5, fill=NA),
        panel.grid = element_blank(), 
        axis.title = element_text(color='black',size = 14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color = 'black'),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.text.x = element_text(colour = 'black', size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.6, 0.75)) +
  coord_flip()
# check the spatial autocorrelation
#determine the distance matrix
library(geodist)
x <- tibble::tibble (x = env_div_rf$longitude,
                     y = env_div_rf$latitude)
distance_matrix <- geodist (x)/100000

#names of the response variable and the predictors
dependent.variable.name <- "LCBD"
predictor.variable.names <- colnames(env_div_rf)[4:9]

#coordinates of the cases
xy <- env_div_rf[, c('longitude', 'latitude')]

#distance matrix
distance.matrix <- distance_matrix

#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 1, 2, 4, 8, 11)

#random seed for reproducibility
random.seed <- 1
#Fitting a non-spatial Random Forest model with rf()
model.spatial <- spatialRF::rf_spatial(
  data = env_div_rf,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  verbose = FALSE
)
#shows the Moran’s I of the residuals of the spatial model
Moran_residual_plot_TP <- spatialRF::plot_moran(model.spatial, verbose = FALSE)

#############################################
## random forest analysis for Pan-Arctic
env_div_rf_pa <- env_div %>%
  filter(Region == 'Pan-Arctic') %>%
  dplyr::select(c('latitude', 'longitude','LCBD',  'MAT', 
                  'MAP', 'DOC', 'SUVA254', 'a320', 'pH'))
set.seed(123)
out.rf.pa <- a3(LCBD ~ . + 0, data = (env_div_rf_pa)[3:9], randomForest, 
                model.args = list(ntree = 999, num.rep = 999), p.acc = 0.001)
out.rf.pa
pa.lcbd.rf <- rfPermute(LCBD ~ ., data = (env_div_rf_pa)[3:9], ntree = 999, num.rep = 999,
                     importance = TRUE, na.action = na.omit)
## Show "importance" of variables: higher value mean more important:
mean(pa.lcbd.rf$rf$rsq)
sd(pa.lcbd.rf$rf$rsq)
pa.lcbd.rf$pval

impor.dat <-data.frame(variables = rownames(importance(pa.lcbd.rf)), 
                       IncMSE = importance(pa.lcbd.rf)[,1])

library(dplyr)
impor_pa_plot <- impor.dat %>% 
  arrange(IncMSE)  %>%
  mutate(variables=factor(variables, levels = variables)) %>%
  ggplot(aes(x = variables, y = IncMSE)) +
  geom_bar(stat="identity", fill= rev(c('#eb8f70', "#da5724", '#3b86bc', '#915c83',  '#c6dfed', '#f4c2c2')),
           colour = NA, width = 0.50) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
  ylab('Increases in MSE (%)') + xlab(NULL) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", size = 0.5, fill = NA),
        panel.grid=element_blank(), 
        axis.title = element_text(color = 'black', size = 14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color = 'black'),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(colour = 'black', size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.6, 0.75)) +
  coord_flip()
# check the spatial autocorrelation
#determine the distance matrix
env_div_rf <- na.omit(env_div_rf)
library(geodist)
x <- tibble::tibble (x = env_div_rf$longitude,
                     y = env_div_rf$latitude)
distance_matrix <- geodist (x)/100000

#names of the response variable and the predictors
dependent.variable.name <- "LCBD"
predictor.variable.names <- colnames(env_div_rf)[4:9]

#coordinates of the cases
xy <- env_div_rf[, c('longitude', 'latitude')]

#distance matrix
distance.matrix <- distance_matrix

#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 10, 20, 40, 80, 130)

#random seed for reproducibility
random.seed <- 1
#Fitting a non-spatial Random Forest model with rf()
model.spatial <- spatialRF::rf_spatial(
  data = env_div_rf,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  verbose = FALSE
)
#shows the Moran’s I of the residuals of the spatial model
Moran_residual_plot_PA <- spatialRF::plot_moran(model.spatial, verbose = FALSE)

## combine plots
rf_plot <- cowplot::plot_grid(impor_meta_plot, impor_tp_plot, impor_pa_plot,
                              labels = c('A', 'B', 'C'), ncol = 3, 
                              label_x = .01, label_y = 1, 
                              hjust = 0, label_size = 14, align = "v")
rf_plot
#ggsave(rf_plot, file = "./meta_analysis/results/figs/rf_plot.pdf",
#       width = 13, height = 6, units = 'in', device='pdf', dpi=300)

Moran_residual_plot <- cowplot::plot_grid(Moran_residual_plot_NH, 
                                          Moran_residual_plot_TP, 
                                          Moran_residual_plot_PA,
                                          labels = c('A', 'B', 'C'), ncol = 1, 
                                          label_x = .01, label_y = 1, 
                                          hjust = 0, label_size = 14, align = "v")
Moran_residual_plot

# correlation between diversity and environmental variables
library(lme4)
library(lmerTest)
library(multcomp)
ggplot(meta_diversity, aes(x = MAT, y = Chao1)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = MAP, y = Chao1)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = DOC, y = Chao1)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = SUVA254, y = Chao1)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = a320, y = Chao1)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = pH, y = Chao1)) + geom_point(alpha = 0.3)

ggplot(meta_diversity, aes(x = MAT, y = Shannon)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = MAP, y = Shannon)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = DOC, y = Shannon)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = SUVA254, y = Shannon)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = a320, y = Shannon)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = pH, y = Shannon)) + geom_point(alpha = 0.3)

ggplot(meta_diversity, aes(x = MAT, y = Simpson)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = MAP, y = Simpson)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = DOC, y = Simpson)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = SUVA254, y = Simpson)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = a320, y = Simpson)) + geom_point(alpha = 0.3)
ggplot(meta_diversity, aes(x = pH, y = Simpson)) + geom_point(alpha = 0.3)


mode1 <- lmer(log(Chao1) ~ Region + (1|Site), meta_diversity)
summary(mode1)
MuMIn::r.squaredGLMM(mode1)
mode2 <- lmer(Shannon ~ Region + (1|Site), meta_diversity)
summary(mode2)
MuMIn::r.squaredGLMM(mode2)
mode3 <- lmer(log(Simpson) ~ Region + (1|Site), meta_diversity)
summary(mode3)
MuMIn::r.squaredGLMM(mode3)
ggplot(meta_diversity, aes(x=MAT, y=Simpson)) + geom_point(alpha = 0.3)
mode4 <- lmer(log(Chao1) ~ MAT + (1|Site), meta_diversity)
summary(mode4)
MuMIn::r.squaredGLMM(mode4)
mode5 <- lmer(Shannon ~ MAT + (1|Site), meta_diversity)
summary(mode5)
MuMIn::r.squaredGLMM(mode5)
mode6 <- lmer(log(Simpson) ~ MAT + (1|Site), meta_diversity)
summary(mode6)
MuMIn::r.squaredGLMM(mode6)
ggplot(meta_diversity, aes(x=MAP, y=Chao1)) + geom_point(alpha = 0.3)
mode7 <- lmer(log(Chao1) ~ MAP + (1|Site), meta_diversity)
summary(mode7)
MuMIn::r.squaredGLMM(mode7)
mode8 <- lmer(Shannon ~ MAP + (1|Site), meta_diversity)
summary(mode8)
MuMIn::r.squaredGLMM(mode8)
mode9 <- lmer(log(Simpson) ~ MAP + (1|Site), meta_diversity)
summary(mode9)
MuMIn::r.squaredGLMM(mode9)
ggplot(meta_diversity, aes(x=DOC, y=Chao1)) + geom_point(alpha = 0.3)
mode10 <- lmer(log(Chao1) ~ DOC + (1|Site), meta_diversity)
summary(mode10)
MuMIn::r.squaredGLMM(mode10)
mode11 <- lmer(Shannon ~ DOC + (1|Site), meta_diversity)
summary(mode11)
MuMIn::r.squaredGLMM(mode11)
mode12 <- lmer(log(Simpson) ~ DOC + (1|Site), meta_diversity)
summary(mode12)
MuMIn::r.squaredGLMM(mode12)
mode13 <- lmer(log(Chao1) ~ SUVA254 + (1|Site), meta_diversity)
summary(mode13)
MuMIn::r.squaredGLMM(mode13)
mode14 <- lmer(Shannon ~ SUVA254 + (1|Site), meta_diversity)
summary(mode14)
MuMIn::r.squaredGLMM(mode14)
mode15 <- lmer(log(Simpson) ~ SUVA254 + (1|Site), meta_diversity)
summary(mode15)
MuMIn::r.squaredGLMM(mode15)
mode16 <- lmer(log(Chao1) ~ a320 + (1|Site), meta_diversity)
summary(mode16)
MuMIn::r.squaredGLMM(mode16)
mode17 <- lmer(Shannon ~ a320 + (1|Site), meta_diversity)
summary(mode17)
MuMIn::r.squaredGLMM(mode17)
mode18 <- lmer(log(Simpson) ~ a320 + (1|Site), meta_diversity)
summary(mode18)
MuMIn::r.squaredGLMM(mode18)
mode19 <- lmer(log(Chao1) ~ pH + (1|Site), meta_diversity)
summary(mode19)
MuMIn::r.squaredGLMM(mode19)
mode20 <- lmer(Shannon ~ pH + (1|Site), meta_diversity)
summary(mode20)
MuMIn::r.squaredGLMM(mode20)
mode21 <- lmer(log(Simpson) ~ pH + (1|Site), meta_diversity)
summary(mode21)
MuMIn::r.squaredGLMM(mode21)


# test the relationship between LCBD and MAP across the Northern Hemisphere
library(adespatial)
## total community
beta_meta_div <- beta.div(t(as.matrix(otu_table(meta_physeq))), 
                          method = "hellinger", sqrt.D = FALSE, samp = TRUE, 
                          nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)

env_div <- data.frame(LCBD = beta_meta_div$LCBD, sample_data(meta_physeq))

## calculate the average LCBD for each site
library(dplyr)
env_div_agg_meta <-  env_div %>% 
  dplyr::select(c(1, 5, 6, 7, 8, 19)) %>%
  group_by(Sitegroup, Site, Sitegroup1) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE))

# ggplot(env_div, aes(x=pH, y=LCBD)) + geom_point(alpha = 0.3)
# lm1 <- lm(LCBD ~  pH + I(pH^2), data = env_div)
# summary(lm1)
# ggplot(env_div, aes(x=DOC, y=LCBD)) + geom_point(alpha = 0.3)
# lm2 <- lm(LCBD ~  DOC, data = env_div)
# summary(lm2)
# ggplot(env_div, aes(x=a320, y=LCBD)) + geom_point(alpha = 0.3)
# lm3 <- lm(LCBD ~  a320 + I(a320^2), data = env_div)
# summary(lm3)
# ggplot(env_div, aes(x=SUVA254, y=LCBD)) + geom_point(alpha = 0.3)
# lm4 <- lm(LCBD ~  SUVA254, data = env_div)
# summary(lm4)
# ggplot(env_div, aes(x=log(MAT+15), y=LCBD)) + geom_point(alpha = 0.3)
# lm5 <- lm(LCBD ~ log(MAT+15), data = env_div)
# summary(lm5)
# ggplot(env_div, aes(x=MAP, y=LCBD)) + geom_point(alpha = 0.3)

## random forest analysis for Northern Hemisphere
### Regression:
library(randomForest)
library(rfPermute)
env_div_rf <- env_div %>%
  dplyr::select(c('latitude', 'longitude', 'LCBD',  'MAT', 'MAP', 'DOC', 'SUVA254', 'a320', 'pH'))
set.seed(123)
meta.lcbd.rf <- rfPermute(LCBD ~ ., data = (env_div_rf)[3:9], ntree = 999, num.rep = 999,
                     importance = TRUE, na.action = na.omit)

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

#determine the distance matrix
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

## random forest analysis for Tibetan Plateau
env_div_rf <- env_div %>%
  filter(Region == 'Tibetan Plateau') %>%
  dplyr::select(c('latitude', 'longitude','LCBD',  'MAT', 'MAP', 'DOC', 'SUVA254', 'a320', 'pH'))
set.seed(123)
tp.lcbd.rf <- rfPermute(LCBD ~ ., data = (env_div_rf)[3:9], ntree = 999, num.rep = 999,
                     importance = TRUE, na.action = na.omit)

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
## random forest analysis for Pan-Arctic
env_div_rf <- env_div %>%
  filter(Region == 'Pan-Arctic') %>%
  dplyr::select(c('latitude', 'longitude', 'LCBD',  'MAT', 'MAP', 'DOC', 'SUVA254', 'a320', 'pH'))
set.seed(123)
pa.lcbd.rf <- rfPermute(LCBD ~ ., data = (env_div_rf)[3:9], ntree = 999, num.rep = 999,
                     importance = TRUE, na.action = na.omit)

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
## Show "importance" of variables: higher value mean more important:
mean(meta.lcbd.rf$rf$rsq)
sd(meta.lcbd.rf$rf$rsq)
meta.lcbd.rf$pval

mean(tp.lcbd.rf$rf$rsq)
sd(tp.lcbd.rf$rf$rsq)
tp.lcbd.rf$pval

mean(pa.lcbd.rf$rf$rsq)
sd(pa.lcbd.rf$rf$rsq)
pa.lcbd.rf$pval
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




#alpine dataset analysis

## read tibet plateau dataset
asv.table <- read.csv('../tp_data/asv.table.csv',sep=",", row.names=1)
asv.table <- as.matrix(asv.table)

### read in taxonomy
taxonomy <- read.csv('../tp_data/taxonomy.csv',sep=",",row.names=1)
taxonomy <- as.matrix(taxonomy)

### read in metadata
metadata <- read.csv("../tp_data/sample_data.csv", row.names=1, header = T)

### read in tree
total.tree <- read_tree('../tp_data/tree.nwk')

### read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "../tp_data/ref.seqs.fasta",
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

### import as phyloseq objects
asv.table <- otu_table(asv.table, taxa_are_rows = TRUE)
taxonomy <- tax_table(taxonomy)
metadata <- sample_data(metadata)

#### merge into one phyloseq object
tp_physeq <- phyloseq(asv.table, taxonomy, metadata, total.tree, ref_seqs)
tp_physeq

# beta diversity analysis
## local contribution to beta diversity (LCBD) analysis for total community
library(adespatial)
beta_tax_div <- beta.div(t(as.matrix(otu_table(tp_physeq))), 
                         method = "hellinger", sqrt.D = FALSE, samp = TRUE, 
                         nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)
env_div <- data.frame(LCBD = beta_tax_div$LCBD, sample_data(tp_phylo))
### calculate the average LCBD for each site
library(dplyr)
env_div_agg <-  env_div %>% 
  dplyr::select(-c(2:6, 9, 12, 17)) %>%
  group_by(Site, Sitegroup1) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE))
#write.csv(env_div_agg, file = '../tp_data/env_div_agg.csv')

### determine the relationships between LCBD and envs using linear regression modes
vars <- c("MAT", "MAP", "DOC", "S275_295", "SUVA254", "a320", "BIX", "HIX",
          "TN", "NH4_N", "DO", "pH", "Conductivity", "Salinity",
          "K", "Ca", "Na", "Mg")
mode <- lapply(vars, function(x) {
  lm(substitute(LCBD ~ i, list(i = as.name(x))), data = env_div_agg)})
sum.mode <- lapply(mode, broom::glance)
### normality test using Shapiro-Wilk test 
res <- lapply(mode, residuals)
norm_test <- lapply(res, shapiro.test)
norm_results <- data.frame(
  variables = vars, 
  w = sapply(norm_test, "[[", "statistic"), 
  pvalue = sapply(norm_test, "[[", "p.value")
)
norm_results

### extract the standardized regression coefficients
sd.coeff <- lapply(mode, QuantPsyc::lm.beta)
### arrange the table for plot
LCBD <- c(rep('LCBD', length(vars)))
sd.coeff <- sapply(sd.coeff, function(x){as.numeric(x[1])})
r.squared <- sapply(sum.mode, "[[", "r.squared")
adj.r.squared <- sapply(sum.mode, "[[", "adj.r.squared")
pvalue <- sapply(sum.mode, "[[", "p.value")
sig <- cut(pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
results <- data.frame(vars, LCBD, sd.coeff, r.squared, adj.r.squared, pvalue, sig)
results

#model selection
library(MASS)
library(glmulti)
A1 <- glmulti(LCBD ~ MAP + MAT + S275_295 + SUVA254 + a320 + BIX + HIX +
                TN + Conductivity + Salinity + Mg + K + Na, data = env_div_agg,
              level=1, fitfunction=lm, crit="aicc", confsetsize= 2^13, plotty = F, trace = 0)
top <- weightable(A1)
###  models with values more than 2 units away are considered substantially 
### less plausible than those with AICc values closer to that of the best model. 
### refrence:Anderson, D. R. (2007). Model based inference in the life sciences: A primer on evidence. New York: Springer. 
top_1 <- top[top$aicc <= min(top$aicc) + 2,] # 
top_1

modes_inf <- NULL
for(i in 1:nrow(top_1)){
  rse_sum <- summary(A1@objects[[i]])
  adj.r.squared <- rse_sum$adj.r.squared # obtain the adjust r squared
  multicollinearity <- any(car::vif(A1@objects[[i]]) > 2) # check the multicollinearity
  tmp <- data.frame(adj.r.squared, multicollinearity)
  if(is.null(modes_inf)){
    modes_inf<-tmp
  } else {
    modes_inf <- rbind(modes_inf,tmp)
  } 
}
modes_inf <- cbind(top_1, modes_inf)
modes_inf

vpa.mod <- varpart(env_div_agg$LCBD, ~ env_div_agg$HIX,
                   ~ env_div_agg$MAP)
plot(vpa.mod)

### heatmap using standardized regression coefficients to explore the relationship between the LCBD and environment factors
results_plot_data <- data.frame(group = c(rep('Total community', nrow(results))),results)

results_plot_data$vars <- factor(results_plot_data$vars,levels = rev(vars))
results_plot_data$group <- factor(results_plot_data$group)
p_env_div <- ggplot(aes(x=LCBD, y=vars, fill=sd.coeff), data=results_plot_data) +
  geom_tile() +
  scale_fill_gradient2(low='#1b9e77', mid='white', high='#d95f02') +
  geom_text(aes(label=sig), color="black", size=6) +
  labs(y=NULL, x=NULL, fill='Standardized regression coefficients') +
  facet_wrap( .~ group, ncol = 2) +
  theme_bw()+
  theme(legend.position="bottom", 
        panel.border = element_blank(),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
###  plot the linear regression relationships between the LCBD and best explained variables
#### total community
p_linear <- env_div_agg %>%
  dplyr::select(LCBD, MAP, HIX) %>%
  tidyr::gather(varibales, value, MAP:HIX, factor_key=TRUE) %>%
  ggplot(aes(value, LCBD)) +
  geom_point(size=3.5, alpha=0.8, aes(colour = as.factor(varibales))) +
  geom_smooth(method="lm", size=1, se=T, colour='black') +
  scale_color_manual(values = c('#1b9e77', '#d95f02')) +
  scale_y_continuous(limits = c(0, 0.01)) +
  facet_wrap( .~ varibales, scales="free_x", ncol = 2) +
  ylab('LCBD')+xlab('Values') +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        strip.text = element_text(size = 14),
        legend.position='none')

p_linear

## venn plot
library("VennDiagram")
venn.plot <- draw.pairwise.venn(area1 = 0.27, area2 = 0.31, cross.area = 0.21,
                                category=c('DOM proporties', 'Climate elements'),
                                fill = c( '#d95f02', '#1b9e77'), scaled = 0, 
                                ind = FALSE, cat.col = c(rep('black',2)), 
                                cat.cex = 1.2, cat.dist =  c(0.02, 0.02),
                                cat.pos = c(30, 330), margin = 0.05, lty = 'blank')
ggdraw() +
  draw_plot(p_env_div, x = 0, y = 0, width = 0.3, height = 1) +
  draw_plot(p_linear, x = 0.3, y = 0.5, width = 0.7, height = 0.5) +
  draw_plot(venn.plot, x = 0.3, y = 0, width = 0.7, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 14,
                  x = c(0, 0.3, 0.3), y = c(1, 1, 0.7))

## partial mantel test function for Tibatan Plateau and Pan Arctic
partial.mantel.fun <- function(phylo) {
  env.table <- data.frame(sample_data(phylo))
  #env.table <- env.table[complete.cases(env.table), ]
  otu_table <- as.matrix(t(otu_table(phylo)))
  otu_table_hel <- decostand(otu_table, 'hellinger')
  otu_table_hel_dist <- vegdist(otu_table_hel, 'bray',upper=F)
  df <- NULL
  vars <- c("DOC", "SUVA254", "a320", "MAP", "MAT", "pH")
  for (x in vars) {
    x.dist <- vegdist(scale(env.table[,x]), 'euclidean', na.rm = T)
    z.dist <- vegdist(scale(env.table[ , setdiff(vars, x)]), 'euclidean', na.rm = T)
    mode <- mantel.partial(x.dist, otu_table_hel_dist, z.dist, 
                           method = "pearson", permutations = 999, na.rm = T)
    r <- mode$statistic
    p <- mode$signif
    tmp <- data.frame(env = x, r = r, p.value = p)
    if(is.null(df))
      df <- tmp
    else
      df <- rbind(df ,tmp)
  }
  return(df)
}

all.total.commun.par.mant <- partial.mantel.fun(meta_physeq)
## PLOT
## devtools::install_github('hannet91/ggcor')
library(ggcor)
set.seed(123456)

all.par.man.tibble <- tibble(spec = c(rep('all.total.commun.par.mant', nrow(all.total.commun.par.mant))), 
                             all.total.commun.par.mant)
#par.man.tibble <- tibble(spec = c(rep('Total community composition', nrow(total.commun.par.mant))), 
#                         total.commun.par.mant)
vars <- c("MAT", "MAP", "DOC", "SUVA254", "a320", "pH")
env.table <- sample_data(meta_physeq)[ , vars]

mantel02 <- all.par.man.tibble %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.20, 0.4, Inf), 
                 labels = c("<0.20", "0.30-0.40", ">0.40"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"),
                       right = T))
quickcor(env.table[complete.cases(env.table),], type = "upper") + geom_square() + 
  add_link(mantel02, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE) +
  scale_color_manual(values = c('#d95f02', 'grey', '#1b9e77', '#3C5488FF')) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  geom_diag_label() + remove_axis("x")

## partial mantel test function for Tibatan Plateau
## partial mantel test function
partial.mantel.fun <- function(phylo) {
  env.table <- data.frame(sample_data(phylo))
  otu_table <- as.matrix(t(otu_table(phylo)))
  otu_table_hel <- decostand(otu_table, 'hellinger')
  otu_table_hel_dist <- vegdist(otu_table_hel, 'bray',upper=F)
  df <- NULL
  vars <- c("DOC", "TN", "NH4_N", "S275_295", "SUVA254",
            "a320", "MAP", "MAT", "DO", "pH", "Conductivity", 
            "Salinity", "Ca", "Mg", "K", "Na", "BIX", "HIX")
  for (x in vars) {
    x.dist <- vegdist(scale(env.table[,x]), 'euclidean')
    z.dist <- vegdist(scale(env.table[ , setdiff(vars, x)]), 'euclidean')
    mode <- mantel.partial(x.dist, otu_table_hel_dist, z.dist, 
                           method = "pearson", permutations = 999)
    r <- mode$statistic
    p <- mode$signif
    tmp <- data.frame(env = x, r = r, p.value = p)
    if(is.null(df))
      df <- tmp
    else
      df <- rbind(df ,tmp)
  }
  return(df)
}
tp.commun.par.mant <- partial.mantel.fun(tp_physeq)
## PLOT
## devtools::install_github('hannet91/ggcor')
library(ggcor)
set.seed(123456)

tp.par.man.tibble <- tibble(spec = c(rep('tp.commun.par.mant', nrow(tp.commun.par.mant))), 
                             tp.commun.par.mant)
#par.man.tibble <- tibble(spec = c(rep('Total community composition', nrow(total.commun.par.mant))), 
#                         total.commun.par.mant)
vars <- c("MAT", "MAP", "DOC", "S275_295", "SUVA254", "a320", "BIX", "HIX",
          "TN", "NH4_N", "DO", "pH", "Conductivity", "Salinity",
          "K", "Ca", "Na", "Mg")
env.table <- sample_data(tp_physeq)[ , vars]

mantel02 <- tp.par.man.tibble %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.2, 0.5, Inf), 
                 labels = c("<0.20", "0.20-0.5", ">0.50"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"),
                       right = T))
quickcor(env.table, type = "upper") + geom_square() + 
  add_link(mantel02, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE) +
  scale_color_manual(values = c('#d95f02', '#1b9e77', '#3C5488FF', 'grey')) +
  scale_size_manual(values = c(0.5, 2)) +
  geom_diag_label() + remove_axis("x")

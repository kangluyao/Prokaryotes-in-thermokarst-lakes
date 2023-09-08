# test the difference of the environmental factors between TP and PA
##standard error function
stderr <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
meta_env_sum <- metadata %>% dplyr::select(c('Region', 'DOC', 'SUVA254', 'a320',
                                             'MAP', 'MAT', 'pH')) %>% group_by(Region) %>%
  dplyr::summarise_all(list(mean = mean, sd = sd, se = stderr), na.rm = TRUE)
# write.csv(meta_env, file = 'E:/thermokast_lakes/water_microbes/meta_analysis/results/tables/meta_env.csv')
# test the difference
library(lme4)
library(lmerTest)
library(multcomp)
par(mfrow = c(2, 3))
hist(metadata$DOC)
hist(metadata$SUVA254)
hist(metadata$a320)
hist(metadata$MAP)
hist(metadata$MAT)
hist(metadata$pH)
par(mfrow = c(1, 1))

par(mfrow = c(2, 3))
hist(log(metadata$DOC))
hist(log(metadata$SUVA254))
hist(log(metadata$a320))
hist(log(metadata$MAP))
hist(log(metadata$MAT+15))
hist(log(metadata$pH))
par(mfrow = c(1, 1))

env_mode1 <- lmer(log(DOC) ~ Region + (1|Site), metadata)
summary(env_mode1)
env_mode2 <- lmer(log(SUVA254) ~ Region + (1|Site), metadata)
summary(env_mode2)
env_mode3 <- lmer(log(a320) ~ Region + (1|Site), metadata)
summary(env_mode3)
env_mode4 <- lmer(log(MAP) ~ Region + (1|Site), metadata)
summary(env_mode4)
env_mode5 <- lmer(MAT ~ Region + (1|Site), metadata)
summary(env_mode5)
env_mode6 <- lmer(log(pH) ~ Region + (1|Site), metadata)
summary(env_mode6)

# test the difference of the community structure between TP and PA
#PERMANOVA analysis and NMDS plot
meta.ord <- ordinate(meta_physeq, "NMDS", "bray")
meta.ord
NMDS_plot <- plot_ordination(meta_physeq, meta.ord, type = "samples", color = "Region") + 
  geom_point(size = 2.5) + scale_color_manual(values = c("#d95f02", "#1b9e77")) + 
  #stat_ellipse(type = 'norm', linetype = 1) + "#1b9e77")) + 
  #stat_ellipse(type = 'norm', linetype = 1) +
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.4, aes(fill = Region)) + 
  theme_bw() + 
  theme(legend.position = c(0.85, 0.88), 
        panel.grid = element_blank(), 
        axis.title = element_text(color = "black", size = 9), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = 8), 
        legend.title = element_blank(),
        legend.text = element_text(size = 8), 
        legend.key = element_blank(), 
        legend.key.size = unit(1, "line"),
        legend.background = element_rect(colour = "white"))

metadata <- as(sample_data(meta_physeq), "data.frame")
adonis2(phyloseq::distance(meta_physeq, method = "bray") ~ Region,
        data = metadata)

#diff taxa using LEfSe with microeco package
library(microeco)
library(file2meco)
meco_df <- file2meco::phyloseq2meco(meta_physeq)
#calculate the abundance table
m1 <- meco_df$cal_abund()
ps_genus <- phyloseq::tax_glom(meta_physeq, taxrank = 'Genus')
meco_genus_df <- phyloseq2meco(ps_genus)
#calculate the abundance table
m1_genus <- meco_genus_df$cal_abund()

m1_genus <- trans_diff$new(dataset = meco_genus_df, method = "lefse", 
                           group = "Region", alpha = 0.01, 
                           lefse_subgroup = NULL)
#t1$res_lefse is the LEfSe result
#t1$res_abund is the abundance information
lefse_plot <- m1_genus$plot_diff_bar(threshold = 4, 
                                     color_values = c('#d95f02', '#1b9e77'), 
                                     axis_text_y = 9,
                        plot_vertical = TRUE, width = 0.5) +
  theme(legend.position = c(0.85, 0.15), 
        panel.grid = element_blank(), 
        axis.title = element_text(color = "black", size = 9), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = 8), 
        legend.title = element_blank(),
        legend.text = element_text(size = 8), 
        legend.key = element_blank(), 
        legend.key.size = unit(1, "line"),
        legend.background = element_rect(colour = "white"))

# Combining Plots
library(cowplot)
NMDS_lefse_plots <- ggdraw() +
  draw_plot(NMDS_plot, x = 0, y = 0, width = 0.48, height = 1) +
  draw_plot(lefse_plot, x = 0.52, y = 0, width = 0.48, height = 1) +
  draw_plot_label(label = c("(a)", "(b)"), size = 9,
                  x = c(0, 0.5), y = c(1, 1))

NMDS_lefse_plots 
# ggsave(NMDS_lefse_plots,
#        file = "E:/thermokast_lakes/water_microbes/meta_analysis/results/figs/2022.08.23/NMDS_lefse.pdf",
#        width = 8.9, height = 4.5, units = 'in', device='pdf', dpi=300)
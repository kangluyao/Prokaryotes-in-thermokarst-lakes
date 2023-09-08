# function for checking the outlier
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var), eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
  boxplot(var_name, main = "With outliers")
  hist(var_name, main = "With outliers", xlab = NA, ylab = NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main = "Without outliers")
  hist(var_name, main = "Without outliers", xlab = NA, ylab = NA)
  title("Outlier Check", outer = TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1)/sum(!is.na(var_name)) * 100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  response <- readline(prompt = "Do you want to remove outliers and to replace with NA? [yes/no]: ")
  if (response == "y" | response == "yes") {
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "n")
    return(invisible(dt))
  } else {
    cat("Nothing changed", "n")
    return(invisible(var_name))
  }
}

# function for arranging the OTU table at the order and genus level
arrange.tab <- function(phylo, N, taxrank, vect) {
  phylo.rel <- transform_sample_counts(phylo, function(x) x / sum(x))
  subphylo.rel <- tax_glom(phylo.rel, taxrank, NArm = TRUE) #here, the unclassified taxa are removed with the option NArm = TRUE, cause the Unclassified taxa at the genus level are not shown in the Fig S3. 
  ra.tab <- otu_table(subphylo.rel)
  MRA <- rowMeans(ra.tab)
  group <- tax_table(subphylo.rel)[,vect]
  mra.tab <- data.frame(group,MRA)
  colnames(mra.tab) <- c('level1', 'level2', 'MRA')
  #arrange the class table
  mra.tab_level1 = mra.tab %>% group_by(level1) %>% 
    summarise(sum_MRA = sum(MRA)) %>% 
    arrange(desc(sum_MRA))
  top_N_level1 = mra.tab_level1[1:N, ]$'level1'
  top_N_tab = mra.tab[mra.tab$'level1' %in% top_N_level1, ]
  mra.tab_level2 = top_N_tab %>% group_by(level2) %>% 
    summarise(sum_MRA = sum(MRA)) %>% 
    arrange(desc(sum_MRA))
  order_level2 = mra.tab_level2$'level2'
  top_N_tab$'level1' = factor(top_N_tab$'level1', ordered = T, levels = top_N_level1)
  top_N_tab$'level2' = factor(top_N_tab$'level2', ordered = T, levels = rev(order_level2))
  top_N_tab$MRA = top_N_tab$MRA*100
  return(top_N_tab)
}
# env table
DOM_env <-  metadata %>% dplyr::select(c('DOC', 'SUVA254', 'a320', 'MAP', 'MAT', 'pH'))
DOM_env[1:5, 1:3]
# extract the most abundant 10% of taxa in at least half of the samples at the order level
phylo.rel <- transform_sample_counts(meta_physeq, function(x) x / sum(x))
subphylo.rel <- tax_glom(phylo.rel, taxrank = 'Order', NArm = F)
f1<- filterfun_sample(topf(0.9))
wh1 <- genefilter_sample(subphylo.rel, f1, A=(1/2*nsamples(subphylo.rel)))
sum(wh1)
meta.com.ord <- prune_taxa(wh1, subphylo.rel)
domin_tax_names <- as.character(tax_table(meta.com.ord)[, 4])

# The composition of the domiant taxa at the genus level
top10_order_meta <- arrange.tab(meta_physeq, 10, 'Genus', c(4,6))
top_domin_order_meta <- subset(top10_order_meta, level1 %in% domin_tax_names)
mra.tab_level1 = top_domin_order_meta %>% group_by(level1) %>% 
  summarise(sum_MRA = sum(MRA)) %>% 
  arrange(desc(sum_MRA))
mra.tab_level1
mra.tab_level2 = top_domin_order_meta %>% group_by(level2) %>% 
  summarise(sum_MRA = sum(MRA)) %>% 
  arrange(desc(sum_MRA))
mra.tab_level2 [1:20, ]
order_level2 = mra.tab_level2$'level2'
#stack bar plot
taxa_domin_genus_barplot <- ggplot(top_domin_order_meta, aes(fill=level2, y=MRA, x=level1)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual('Genus', breaks = order_level2[1:25], 
                    values = rep(c(mycols[1:12], mycols[-c(1:12)]), 20)[1:nrow(top_domin_order_meta)]) + 
  labs(x = 'Orders', y = 'Mean relative abundance (%)') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +
  theme_classic()+
  theme(legend.position = c(0.6, 0.7),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=9),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=8),
        axis.text.x = element_text(colour='black', size = 8, angle = 45, hjust = 1),
        legend.text=element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "white"))
taxa_domin_genus_barplot

# check the outliers of data and prepare the data for plot
order_table <- otu_table(meta.com.ord)
order_table[1:7, 1:5]
test.dat <- data.frame(Abund = colSums(order_table), DOM_env)

# Frequence distribution of the mean relative abundance of dominant taxa
plot_df <- data.frame(x1 = rep("abundance", ncol(order_table)), Abund = colSums(order_table) *100)
source("https://raw.githubusercontent.com/samhforbes/PupillometryR/master/R/geom_flat_violin.R")
rel_abun_domin_plot <- ggplot(plot_df, aes(x1, Abund)) + 
  geom_flat_violin(aes(fill = "tomato3"), position = position_nudge(x = .25), color = "black") + 
  geom_jitter(aes(color = "tomato3"), width = 0.1) + 
  geom_boxplot(width = .1, position = position_nudge(x = 0.25), fill = "white", size = 0.5) + 
  labs(x = NULL, y = "Mean relative abundance of dominant taxa (%)") +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=9),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text=element_text(size=8),
        legend.position = "none")
rel_abun_domin_plot

# Combining Plots
library(cowplot)
domian_genera_plot <- ggdraw() +
  draw_plot(rel_abun_domin_plot, x = 0, y = 0.18, width = 0.32, height = 0.82) +
  draw_plot(taxa_domin_genus_barplot, x = 0.38, y = 0, width = 0.62, height = 1) +
  draw_plot_label(label = c("(a)", "(b)"), size = 9,
                  x = c(0, 0.35), y = c(1, 1))
domian_genera_plot

# ggsave(domian_genera_plot,
#        file = "E:/thermokast_lakes/water_microbes/meta_analysis/results/revision/domian_genera_plot11.pdf",
#        width = 7.6, height = 4.5, units = 'in', device='pdf', dpi=300)


# Test the correlation between the relative abundance of dominant taxa and environmental variables
outlierKD(test.dat, Abund)
y
outlierKD(test.dat, pH)
y
outlierKD(test.dat, MAP)
y
outlierKD(test.dat, MAT)
y
outlierKD(test.dat, SUVA254)
y
outlierKD(test.dat, DOC)
y
outlierKD(test.dat, a320)
y

# set the environmmental vatiables
vars<-c('DOC', 'SUVA254', 'a320', 'MAP', 'MAT', 'pH')
# spearman test and plot
main_theme <- theme_bw() + 
  theme(axis.ticks.length = unit(0.4, "lines"), 
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour = "black", size = 12), 
        axis.title.y = element_text(colour = "black", size = 12), 
        axis.text.y = element_text(colour = "black", size = 10), 
        axis.text.x = element_text(colour = "black", size = 10), 
        panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid = element_blank())

# spearman
test.dat %>% pivot_longer(cols = -c(Abund), names_to = "env_name", values_to = 'value') %>%
  mutate(env_name = factor(env_name, levels = c('MAT', 'MAP', 'DOC', 'SUVA254', 'a320', 'pH'))) %>%
  mutate(Abund = Abund * 100) %>%
  ggplot(aes(x = value, y = Abund)) +
  geom_point(shape = 19, size = 1, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(value, Abund, label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                   cor.coef.name = "rho", p.accuracy = 0.001, r.accuracy = 0.01,
                   method = "spearman", label.x.npc = 0.3, label.y.npc = 0.05, size = 4) +
  xlab('Environmental variables') +
  ylab('Relative abundance of abundant taxa (%)') +
  facet_wrap( ~ env_name, scales = 'free_x') +
  main_theme

# test the relationship between the relative abundance of Burkholderiales and DOM
tax_env_table <- data.frame(tax_table(subphylo.rel)[, 4], otu_table(subphylo.rel)) %>% 
    mutate(MRA = rowMeans(dplyr::select(., rownames(sample_data(subphylo.rel))))) %>%
    arrange(desc(MRA)) %>% dplyr::top_n(5, MRA) %>%
    dplyr::select(., -c('MRA')) %>% 
    bind_rows(summarise_all(., ~if(is.numeric(.)) 1-sum(.) else "Others")) %>%
    remove_rownames %>% column_to_rownames(var = "Order") %>%
    t() %>% data.frame(., DOM_env)
# cheak the outliers
outlierKD(tax_env_table, Burkholderiales)
y
outlierKD(tax_env_table, SUVA254)
y
outlierKD(tax_env_table, DOC)
y
outlierKD(tax_env_table, a320)
y
# relationship between the relative abundance of Burkholderiales and DOM properties
tax_env_table %>% dplyr::select('Burkholderiales', 'DOC', 'SUVA254', 'a320') %>%
  pivot_longer(cols = -c(Burkholderiales), names_to = "env_name", values_to = 'value') %>%
  mutate(env_name = factor(env_name, levels = c('DOC', 'SUVA254', 'a320'))) %>%
  mutate(Burkholderiales = Burkholderiales * 100) %>%
  ggplot(aes(x = value, y = Burkholderiales)) +
  geom_point(shape = 19, size = 1, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(value, Burkholderiales, label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                   cor.coef.name = "rho", p.accuracy = 0.001, r.accuracy = 0.01,
                   method = "spearman", label.x.npc = 0.3, label.y.npc = 0.05, size = 4) +
  xlab('DOM properties') +
  ylab('Relative abundance of Burkholderiales (%)') +
  facet_wrap( ~ env_name, scales = 'free_x') +
  main_theme
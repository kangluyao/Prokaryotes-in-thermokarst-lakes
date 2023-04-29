#diff class using microeco package
library(microeco)
library(file2meco)
meco_df <- phyloseq2meco(meta_physeq)

# Functional analysis
# create object of trans_func
m2 <- trans_func$new(meco_df)
# mapping the taxonomy to the database
# the function can recognize prokaryotes or fungi automatically.
m2$cal_spe_func()
# return m2$res_spe_func, 1 represent function exists, 0 represent no or cannot confirmed.
m2$res_spe_func[1:5, 1:2]
m2$func_group_list
#carbon cycling taxa analysis
# If you want to change the group list, reset the list t2$func_group_list
C_cycle_process <- c("chitinolysis", "cellulolysis", "xylanolysis", "ligninolysis", 
                     'aromatic_hydrocarbon_degradation',
                     'aromatic_compound_degradation', 'hydrocarbon_degradation')
m2$func_group_list$`C-cycle` <- C_cycle_process
# use show_prok_func to see the detailed information of prokaryotic traits
m2$show_prok_func("methanotrophy")
# calculate the percentages for communities
m2$cal_spe_func_perc(abundance_weighted = FALSE)
m2$res_spe_func_perc[1:5, 1:2]

meta_c_cycle <- data.frame(m2$res_spe_func_perc[ ,colnames(m2$res_spe_func_perc) %in% c(m2$func_group_list$`C-cycle`)], 
                           m2$sample_table[, c('Site', 'Region')])
vars <- colnames(m2$res_spe_func_perc[ ,colnames(m2$res_spe_func_perc) %in% c(m2$func_group_list$`C-cycle`)])

model_df <- meta_c_cycle %>% dplyr::select(vars) %>%
  gather("C_process", "propertion")

multi_aov <- function(Y, X, paired = FALSE, p.adj = "none") { #X and Y are vectors
  aa <- levels(as.factor(X))
  an <- as.character(c(1:length(aa)))
  tt1 <- matrix(nrow = length(aa), ncol = 6)    
  for (i in 1:length(aa))
  {
    temp <- Y[X == aa[i]]
    tt1[i, 1] <- mean(temp, na.rm = TRUE)
    tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
    tt1[i, 3] <- sd(temp, na.rm = TRUE)
    tt1[i, 4] <- min(temp, na.rm = TRUE)
    tt1[i, 5] <- max(temp, na.rm = TRUE)
    tt1[i, 6] <- length(temp)
  }
  tt1 <- data.frame(aa, tt1)
  colnames(tt1) <- c("group", "mean", "se", "sd", "min", "max", "n")
  require(agricolae)
  Xn <- factor(X, labels = an)
  sig <- "ns"
  model <- aov(Y ~ Xn)    
  if (paired == TRUE & length(aa) == 2)
  {
    coms <- t.test(Y ~ Xn, paired = TRUE)
    pp <- coms$p.value
  }    else
  {
    pp <- anova(model)$Pr[1]
  }    
  if (pp <= 0.1)
    sig <- "."
  if (pp <= 0.05)
    sig <- "*"
  if (pp <= 0.01)
    sig <- "**"
  if (pp <= 0.001)
    sig <- "***"
  if(pp <= 0.05) {
    comp <- LSD.test(model,
                     "Xn",
                     alpha = 0.05,
                     p.adj = "none",
                     group = TRUE)
    gror <- comp$groups[order(rownames(comp$groups)), ]
    tt1$sig <- gror$groups
  }
  list(comparison = tt1, p.value = pp)
}
plot_df <- multi_aov(model_df$propertion, model_df$C_process, paired = FALSE, p.adj = "none")

# plot
ggplot(plot_df$comparison,aes(x = group, y = mean)) +
  geom_bar(position = 'dodge', stat = 'identity', 
           fill= rev(c("#da5724", '#3b86bc', '#915c83', '#eb8f70', '#c6dfed', '#f4c2c2')), 
           colour = 'black', width = 0.7) +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width=.2,  position = position_dodge(0.7)) +
  scale_x_discrete(limits = c('cellulolysis', 'xylanolysis', 'chitinolysis',   
                              'aromatic_hydrocarbon_degradation',
                              'aromatic_compound_degradation', 'hydrocarbon_degradation')) +
  labs(x = 'Carbon degradation', y = 'Mean proportion (%)') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  geom_text(aes(label = sig, y = mean + se + 0.02*max(mean)),
            position = position_dodge(0.7), vjust = 0)+
  theme_bw()+
  theme(legend.position = NULL,
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1)
  )

# test the significance of the functional taxa among the region with LMM 
library(lme4)
library(lmerTest)
mode <- lapply(vars, function(x) {
  lmer(substitute(i ~ Region + (1|Site), list(i = as.name(x))), data = meta_c_cycle)})
summary.model <- function(model){
  F.value <- anova(model)$'F value'
  p.value <- anova(model)$'Pr(>F)'
  p.stars <- function(p.values) {
    unclass(symnum(p.values, corr = FALSE, 
                   na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")))}
  sig <- p.stars(p.value)
  results<-data.frame(F.value, p.value, sig)
  return(results)
}
df <- NULL
for(i in 1:length(vars)) {
  tmp <- summary.model(mode[[i]])
  if (is.null(df)){
    df <- tmp
  } else {
    df <- rbind(df, tmp)
  }
}
result_carbon <-data.frame(vars, df)
result_carbon

### arrange data for plot
library(tidyverse)
library(reshape)
melt_df <- melt(meta_c_cycle[, !colnames(meta_c_cycle) %in% 'Site'], id.vars = c('Region'))
#piechart for all samples
all_df <- melt_df %>% 
  dplyr::select(c(variable, value)) %>%
  group_by(variable) %>% 
  dplyr::summarize(value = sum(value)) %>%
  arrange(desc(value))  %>%
  mutate(prop = value / sum(melt_df$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
  mutate(variable=factor(variable, levels=variable))


library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
# Create a ggplot with 11 colors 
#barplot
new_df <- ddply(melt_df, c('Region','variable'), summarise,
                mean = mean(value), sd = sd(value),
                sem = sd(value)/sqrt(length(value)))
sig <- c('b', 'b', 'a', 'a', 'b', 'a', 
         'a', 'a', 'a', 'a', 'a', 'b')
new_df <- cbind(new_df ,sig)
new_df$Region <- factor(new_df$Region, levels = c('Tibetan Plateau', 'Pan-Arctic'))
new_df
#plot
carbon_fun_plot <- ggplot(new_df,aes(x = variable, y = mean, fill = Region)) +
  geom_bar(position = 'dodge', stat = 'identity', colour = 'black', width = 0.7) +
  scale_fill_manual(values=c('#1b9e77', '#d95f02')) +
  geom_errorbar(aes(ymin = mean, ymax = mean + sem), width=.2,  position = position_dodge(0.7)) +
  scale_x_discrete(limits = c('cellulolysis', 'xylanolysis', 'chitinolysis',   
                              'aromatic_hydrocarbon_degradation',
                              'aromatic_compound_degradation', 'hydrocarbon_degradation')) +
  labs(x = 'Carbon cycle', y = 'Mean relative abundance (%)', fill = "Region") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  geom_text(aes(label = sig, y = mean + sem + 0.02*max(mean)), 
            position = position_dodge(0.7), vjust = 0)+
  theme_bw()+
  theme(legend.position = c(0.7,0.88),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
carbon_fun_plot


#extract otus involved in carbon cycling
meta_c_cycle_otus <- m2$res_spe_func[ ,c(m2$func_group_list$`C-cycle`)]

meta_carbon_phy <- subset_taxa(meta_physeq, OTU %in% meta_c_cycle_otus)
meta_carbon_phy_rel <- subset_taxa(meta_physeq_rel, OTU %in% meta_c_cycle_otus)

meta_meco_df <- phyloseq2meco(meta_carbon_phy)
m3 <- meta_meco_df$cal_abund()
#trans_diff class, The third approach is rf, which depends on the random forest[14, 15] and the non-parametric test. 
# use Genus level for parameter rf_taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
m3 <- trans_diff$new(dataset = meta_meco_df, method = "rf", 
                     group = "Region", rf_taxa_level = "Genus")

# m3$res_rf is the result stored in the object
# plot the result
res_abund <- data.frame(genus = sapply(strsplit(m3$res_abund$Taxa, "\\|"), "[[", 6),
                        m3$res_abund[, c(2:6)])
res_rf <- data.frame(genus = sapply(strsplit(m3$res_rf$Taxa, "\\|"), "[[", 6),
                     m3$res_rf[, c(2:3)])
top_genus <- res_rf$genus[1:10]

res_abund_plot <- res_abund %>% filter(genus %in% top_genus) %>% 
  mutate(Group = factor(Group, levels = c('Tibetan Plateau', 'Pan-Arctic'))) %>%
  ggplot(aes(x = genus, y = 100*Mean, fill = Group))+
  geom_bar(position = 'dodge', stat = 'identity', colour = 'black', width = 0.7)+
  scale_fill_manual(values=c('#1b9e77', '#d95f02'))+
  geom_errorbar(aes(ymin = 100*Mean, ymax = 100*Mean + 100*SE), width=.2,  position = position_dodge(0.7))+
  scale_x_discrete(limits = top_genus) +
  labs(x = 'Genus', y = 'Relative abundance (%)', fill = "Group") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +
  theme_bw()+
  theme(legend.position = c(0.15,0.9),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
res_abund_plot
library(grid)
vie <- viewport(width=0.5, height=0.45, x=0.68, y=0.75)
res_rf_plot <- res_rf %>% filter(genus %in% top_genus) %>% 
  mutate(genus = factor(genus, levels = top_genus)) %>%
  ggplot(aes(x = genus, y = MeanDecreaseGini)) +
  geom_bar(position="stack", stat="identity", alpha=0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20))+
  labs(x = NULL, y = "MeanDecreaseGini")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size=12, colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank())
print(res_rf_plot, vp = vie)
ggsave(file = "./meta_analysis/results/figs/res_rf_plot.pdf",
       width = 12, height = 10, units = 'in', device='pdf', dpi=300)
carbon_plot <- ggdraw() +
  draw_plot(carbon_fun_plot, x = 0, y = 0.45, width = 1, height = 0.55) +
  draw_plot(carbon_rf_plot, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("(a)", "(b)"), size = 14,
                  x = c(0, 0), y = c(1, 0.45))

carbon_plot
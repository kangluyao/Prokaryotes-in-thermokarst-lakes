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
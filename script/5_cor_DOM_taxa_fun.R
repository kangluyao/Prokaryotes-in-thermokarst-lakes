##prune out Order below 1% in each sample and prevalence lower than 50/100 at Order level
meta_physeq_rel <- microbiome::transform(meta_physeq, "compositional")
meta.com.ord <- microbiome::aggregate_rare(meta_physeq_rel, level = "Order", 
                                           detection = 1/100, prevalence = 50/100)


dat_pr_high = filter_taxa(meta_physeq, function(x) {
  (sum(x > 1) > 306 * 0.1) & (mean(x/sample_sums(meta_physeq)) > 0.001)
}, prune = T)
order_table <- otu_table(meta.com.ord)
order_table[1:9, 1:5]

#plot
main_theme <- theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                    panel.grid = element_blank())+ 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour='black', size=14),
        axis.title.y = element_text(colour='black', size=14),
        axis.text.y = element_text(colour='black', size=12),
        axis.text.x = element_text(colour = "black", size = 12))

test.dat <- data.frame(Other = as.vector(order_table[c('Other'), ]), DOM_env)
outlierKD(test.dat, Other)
y
outlierKD(test.dat, DOC)
y
p_doc <- ggplot(test.dat, aes(x = DOC, y = log(Other + 1))) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(DOC, log(Other + 1), label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254')+ylab('Burkholderiales') +
  # scale_y_continuous(expand = c(0,0), limits = c(0.2, 1)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  main_theme
p_doc

test.dat <- data.frame(Other = as.vector(order_table[c('Other'), ]), DOM_env)
outlierKD(test.dat, Other)
y
outlierKD(test.dat, SUVA254)
y
p_suv254 <- ggplot(test.dat, aes(x = log(SUVA254), y = log(Other + 1))) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(log(SUVA254), log(Other + 1), label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254') + ylab('Burkholderiales') +
  scale_y_continuous(expand = c(0,0), limits = c(-0.25, 0.75)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  main_theme
p_suv254

test.dat <- data.frame(Other = as.vector(order_table[c('Other'), ]), DOM_env)
outlierKD(test.dat, Other)
y
outlierKD(test.dat, a320)
y
p_a320 <- ggplot(test.dat, aes(x = log(a320), y = log(Other + 1))) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(log(a320), log(Other + 1), label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254') + ylab('Burkholderiales') +
  scale_y_continuous(expand = c(0,0), limits = c(-0.25, 0.75)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  main_theme
p_a320

test.dat <- data.frame(Other = as.vector(order_table[c('Other'), ]), DOM_env)
outlierKD(test.dat, Other)
y
outlierKD(test.dat, pH)
y
p_pH <- ggplot(test.dat, aes(x = pH, y = log(Other + 1))) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(pH, log(Other + 1), label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254')+ylab('Burkholderiales') +
  scale_y_continuous(expand = c(0,0), limits = c(-0.25, 0.75)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  main_theme
p_pH

test.dat <- data.frame(Other = as.vector(order_table[c('Other'), ]), DOM_env)
outlierKD(test.dat, Other)
y
outlierKD(test.dat, MAP)
y
p_MAP <- ggplot(test.dat, aes(x = MAP, y = log(Other+1))) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(MAP, log(Other+1), label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254')+ylab('Burkholderiales') +
  scale_y_continuous(expand = c(0,0), limits = c(-0.5, 1)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  main_theme
p_MAP

test.dat <- data.frame(Other = as.vector(order_table[c('Other'), ]), DOM_env)
outlierKD(test.dat, Other)
y
outlierKD(test.dat, MAT)
y
p_MAT <- ggplot(test.dat, aes(x = MAT, y = log(Other + 1))) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(MAT, log(Other+1), label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254')+ylab('Burkholderiales') +
  scale_y_continuous(expand = c(0,0), limits = c(-0.5, 1)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  main_theme
p_MAT


allmean <- rowMeans(order_table)
regiontype <- as.factor(sample_data(meta.com.ord)$Region)
table(regiontype)
mean_region <- data.frame(Pan_Arctic = rowMeans(otu_table(subset_samples(meta.com.ord, Region == 'Pan-Arctic'))),
                          Tibetan_Plateau = rowMeans(otu_table(subset_samples(meta.com.ord, Region == 'Tibetan Plateau'))))
## arrange the data
rel_abun_dat_ord <- data.frame(Order = rownames(mean_region), mean_region, All_mean = allmean)
rel_abun_dat_ord <- dplyr::arrange(rel_abun_dat_ord, desc(All_mean))
Order_levels <- rel_abun_dat_ord$Order[!rel_abun_dat_ord$Order %in% c('Other', "uncultured")]

tax_table <- as.data.frame(t(data.frame(order_table))) %>% dplyr::select(Order_levels)
ncol(tax_table)
DOM_env <-  metadata %>% dplyr::select(c('DOC', 'SUVA254', 'a320', 'MAP', 'MAT', 'pH'))
DOM_env[1:5, 1:3]

#Now calculate the correlation between individual Taxa and the environmental data
x <- log(tax_table)
y <- scale(DOM_env)
groups <- c(rep('All', 306))
#You can use kendall, spearman, or pearson below:
method <- "spearman"
df <- NULL
for (i in colnames(x)) {
  for (j in colnames(y)) {
    for (k in unique(groups)) {
      a <- x[groups == k, i, drop = F]
      b <- y[groups == k, j, drop = F]
      tmp <- c(i, j, cor(a[complete.cases(b), ], b[complete.cases(b), ], use = "everything", method = method), cor.test(a[complete.cases(b),
      ], b[complete.cases(b), ], method = method)$p.value, k)
      if (is.null(df)) {
        df <- tmp
      } else {
        df <- rbind(df, tmp)
      }
    }
  }
}

df <- data.frame(row.names = NULL, df)
colnames(df) <- c("Index", "Env", "Correlation", "Pvalue", "Type")
df$Pvalue <- as.numeric(as.character(df$Pvalue))
df$AdjPvalue <- rep(0, dim(df)[1])
df$Correlation <- as.numeric(as.character(df$Correlation))
#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
# 1 -> donot adjust
# 2 -> adjust Env + Type (column on the correlation plot)
# 3 -> adjust Taxa + Type (row on the correlation plot for each type)
# 4 -> adjust Taxa (row on the correlation plot)
# 5 -> adjust Env (panel on the correlation plot)
adjustment_label <- c("NoAdj", "AdjEnvAndType", "AdjTaxaAndType", "AdjTaxa", "AdjEnv")
adjustment <- 5

if (adjustment == 1) {
  df$AdjPvalue <- df$Pvalue
} else if (adjustment == 2) {
  for (i in unique(df$Env)) {
    for (j in unique(df$Type)) {
      sel <- df$Env == i & df$Type == j
      df$AdjPvalue[sel] <- p.adjust(df$Pvalue[sel], method = "BH")
    }
  }
} else if (adjustment == 3) {
  for (i in unique(df$Taxa)) {
    for (j in unique(df$Type)) {
      sel <- df$Taxa == i & df$Type == j
      df$AdjPvalue[sel] <- p.adjust(df$Pvalue[sel], method = "BH")
    }
  }
} else if (adjustment == 4) {
  for (i in unique(df$Taxa)) {
    sel <- df$Taxa == i
    df$AdjPvalue[sel] <- p.adjust(df$Pvalue[sel], method = "BH")
  }
} else if (adjustment == 5) {
  for (i in unique(df$Env)) {
    sel <- df$Env == i
    df$AdjPvalue[sel] <- p.adjust(df$Pvalue[sel], method = "BH")
  }
}
#Now we generate the labels for signifant values
df$Significance <- cut(df$Pvalue, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ""))

#ignore NAs
df <- df[complete.cases(df),]

#We want to reorganize the Env data based on they appear
df$Env <- factor(df$Env, levels = rev(c('DOC', 'SUVA254', 'a320', 'MAP', 'MAT', 'pH')))
df$Index <- factor(df$Index, levels = c(Order_levels))
#df$Taxa<-factor(df$Taxa,levels=rev(taxa_list))
#df$Correlation[which(abs(df$AdjPvalue) >= 0.05)]<-0
#We use the function to change the labels for facet_grid in ggplot2
#

#correlation plot
p_env_taxa_cor <- ggplot(aes(x = Env, y = Index, col = Correlation), data = df) +
  geom_point(aes(size = abs(Correlation)))+ 
  scale_color_gradient2(low = '#1b9e77', mid = 'white', high = '#d95f02') +
  theme_bw(base_size = 12)+
  geom_text(aes(label = Significance), color = "black", size = 5)+
  scale_size(range = c(1, 8),breaks = c(0.1, 0.2, 0.4, 0.6, 0.8), limits = c(0, 0.7))+ 
  #theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = 'Environmental factors', y = '', col = "Spearman's r")+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color = 'black', size = 14),
        axis.ticks.length = unit(0.25, "lines"), axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(colour = 'black', size = 12, angle = 45, hjust = 1),
        legend.position = 'right',
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 9),
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white")) +
  coord_flip()

print(p_env_taxa_cor)

tax_env_table <-  data.frame(tax_table, DOM_env)
# check the outlier
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

tax_env_table <- data.frame(tax_table, DOM_env)
outlierKD(tax_env_table, Burkholderiales)
y
outlierKD(tax_env_table, Burkholderiales)
y
outlierKD(tax_env_table, Burkholderiales)
y
outlierKD(tax_env_table, SUVA254)
y
outlierKD(tax_env_table, DOC)
y
outlierKD(tax_env_table, a320)
y
# relationship between the relative abundance of Burkholderiales and DOM properties
nrow(tax_env_table)
p_doc <- ggplot(tax_env_table, aes(x = DOC, y = Burkholderiales)) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(DOC, Burkholderiales, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254')+ylab('Burkholderiales') +
  # scale_y_continuous(expand = c(0,0), limits = c(0.2, 1)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank())+ 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour='black', size=14),
        axis.title.y = element_text(colour='black', size=14),
        axis.text.y = element_text(colour='black', size=12),
        axis.text.x = element_text(colour = "black", size = 12))

p_suv254 <- ggplot(tax_env_table, aes(x = SUVA254, y = Burkholderiales)) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(SUVA254, Burkholderiales, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254')+ylab('Burkholderiales') +
  # scale_y_continuous(expand = c(0,0), limits = c(0.2, 1)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank())+ 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour='black', size=14),
        axis.title.y = element_text(colour='black', size=14),
        axis.text.y = element_text(colour='black', size=12),
        axis.text.x = element_text(colour = "black", size = 12))

p_a320 <- ggplot(tax_env_table, aes(x = a320, y = Burkholderiales)) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(a320, Burkholderiales, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254')+ylab('Burkholderiales') +
  # scale_y_continuous(expand = c(0,0), limits = c(0.2, 1)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank())+ 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour='black', size=14),
        axis.title.y = element_text(colour='black', size=14),
        axis.text.y = element_text(colour='black', size=12),
        axis.text.x = element_text(colour = "black", size = 12))

# Combining Plots
env_tax_plot <- ggdraw() +
  draw_plot(p_env_taxa_cor, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(p_doc, x = 0, y = 0, width = 1/3, height = 0.5) +
  draw_plot(p_suv254, x = 1/3, y = 0, width = 1/3, height = 0.5) +
  draw_plot(p_a320, x = 2/3, y = 0, width = 1/3, height = 0.5) +
  draw_plot_label(label = c("A", "B", 'C', 'D'), size = 14,
                  x = c(0, 0, 1/3, 2/3), y = c(1, 0.5, 0.5, 0.5))
env_tax_plot


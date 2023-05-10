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
# env table
DOM_env <-  metadata %>% dplyr::select(c('DOC', 'SUVA254', 'a320', 'MAP', 'MAT', 'pH'))
DOM_env[1:5, 1:3]
# extract the most abundant 10% of taxa in at least half of the samples at the order level
# subphylo <- tax_glom(meta_physeq, taxrank = 'Order')
# subphylo.rel <- transform_sample_counts(subphylo, function(x) x / sum(x) )
f1<- filterfun_sample(topf(0.9))
wh1 <- genefilter_sample(subphylo.rel, f1, A=(1/2*nsamples(subphylo.rel)))
sum(wh1)
meta.com.ord <- prune_taxa(wh1, subphylo.rel)

# check the outliers of data and prepare the data for plot
order_table <- otu_table(meta.com.ord)
order_table[1:7, 1:5]
test.dat <- data.frame(Abund = colSums(order_table), DOM_env)
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

# linear regression test
summary.is <- function(vars, df) {
  require(dplyr)
  model <- lapply(vars, function(x) {
    lm(substitute(Abund ~ i, list(i = as.name(x))), data = df)
  })
  r.squre <- round(as.vector(unlist(lapply(model, function(x) summary(x)$r.squared))), 3)
  p.value <- round(as.vector(unlist(lapply(model, function(x) anova(x)$"Pr(>F)"[1]))), 3)
  p.stars = function(p.values) {
    unclass(symnum(p.values, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " ")))
  }
  sig <- as.vector(unlist(lapply(p.value, p.stars)))
  normality <- round(as.vector(unlist(lapply(model, function(x) shapiro.test(residuals(x))$p.value))), 3)
  summary.table <- data.frame(vars, r.squre, p.value, sig, normality)
  return(summary.table)
}
summary.is(vars, test.dat)

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
# other regression models
source("https://raw.githubusercontent.com/kangluyao/source_R/main/multivariables_regression_plot.R")
# linear model
vars = c('MAP', 'MAT', 'pH', 'DOC', 'SUVA254', 'a320')
plot_cor_df(
  data = test.dat,
  y.name = 'Abund',
  x.names = vars,
  method = "spearman",
  ncol = 3)

# quadratic model
plot_poly_df(
  data = test.dat,
  y.name = 'Abund',
  x.names = vars,
  poly = 2,
  ncol = 3)

# generalized additive model
plot_gam_df(
  data = test.dat,
  y.name = 'Abund',
  x.names = vars,
  k = 2,
  ncol = 3)

plot_loess_df(
  data = test.dat,
  y.name = 'Abund',
  x.names = vars,
  ncol = 3)

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

# determine the composition of the dominant taxa
tax_table(meta.com.ord) 
  
  
  
  
  
  
  



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


#alpha diversity
tp_diversity <- estimate_richness(tp_physeq, measures = c("Chao1", 'Shannon', 'Simpson'))
tp_diversity <- tp_diversity[ ,c("Chao1", 'Shannon', 'Simpson')]

tp_diversity <- cbind(diversity, Region = meta_table$Region, 
                        Site = meta_table$Site)

#phylogenetic diversity
library(picante)
df.pd <- pd(t(otu.table), total.tree, include.root=T)

#functional diversity
fun.table <- read.table('./data/total_taxa_data/Pred/functional_prediction.txt', 
                        sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
fun.table <-data.frame(t(fun.table[,-189]))
R.fun <- specnumber(fun.table)
H.fun <- diversity(fun.table, 'shannon')
S.fun <- diversity(fun.table, 'simpson')

#metabolic diversity
#pathway 
path.table <- read.table('./data/total_taxa_data/Pred/pathway_prediction.txt', 
                         sep = '\t', header = T, row.names = 1, stringsAsFactors = F, quote = "")
path.table <- data.frame(t(path.table[,!colnames(path.table) %in% c('level1','level2', 'level3')]))
R.path <- specnumber(path.table)
H.path <- diversity(path.table, 'shannon')
S.path <- diversity(path.table, "simpson")


diver_table <- data.frame(alpha_diversity, df.pd, R.fun, H.fun, S.fun, R.path, H.path, S.path)
write.csv(diver_table,file = './results/tables/total_α_diver_table.csv')

env.table <- as(sample_data(water_physeq), 'data.frame')
sub.envtable <- cbind(env.table[,1:3], 
                      diver_table[,!colnames(diver_table) %in% c('chao1', 'se.chao1', 'se.ACE', 'ACE', 'InvSimpson','Fisher')])
geo.diver <- aggregate(sub.envtable[,-1], by = list(sub.envtable$Site), mean)

#alpha diversity distribution across geographic
library(scales)
library(ggplot2)
mytheme<-theme_bw()+
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size=12),
        legend.position = c(0.8, 0.2), 
        panel.grid = element_blank())

rescale(c(2500,1500,500),to=c(3,7))
bubble.chao1<-ggplot(data = geo.diver, aes(x = Longitude,y = Latitude,colour=factor(sign(Chao1)), size = Chao1)) +
  geom_point(aes(alpha=0.3)) +
  coord_map("polyconic") +
  xlab('Longitude')+ylab('Latitude')+
  scale_color_manual(values = c("dodgerblue","#E69F00"),guide = FALSE)+
  scale_size_continuous(breaks = c(500,1000,1000),
                        labels = c(500,1500,2500),range = c(3,7))+
  guides(size = guide_legend(override.aes = list(colour = list("dodgerblue","dodgerblue","dodgerblue"),alpha=0.3,
                                                 size=c(3,5,7))))+
  mytheme

bubble.Shannon <- ggplot()+
  geom_point(data = geo.diver, aes(x = Longitude,y = Latitude, size = Shannon),
             shape = 21, colour = "dodgerblue", 
             fill = "dodgerblue",alpha=0.3)+ 
  scale_size_area(max_size=7)+         
  coord_map("polyconic") +
  xlab('Longitude') + ylab('Latitude') +
  mytheme

rescale(c(1,0.85,0.7),to=c(3,7))
bubble.Simpson <- ggplot(data = geo.diver, aes(x = Longitude,y = Latitude,colour=factor(sign(Simpson)), size=Simpson)) +
  geom_point(aes(alpha=0.3)) +
  coord_map("polyconic") +
  xlab('Longitude')+ylab('Latitude')+
  scale_color_manual(values = c("dodgerblue","#E69F00"),guide = FALSE)+
  scale_size_continuous(breaks = c(0.7,0.15,0.15),
                        labels = c(0.7,0.85,1),range = c(3,7))+
  guides(size = guide_legend(override.aes = list(colour = list("dodgerblue","dodgerblue","dodgerblue"),alpha=0.5,
                                                 size=c(3,5,7))))+
  mytheme

rescale(c(150,100,50),to=c(3,7))
bubble.PD <- ggplot(data = geo.diver, aes(x = Longitude,y = Latitude,colour=factor(sign(PD)), size=PD)) +
  geom_point(aes(alpha=0.3)) +
  coord_map("polyconic") +
  xlab('Longitude')+ylab('Latitude')+
  scale_color_manual(values = c("dodgerblue","#E69F00"),guide = FALSE)+
  scale_size_continuous(breaks = c(50,50,50),
                        labels = c(50,100,150),range = c(3,7))+
  guides(size = guide_legend(override.aes = list(colour = list("dodgerblue","dodgerblue","dodgerblue"),alpha=0.5,
                                                 size=c(3,5,7))))+
  mytheme

rescale(c(8600,7800,7000),to=c(3,7))
bubble.R.fun <- ggplot(data = geo.diver, aes(x = Longitude,y = Latitude,colour=factor(sign(R.fun)), size=R.fun)) +
  geom_point(aes(alpha=0.3)) +
  coord_map("polyconic") +
  xlab('Longitude')+ylab('Latitude')+
  scale_color_manual(values = c("dodgerblue","#E69F00"),guide = FALSE)+
  scale_size_continuous(breaks = c(7000,800,800),
                        labels = c(7000,7800,8600),range = c(3,7))+
  guides(size = guide_legend(override.aes = list(colour = list("dodgerblue","dodgerblue","dodgerblue"),alpha=0.5,
                                                 size=c(3,5,7))))+
  mytheme

rescale(c(8,7.7,7.5),to=c(3,7))
bubble.H.fun <- ggplot(data = geo.diver, aes(x = Longitude,y = Latitude,colour=factor(sign(H.fun)), size=H.fun)) +
  geom_point(aes(alpha=0.3)) +
  coord_map("polyconic") +
  xlab('Longitude')+ylab('Latitude')+
  scale_color_manual(values = c("dodgerblue","#E69F00"),guide = FALSE)+
  scale_size_continuous(breaks = c(7.5,0.2,0.3),
                        labels = c(7.5,7.7,8),range = c(3,7))+
  guides(size = guide_legend(override.aes = list(colour = list("dodgerblue","dodgerblue","dodgerblue"),alpha=0.5,
                                                 size=c(3,4.6,7))))+
  mytheme

rescale(c(350,300,250),to=c(3,7))
bubble.R.path <- ggplot(data = geo.diver, aes(x = Longitude,y = Latitude,colour=factor(sign(R.path)), size=R.path)) +
  geom_point(aes(alpha=0.3)) +
  coord_map("polyconic") +
  xlab('Longitude')+ylab('Latitude')+
  scale_color_manual(values = c("dodgerblue","#E69F00"),guide = FALSE)+
  scale_size_continuous(breaks = c(250,50,50),
                        labels = c(250,300,350),range = c(3,7))+
  guides(size = guide_legend(override.aes = list(colour = list("dodgerblue","dodgerblue","dodgerblue"),alpha=0.5,
                                                 size=c(3,5,7))))+
  mytheme

rescale(c(4.3,4.15,4),to=c(3,7))
bubble.H.path <- ggplot(data = geo.diver, aes(x = Longitude,y = Latitude,colour=factor(sign(H.path)), size=H.path)) +
  geom_point(aes(alpha=0.3)) +
  coord_map("polyconic") +
  xlab('Longitude')+ylab('Latitude')+
  scale_color_manual(values = c("dodgerblue","#E69F00"),guide = FALSE)+
  scale_size_continuous(breaks = c(4,0.15,0.15),
                        labels = c(4,4.15,4.3),range = c(3,7))+
  guides(size = guide_legend(override.aes = list(colour = list("dodgerblue","dodgerblue","dodgerblue"),alpha=0.5,
                                                 size=c(3,5,7))))+
  mytheme

##combine multip plots into one
library(cowplot)
plot_grid(bubble.chao1, bubble.Shannon, bubble.PD,
          bubble.R.fun, bubble.H.fun,
          labels = 'auto', nrow = 2, label_x = .01, label_y = 0.99,
          hjust = 0, label_size=14,align = "v")

#compare the diversity between different vegetation type and permafrost type using kruskal.wallis test
veg_cat_diver <- data.frame(vegetable_type = env.table[ , c(5)], 
                            diver_table[,!colnames(diver_table) %in% 
                                          c('se.chao1', 'se.ACE', 'Observed', 'ACE', 'Simpson',
                                            'InvSimpson','Fisher', 'SR', 'S.fun', 'R.path',
                                            'H.path', 'S.path')])
hist(veg_cat_diver$Chao1)
hist(veg_cat_diver$Shannon)
hist(veg_cat_diver$PD)
hist(veg_cat_diver$R.fun)
hist(veg_cat_diver$H.fun)


kruskal.test(Chao1 ~ vegetable_type, data = veg_cat_diver)
pairwise.wilcox.test(veg_cat_diver$Chao1, veg_cat_diver$vegetable_type,
                     p.adjust.method = "BH")
kruskal.test(Shannon ~ vegetable_type, data = veg_cat_diver)
pairwise.wilcox.test(veg_cat_diver$Shannon, veg_cat_diver$vegetable_type,
                     p.adjust.method = "BH")
kruskal.test(PD ~ vegetable_type, data = veg_cat_diver)
pairwise.wilcox.test(veg_cat_diver$PD, veg_cat_diver$vegetable_type,
                     p.adjust.method = "BH")
kruskal.test(R.fun ~ vegetable_type, data = veg_cat_diver)
pairwise.wilcox.test(veg_cat_diver$R.fun, veg_cat_diver$vegetable_type,
                     p.adjust.method = "BH")
kruskal.test(H.fun ~ vegetable_type, data = veg_cat_diver)
pairwise.wilcox.test(veg_cat_diver$H.fun, veg_cat_diver$vegetable_type,
                     p.adjust.method = "BH")
veg_cat_diver_long <- melt(veg_cat_diver, id = (c('vegetable_type')))

diversity.veg <- ggplot(veg_cat_diver_long, aes(x = vegetable_type, y = value, fill = vegetable_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = vegetable_type), size = 1.4, alpha = 0.2) +
  facet_wrap(~ variable, nrow = 2, scale = 'free')+
  theme_bw() +
  theme(legend.position="none", 
        strip.text = element_text(size = rel(1.5)),
        strip.background = element_rect(fill = 'lightblue', colour = 'black'),
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12))


Chao1.veg <- ggplot(veg_cat_diver, aes(x = vegetable_type, y = Chao1, fill = vegetable_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="steelblue", size=1.4, alpha=0.2) +
  theme_classic() +
  labs(y="Chao1",x="Vegetable type")+
  theme(legend.position="none",
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12))

Shannon.veg <- ggplot(veg_cat_diver, aes(x = vegetable_type, y = Shannon, fill = vegetable_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="steelblue", size=1.4, alpha=0.2) +
  theme_classic() +
  labs(y="Shannon",x="Vegetable type")+
  theme(legend.position="none",
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12))

Chao1.veg <- ggplot(veg_cat_diver, aes(x = vegetable_type, y = Chao1, fill = vegetable_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="steelblue", size=1.4, alpha=0.2) +
  theme_classic() +
  labs(y="Chao1",x="Vegetable type")+
  theme(legend.position="none",
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12))


##diversity Correlation.R
# ============================================================
# Tutorial on drawing a correlation map using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================
avg_diversity <- read.csv('./results/tables/avg_diversity_table.csv', header = T, row.names = 1, stringsAsFactors = F)
#test correlation between env and diversity
library(lme4)
library(lmerTest)
library(MuMIn)
#Identify, describe, plot, and remove the outliers from the dataset
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
  if(response == "y" | response == "yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "n")
    return(invisible(dt))
  } else{
    cat("Nothing changed", "n")
    return(invisible(var_name))
  }
}
outlierKD(avg_diversity, nitrate_reduction)
yes

avg_diversity <- read.csv('./results/tables/avg_diversity_remove_outliners_table.csv', header = T, row.names = 1, stringsAsFactors = F)
write.csv(avg_diversity, './results/tables/avg_diversity_remove_outliners_table.csv')

outlierKD(avg_diversity, methanotrophy)
yes

vars <- colnames(avg_diversity)[-c(1,31:43)]

mode <- lmer(SR ~ SUVA254 + (1|Region), data = avg_diversity, na.action=na.omit)
summary(mode)$coefficients[2,1]
anova(mode)$Pr
 
vc <- vcov(mode)

# diagonal matrix of standard deviations associated with vcov
S <- sqrt(diag(diag(vc), nrow(vc), nrow(vc)))

# convert vc to a correlation matrix
solve(S) %*% vc %*% solve(S)

summary(mode)$Pr

selectMethod(print , "mer")
so <- summary(mode)
so@vcov@factors$correlation
cov2cor(vcov(mode))

#linner mixed model for matrix to matrix
lmm.mat.cal <- function(y, x, method){
  y <- as.matrix(y)
  x <- as.matrix(x)
  df<-NULL
  for(i in colnames(y)){
    for(j in colnames(x)){
      a <- y[, i, drop = F]
      b <- x[, j, drop = F]
      mode <- lmer(a ~ b + (1|Region), data = avg_diversity, na.action=na.omit)
      coeff <- summary(mode)$coefficients[2,1]
      r.square <- r.squaredGLMM(mode)[1]
      if (coeff>0) r = sqrt(r.square)
      else r = (-1) * sqrt(r.square)
      tmp <- c(i, j, r, r.square, anova(mode)$Pr)
      if(is.null(df)){
        df <- tmp  
      }
      else{
        df <- rbind(df, tmp)
      }    
    }
  }
  df<-data.frame(row.names=NULL,df)
  colnames(df)<-c("Diversity","Env","Correlation","r.square", "Pvalue")
  df$Pvalue<-as.numeric(as.character(df$Pvalue))
  df$AdjPvalue<-rep(0,dim(df)[1])
  df$Correlation<-as.numeric(as.character(df$Correlation))
  #You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
  # 1 -> donot adjust
  # 2 -> adjust Env + Type (column on the correlation plot)
  # 3 -> adjust Diversity + Type (row on the correlation plot for each type)
  # 4 -> adjust Diversity (row on the correlation plot)
  # 5 -> adjust Env (panel on the correlation plot)
  adjustment_label<-c("NoAdj","AdjEnvAndType","AdjDiversityAndType","AdjDiversity","AdjEnv")
  adjustment<-5
  if(adjustment==1){
    df$AdjPvalue<-df$Pvalue
  } else if (adjustment==2){
    for(i in unique(df$Env)){
      for(j in unique(df$Type)){
        sel<-df$Env==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
  } else if (adjustment==3){
    for(i in unique(df$Diversity)){
      for(j in unique(df$Type)){
        sel<-df$Diversity==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
  } else if (adjustment==4){
    for(i in unique(df$Diversity)){
      sel<-df$Diversity==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  } else if (adjustment==5){
    for(i in unique(df$Env)){
      sel<-df$Env==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
  #Now we generate the labels for signifant values
  df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  df$Diversity <-factor(df$Diversity, ordered = T, levels = rev(colnames(y)))
  df$Env <-factor(df$Env, ordered = T, levels = colnames(x))
  return(df)
}

diver_table <- avg_diversity[ ,c("SR", "Shannon", "Simpson", "PD", "R.fun", 'H.fun', "S.fun", 'R.path','H.path', 'S.path',
                                 "methaon_rel_abun", "type1_methano_rel_abun", "type2_methano_rel_abun")]
substrate_table <- avg_diversity[ ,c("DOC", "TN", "TC", "DON", "NH4_N", "NO3_N", "DIN", "S275_295", "SUVA254", 
                                     "a300", "hix", "bix", "fi", "C1...", "C2...", "C3...", 
                                     "C4...", "MAP", "MAT", "Temp", "Depth", "DO", "pH",
                                     "Conductivity", "Salinity", "Ca", "Mg", "K", "Na")]
metabolism_table <- avg_diversity[,44:55]

c("chitinolysis","cellulolysis","fermentation","methanogenesis","methanotrophy",
  "methylotrophy","nitrogen_fixation","aerobic_ammonia_oxidation","nitrification",
  "aerobic_nitrite_oxidation", "nitrate_reduction","nitrate_respiration", "nitrite_respiration")


lmm.matrix <- lmm.mat.cal(diver_table, substrate_table)
lmm.matrix$Correlation[lmm.matrix$AdjPvalue >= 0.05] <- 0

lmm.matrix1 <- lmm.mat.cal(metabolism_table, substrate_table )
lmm.matrix1$Correlation[lmm.matrix1$AdjPvalue >= 0.05] <- 0

#plot
p <- ggplot(aes(x = Env, y = Diversity, fill = Correlation), data = lmm.matrix1)+
  geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C")+
  geom_text(aes(label=Significance), color="black", size = 5.5)+
  labs(y = 'Diversity index', x = 'Environmental factors', fill= 'Correlation') +
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
print(p)

multi.plot <- function(a, b){
  tmp_plot <- ggplot(avg_diversity, aes(x = a, y = b))+
    geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
    geom_smooth(method="lm", size=1.5, se=T,colour='black') +
    #scale_y_continuous(expand = c(0,0))+
    #scale_x_continuous(expand = c(0,0))+
    theme(axis.line=element_line(colour = 'black'))+
    theme_classic()
  return(tmp_plot)
}


DOC_R.fun <- multi.plot(avg_diversity$DOC, avg_diversity$R.fun)
DOC_H.fun <- multi.plot(avg_diversity$DOC, avg_diversity$H.fun)
DOC_S.fun <- multi.plot(avg_diversity$DOC, avg_diversity$S.fun)
DOC_R.path <- multi.plot(avg_diversity$DOC, avg_diversity$R.path)
DOC_H.path <- multi.plot(avg_diversity$DOC, avg_diversity$H.path)
DOC_S.path <- multi.plot(avg_diversity$DOC, avg_diversity$S.path)

TN_R.fun <- multi.plot(avg_diversity$TN, avg_diversity$R.fun)
TN_H.fun <- multi.plot(avg_diversity$TN, avg_diversity$H.fun)
TN_S.fun <- multi.plot(avg_diversity$TN, avg_diversity$S.fun)
TN_R.path <- multi.plot(avg_diversity$TN, avg_diversity$R.path)
TN_H.path <- multi.plot(avg_diversity$TN, avg_diversity$H.path)
TN_S.path <- multi.plot(avg_diversity$TN, avg_diversity$S.path)


a300_R.fun <- multi.plot(avg_diversity$a300, avg_diversity$R.fun)
a300_H.fun <- multi.plot(avg_diversity$a300, avg_diversity$H.fun)
a300_S.fun <- multi.plot(avg_diversity$a300, avg_diversity$S.fun)
a300_R.path <- multi.plot(avg_diversity$a300, avg_diversity$R.path)
a300_H.path <- multi.plot(avg_diversity$a300, avg_diversity$H.path)
a300_S.path <- multi.plot(avg_diversity$a300, avg_diversity$S.path)

fi_R.fun <- multi.plot(avg_diversity$fi, avg_diversity$R.fun)
fi_H.fun <- multi.plot(avg_diversity$fi, avg_diversity$H.fun)
fi_S.fun <- multi.plot(avg_diversity$fi, avg_diversity$S.fun)
fi_R.path <- multi.plot(avg_diversity$fi, avg_diversity$R.path)
fi_H.path <- multi.plot(avg_diversity$fi, avg_diversity$H.path)
fi_S.path <- multi.plot(avg_diversity$fi, avg_diversity$S.path)

C1_R.fun <- multi.plot(avg_diversity$C1..., avg_diversity$R.fun)
C1_H.fun <- multi.plot(avg_diversity$C1..., avg_diversity$H.fun)
C1_S.fun <- multi.plot(avg_diversity$C1..., avg_diversity$S.fun)
C1_R.path <- multi.plot(avg_diversity$C1..., avg_diversity$R.path)
C1_H.path <- multi.plot(avg_diversity$C1..., avg_diversity$H.path)
C1_S.path <- multi.plot(avg_diversity$C1..., avg_diversity$S.path)

C2_R.fun <- multi.plot(avg_diversity$C2..., avg_diversity$R.fun)
C2_H.fun <- multi.plot(avg_diversity$C2..., avg_diversity$H.fun)
C2_S.fun <- multi.plot(avg_diversity$C2..., avg_diversity$S.fun)
C2_R.path <- multi.plot(avg_diversity$C2..., avg_diversity$R.path)
C2_H.path <- multi.plot(avg_diversity$C2..., avg_diversity$H.path)
C2_S.path <- multi.plot(avg_diversity$C2..., avg_diversity$S.path)

C3_R.fun <- multi.plot(avg_diversity$C3..., avg_diversity$R.fun)
C3_H.fun <- multi.plot(avg_diversity$C3..., avg_diversity$H.fun)
C3_S.fun <- multi.plot(avg_diversity$C3..., avg_diversity$S.fun)
C3_R.path <- multi.plot(avg_diversity$C3..., avg_diversity$R.path)
C3_H.path <- multi.plot(avg_diversity$C3..., avg_diversity$H.path)
C3_S.path <- multi.plot(avg_diversity$C3..., avg_diversity$S.path)


library(cowplot)
pdf('./results/figs/env_diversity.pdf', width = 12, height = 12)
multi_plot1 <- plot_grid(DOC_R.fun,DOC_H.fun,DOC_S.fun,TN_R.fun,TN_H.fun,TN_S.fun,
                         DOC_R.path,DOC_H.path,DOC_S.path,TN_R.path,TN_H.path,TN_S.path,
                        labels = 'auto', ncol = 3, nrow = 4, byrow = T,
                        label_x = .01, label_y = 0.99, hjust = 0, label_size=14,align = "v")
print(multi_plot1)
dev.off()

pdf('./results/figs/light_diversity.pdf', width = 20, height = 18)
multi_plot2 <- plot_grid(a300_R.fun,fi_R.fun,C1_R.fun,C2_R.fun,C3_R.fun,
                         a300_H.fun,fi_H.fun,C1_H.fun,C2_H.fun,C3_H.fun,
                         a300_S.fun,fi_S.fun,C1_S.fun,C2_S.fun,C3_S.fun,
                         a300_R.path,fi_R.path,C1_R.path,C2_R.path,C3_R.path,
                         a300_H.path,fi_H.path,C1_H.path,C2_H.path,C3_H.path,
                         a300_S.path,fi_S.path,C1_S.path,C2_S.path,C3_S.path,
                         labels = 'auto', ncol = 5, nrow = 6, byrow = T,
                         label_x = .01, label_y = 0.99, hjust = 0, label_size=14,align = "v")

print(multi_plot2)
dev.off()

View(lmm.matrix)

model1 <- lapply(vars, function(x) {
  lmer(substitute(SR ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model1, summary)


model2 <- lapply(vars, function(x) {
  lmer(substitute(Shannon ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model2, summary)

model3 <- lapply(vars, function(x) {
  lmer(substitute(Simpson ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model3, summary)

model4 <- lapply(vars, function(x) {
  lmer(substitute(PD ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model4, summary)

model5 <- lapply(vars, function(x) {
  lmer(substitute(R.fun ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model5, summary)

model6 <- lapply(vars, function(x) {
  lmer(substitute(H.fun ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model6, summary)

model7 <- lapply(vars, function(x) {
  lmer(substitute(S.fun ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model7, summary)

model8 <- lapply(vars, function(x) {
  lmer(substitute(R.path ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model8, summary)

model9 <- lapply(vars, function(x) {
  lmer(substitute(H.path ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model9, summary)

model10 <- lapply(vars, function(x) {
  lmer(substitute(S.path ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model10, summary)

model11 <- lapply(vars, function(x) {
  lmer(substitute(methaon_rel_abun ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model11, summary)

model12 <- lapply(vars, function(x) {
  lmer(substitute(type1_methano_rel_abun ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model12, summary)

model13 <- lapply(vars, function(x) {
  lmer(substitute(type2_methano_rel_abun ~ i + (1|Region),list(i = as.name(x))), data = avg_diversity)})
lapply(model13, summary)

S275_295_SR <- multi.plot(avg_diversity$S275_295, avg_diversity$SR)
S275_295_Shannon <- multi.plot(avg_diversity$S275_295, avg_diversity$Shannon)
S275_295_Simpson <- multi.plot(avg_diversity$S275_295, avg_diversity$Simpson)
S275_295_PD <- multi.plot(avg_diversity$S275_295, avg_diversity$PD)
S275_295_R.fun <- multi.plot(avg_diversity$S275_295, avg_diversity$R.fun)
S275_295_S.fun <- multi.plot(avg_diversity$S275_295, avg_diversity$S.fun)
S275_295_methaon_rel_abun <- multi.plot(avg_diversity$S275_295, avg_diversity$methaon_rel_abun)
S275_295_type1_methano_rel_abun <- multi.plot(avg_diversity$S275_295, avg_diversity$type1_methano_rel_abun)
S275_295_type2_methano_rel_abun <- multi.plot(avg_diversity$S275_295, avg_diversity$type2_methano_rel_abun)


SUVA254_SR <- multi.plot(avg_diversity$SUVA254, avg_diversity$SR)
SUVA254_Shannon <- multi.plot(avg_diversity$SUVA254, avg_diversity$Shannon)
SUVA254_Simpson <- multi.plot(avg_diversity$SUVA254, avg_diversity$Simpson)
SUVA254_PD <- multi.plot(avg_diversity$SUVA254, avg_diversity$PD)
SUVA254_R.fun <- multi.plot(avg_diversity$SUVA254, avg_diversity$R.fun)
SUVA254_S.fun <- multi.plot(avg_diversity$SUVA254, avg_diversity$S.fun)
SUVA254_methaon_rel_abun <- multi.plot(avg_diversity$SUVA254, avg_diversity$methaon_rel_abun)
SUVA254_type1_methano_rel_abun <- multi.plot(avg_diversity$SUVA254, avg_diversity$type1_methano_rel_abun)
SUVA254_type2_methano_rel_abun <- multi.plot(avg_diversity$SUVA254, avg_diversity$type2_methano_rel_abun)

a300_SR <- multi.plot(avg_diversity$a300, avg_diversity$SR)
a300_Shannon <- multi.plot(avg_diversity$a300, avg_diversity$Shannon)
a300_Simpson <- multi.plot(avg_diversity$a300, avg_diversity$Simpson)
a300_PD <- multi.plot(avg_diversity$a300, avg_diversity$PD)
a300_R.fun <- multi.plot(avg_diversity$a300, avg_diversity$R.fun)
a300_S.fun <- multi.plot(avg_diversity$a300, avg_diversity$S.fun)
a300_methaon_rel_abun <- multi.plot(avg_diversity$a300, avg_diversity$methaon_rel_abun)
a300_type1_methano_rel_abun <- multi.plot(avg_diversity$a300, avg_diversity$type1_methano_rel_abun)
a300_type2_methano_rel_abun <- multi.plot(avg_diversity$a300, avg_diversity$type2_methano_rel_abun)

hix_SR <- multi.plot(avg_diversity$hix, avg_diversity$SR)
hix_Shannon <- multi.plot(avg_diversity$hix, avg_diversity$Shannon)
hix_Simpson <- multi.plot(avg_diversity$hix, avg_diversity$Simpson)
hix_PD <- multi.plot(avg_diversity$hix, avg_diversity$PD)
hix_R.fun <- multi.plot(avg_diversity$hix, avg_diversity$R.fun)
hix_S.fun <- multi.plot(avg_diversity$hix, avg_diversity$S.fun)
hix_methaon_rel_abun <- multi.plot(avg_diversity$hix, avg_diversity$methaon_rel_abun)
hix_type1_methano_rel_abun <- multi.plot(avg_diversity$hix, avg_diversity$type1_methano_rel_abun)
hix_type2_methano_rel_abun <- multi.plot(avg_diversity$hix, avg_diversity$type2_methano_rel_abun)

bix_SR <- multi.plot(avg_diversity$bix, avg_diversity$SR)
bix_Shannon <- multi.plot(avg_diversity$bix, avg_diversity$Shannon)
bix_Simpson <- multi.plot(avg_diversity$bix, avg_diversity$Simpson)
bix_PD <- multi.plot(avg_diversity$bix, avg_diversity$PD)
bix_R.fun <- multi.plot(avg_diversity$bix, avg_diversity$R.fun)
bix_S.fun <- multi.plot(avg_diversity$bix, avg_diversity$S.fun)
bix_methaon_rel_abun <- multi.plot(avg_diversity$bix, avg_diversity$methaon_rel_abun)
bix_type1_methano_rel_abun <- multi.plot(avg_diversity$bix, avg_diversity$type1_methano_rel_abun)
bix_type2_methano_rel_abun <- multi.plot(avg_diversity$bix, avg_diversity$type2_methano_rel_abun)

fi_SR <- multi.plot(avg_diversity$fi, avg_diversity$SR)
fi_Shannon <- multi.plot(avg_diversity$fi, avg_diversity$Shannon)
fi_Simpson <- multi.plot(avg_diversity$fi, avg_diversity$Simpson)
fi_PD <- multi.plot(avg_diversity$fi, avg_diversity$PD)
fi_R.fun <- multi.plot(avg_diversity$fi, avg_diversity$R.fun)
fi_S.fun <- multi.plot(avg_diversity$fi, avg_diversity$S.fun)
fi_methaon_rel_abun <- multi.plot(avg_diversity$fi, avg_diversity$methaon_rel_abun)
fi_type1_methano_rel_abun <- multi.plot(avg_diversity$fi, avg_diversity$type1_methano_rel_abun)
fi_type2_methano_rel_abun <- multi.plot(avg_diversity$fi, avg_diversity$type2_methano_rel_abun)


C1_SR <- multi.plot(avg_diversity$C1..., avg_diversity$SR)
C1_Shannon <- multi.plot(avg_diversity$C1..., avg_diversity$Shannon)
C1_Simpson <- multi.plot(avg_diversity$C1..., avg_diversity$Simpson)
C1_PD <- multi.plot(avg_diversity$C1..., avg_diversity$PD)
C1_R.fun <- multi.plot(avg_diversity$C1..., avg_diversity$R.fun)
C1_S.fun <- multi.plot(avg_diversity$C1..., avg_diversity$S.fun)
C1_methaon_rel_abun <- multi.plot(avg_diversity$C1..., avg_diversity$methaon_rel_abun)
C1_type1_methano_rel_abun <- multi.plot(avg_diversity$C1..., avg_diversity$type1_methano_rel_abun)
C1_type2_methano_rel_abun <- multi.plot(avg_diversity$C1..., avg_diversity$type2_methano_rel_abun)

C2_SR <- multi.plot(avg_diversity$C2..., avg_diversity$SR)
C2_Shannon <- multi.plot(avg_diversity$C2..., avg_diversity$Shannon)
C2_Simpson <- multi.plot(avg_diversity$C2..., avg_diversity$Simpson)
C2_PD <- multi.plot(avg_diversity$C2..., avg_diversity$PD)
C2_R.fun <- multi.plot(avg_diversity$C2..., avg_diversity$R.fun)
C2_S.fun <- multi.plot(avg_diversity$C2..., avg_diversity$S.fun)
C2_methaon_rel_abun <- multi.plot(avg_diversity$C2..., avg_diversity$methaon_rel_abun)
C2_type1_methano_rel_abun <- multi.plot(avg_diversity$C2..., avg_diversity$type1_methano_rel_abun)
C2_type2_methano_rel_abun <- multi.plot(avg_diversity$C2..., avg_diversity$type2_methano_rel_abun)

C3_SR <- multi.plot(avg_diversity$C3..., avg_diversity$SR)
C3_Shannon <- multi.plot(avg_diversity$C3..., avg_diversity$Shannon)
C3_Simpson <- multi.plot(avg_diversity$C3..., avg_diversity$Simpson)
C3_PD <- multi.plot(avg_diversity$C3..., avg_diversity$PD)
C3_R.fun <- multi.plot(avg_diversity$C3..., avg_diversity$R.fun)
C3_S.fun <- multi.plot(avg_diversity$C3..., avg_diversity$S.fun)
C3_methaon_rel_abun <- multi.plot(avg_diversity$C3..., avg_diversity$methaon_rel_abun)
C3_type1_methano_rel_abun <- multi.plot(avg_diversity$C3..., avg_diversity$type1_methano_rel_abun)
C3_type2_methano_rel_abun <- multi.plot(avg_diversity$C3..., avg_diversity$type2_methano_rel_abun)


C4_SR <- multi.plot(avg_diversity$C4..., avg_diversity$SR)
C4_Shannon <- multi.plot(avg_diversity$C4..., avg_diversity$Shannon)
C4_Simpson <- multi.plot(avg_diversity$C4..., avg_diversity$Simpson)
C4_PD <- multi.plot(avg_diversity$C4..., avg_diversity$PD)
C4_R.fun <- multi.plot(avg_diversity$C4..., avg_diversity$R.fun)
C4_S.fun <- multi.plot(avg_diversity$C4..., avg_diversity$S.fun)
C4_methaon_rel_abun <- multi.plot(avg_diversity$C4..., avg_diversity$methaon_rel_abun)
C4_type1_methano_rel_abun <- multi.plot(avg_diversity$C4..., avg_diversity$type1_methano_rel_abun)
C4_type2_methano_rel_abun <- multi.plot(avg_diversity$C4..., avg_diversity$type2_methano_rel_abun)


library(cowplot)
multi_plot <- plot_grid(S275_295_SR,
                       S275_295_Shannon,
                       S275_295_Simpson,
                       S275_295_PD,
                       S275_295_R.fun,
                       S275_295_S.fun,
                       S275_295_methaon_rel_abun,
                       S275_295_type1_methano_rel_abun,
                       S275_295_type2_methano_rel_abun,
                       SUVA254_SR,
                       SUVA254_Shannon,
                       SUVA254_Simpson,
                       SUVA254_PD,
                       SUVA254_R.fun,
                       SUVA254_S.fun,
                       SUVA254_methaon_rel_abun,
                       SUVA254_type1_methano_rel_abun,
                       SUVA254_type2_methano_rel_abun,
                       a300_SR,
                       a300_Shannon,
                       a300_Simpson,
                       a300_PD,
                       a300_R.fun,
                       a300_S.fun,
                       a300_methaon_rel_abun,
                       a300_type1_methano_rel_abun,
                       a300_type2_methano_rel_abun,
                       labels = 'auto', ncol = 9, nrow = 3, byrow = T,
                       label_x = .01, label_y = 0.99, hjust = 0, label_size=14,align = "v")




#linner mixed model plot
p3 <- ggplot()+
  geom_point(data = Q10final,aes(x = OCCa, y = Q2), 
             shape=17, size=5,colour='#00AB94',alpha=1)+
  geom_smooth(data = plot_interval_data,
              aes(x = OCCa,
                  y = pred_R, 
                  ymin = lb_R, 
                  ymax = ub_R,
              ), color = '#00AB94',fill="#00AB94", alpha=0.3,
              stat = "identity") +
  scale_y_continuous(limits = c(1.2, 3.2)) +
  scale_x_continuous(limits = c(-0.6, 0.9))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
  )
print(p3)
plot(Shannon ~ fi, avg_diversity)

mode2<-lmer(Shannon ~ fi + (1|Region), avg_diversity)
r.squaredGLMM (mode2)
library(car)
influencePlot(mode2)
mode3 <- lmer(Shannon ~ fi + (1|Region), avg_diversity[!rownames(avg_diversity) %in% c('Z6', 'Z8', 'Z22', 'Z29'), ])
influencePlot(mode3)
r.squaredGLMM(mode3)
mode4 <- lmer(Shannon ~ fi + (1|Region), avg_diversity[!rownames(avg_diversity) %in% c('Z5', 'Z22', 'Z9', 'Z8', 'Z20', 'Z29'), ])
influencePlot(mode4)
r.squaredGLMM(mode4)
mode5 <- lmer(Shannon ~ fi + (1|Region), avg_diversity[!rownames(avg_diversity) %in% c('Z5', 'Z22', 'Z9', 'Z6', 'Z8', 'Z20', 'Z29',
                                                                                      'Z19', 'Z21', 'Z23', 'Z37', 'Z43'), ])
influencePlot(mode5)

mode5 <- lmer(Shannon ~ fi + (1|Region), avg_diversity[!rownames(avg_diversity) %in% c('Z5', 'Z22', 'Z8', 'Z19', 'Z20', 'Z29', 'Z43', 'Z23', 'Z9'), ])

r.squaredGLMM(mode5)

#test  the colliner relationshape between all factors
library(corrplot)
corr_mat = cor(x, method="pearson")
corr_mat[1:5,1:5]
res <- cor.mtest(x, conf.level = 0.95)

corrplot(corr_mat,addCoef.col = "black", number.digits = 2, type = "upper",
         number.cex = 0.65,tl.col = "black", tl.cex = 1, cl.cex = 1,
         addrect = 2, col = terrain.colors(100))

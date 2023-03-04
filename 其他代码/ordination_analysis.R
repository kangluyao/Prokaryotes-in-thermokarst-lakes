library(vegan)
library(GUniFrac)
library(ggplot2)
library(phangorn)
library(phyloseq)
#taxonomic analysis
total.otu.table <- data.frame(t(otu_table(water.rel)))
total.otu.dis <- vegdist(sqrt(total.otu.table))

#phylogenetic analysis
unifracs <- GUniFrac(t(as.matrix(otu_table(water_physeq))), 
                     midpoint(phy_tree(water_physeq)), alpha = c(0, 0.5, 1))$unifracs
# We can extract a variety of distance matrices with different weightings.
dw <- unifracs[, , "d_1"]  # Weighted UniFrac
du <- unifracs[, , "d_UW"]  # Unweighted UniFrac

#functional analysis
fun.table <- read.table('./data/total_taxa_data/Pred/functional_prediction.txt', 
                        sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
fun.table <-data.frame(t(fun.table[,-189]))
fun.dis <- vegdist(sqrt(fun.table))
#pathway 
path.table <- read.table('./results/tables/res_spe_func_perc.csv', 
                         header = T, row.names = 1, stringsAsFactors = F)
#path.table <- data.frame(t(path.table[,!colnames(path.table) %in% c('level1','level2', 'level3')]))
path.dis <- vegdist(sqrt(t(path.table)))
#read in metadata
env.table <- as(sample_data(water_physeq), 'data.frame')
env.dis <- vegdist(env.table[ , -c(1:5)], "euclid", na.rm = T)

##Test and forward selection of environmental variables
Env_select <- function (Dist_Matrix,Env,Number_Permutations=999) {
  mod1<-capscale(Dist_Matrix~.,Env,add=T)
  mod0<-capscale(Dist_Matrix~1,Env,add=T)
  mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999,trace = F)
  Env.P<-anova.cca(mod,permutations=Number_Permutations)[[4]][1]
  Env_se<-mod$CCA$biplot
  Env_se<-Env[,rownames(Env_se)]
  return(Env_se)
}
env_vars <- data.frame(scale(env.table[ ,-c(1:5)]))
env_se_total <- Env_select (total.otu.dis, env_vars)
env_se_unifrac <- Env_select(du, env_vars)
env_se_fun <- Env_select(fun.dis, env_vars)
env_se_path <- Env_select(path.dis, env_vars)


#==========total microbiomes analysis=================
#toxonomic analysis
ord.total <-  cmdscale(total.otu.dis,  k = 2, eig = T, add = T)
round(ord.total$eig*100/sum(ord.total$eig),1)[c(1,2)]

fit.total <- envfit(ord.total, cbind(env.table[,c('vegetable_type')], env_se_total), perm=999, na.rm = TRUE)
fit.total
data.scores.total = as.data.frame(scores(ord.total))
data.scores.total <- cbind(vegetable_type = env.table[,c('vegetable_type')], data.scores.total)
fit_coord_cont.total = as.data.frame(scores(fit.total, "vectors")) 
fit_coord_cat.total = as.data.frame(scores(fit.total, "factors")) 

gg.total <- ggplot(data = data.scores.total, aes(x = Dim1, y = Dim2)) + 
  geom_point(data = data.scores.total, colour = '#E76BF3', size = 3.5, alpha = 0.8) + 
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               data = fit_coord_cont.total, size =1, alpha = 0.5, colour = "grey30") +
  labs(x = 'Dim1 (11.9%)', y = 'Dim (9.1%)') + 
  geom_text(data = fit_coord_cont.total, aes(x = Dim1, y = Dim2), colour = "black", 
            size = 3.5, label = row.names(fit_coord_cont.total)) + 
  theme(axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))

#phylogenetic analysis
# use vegan's cmdscale function to make a PCoA ordination from a distance matrix.
pcoa_dw <- cmdscale(dw, k = nrow(t(as.matrix(otu_table(water_physeq)))) - 1, eig = TRUE, add = TRUE)
pcoa_du <- cmdscale(du, k = nrow(t(as.matrix(otu_table(water_physeq)))) - 1, eig = TRUE, add = TRUE)
round(pcoa_du$eig*100/sum(pcoa_du$eig),1)[c(1,2)]

fit.phy <- envfit(pcoa_du, env_se_unifrac, perm=999, na.rm = TRUE)
fit.phy
data.scores.phy = as.data.frame(scores(pcoa_du))

fit_coord_cont.phy = as.data.frame(scores(fit.phy, "vectors")) 
fit_coord_cat.phy = as.data.frame(scores(fit.phy, "factors"))

gg.phy <- ggplot(data = data.scores.phy, aes(x = Dim1, y = Dim2)) + 
  geom_point(data = data.scores.phy, colour = '#00BF7D', size = 3.5, alpha = 0.8) + 
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               data = fit_coord_cont.phy, size =1, alpha = 0.5, colour = "grey30") +
  labs(x = 'Dim1 (11.2%)', y = 'Dim (6.7%)') + 
  geom_text(data = fit_coord_cont.phy, aes(x = Dim1, y = Dim2), colour = "black", 
            size = 3.5, label = row.names(fit_coord_cont.phy)) + 
  theme(axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))

#functional analysis
ord.fun <- cmdscale(fun.dis,  k = 2, eig = T, add = T)
round(ord.fun$eig*100/sum(ord.fun$eig),1)[c(1,2)]
fit.fun <- envfit(ord.fun, env_se_fun, perm=999, na.rm = TRUE)
fit.fun
data.scores.fun = as.data.frame(scores(ord.fun))

fit_coord_cont.fun = as.data.frame(scores(fit.fun, "vectors"))
fit_coord_cat.fun = as.data.frame(scores(fit.fun, "factors"))

gg.fun <- ggplot(data = data.scores.fun, aes(x = Dim1, y = Dim2)) + 
  geom_point(data = data.scores.fun, colour = '#F8766D', size = 3.5, alpha = 0.8) + 
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               data = fit_coord_cont.fun, size =1, alpha = 0.5, colour = "grey30") +
  labs(x = 'Dim1 (26.4%)', y = 'Dim (11.2%)') + 
  geom_text(data = fit_coord_cont.fun, aes(x = Dim1, y = Dim2), colour = "black", 
            size = 3.5, label = row.names(fit_coord_cont.fun)) + 
  theme(axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))

#metabolism analysis
ord.path <- cmdscale(path.dis,  k = 2, eig = T, add = T)
round(ord.path$eig*100/sum(ord.path$eig),1)[c(1,2)]
fit.path <- envfit(ord.path, env_se_path, perm=999, na.rm = TRUE)
fit.path
data.scores.path = as.data.frame(scores(ord.path))

fit_coord_cont.path = as.data.frame(scores(fit.path, "vectors"))
fit_coord_cat.path = as.data.frame(scores(fit.path, "factors")) 

gg.path <- ggplot(data = data.scores.path, aes(x = Dim1, y = Dim2)) + 
  geom_point(data = data.scores.path, colour = '#DD5F60', size = 4, alpha = 0.8) + 
  #scale_colour_manual(values = c("orange", "steelblue"))  + 
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               data = fit_coord_cont.path, size =1, alpha = 0.5, colour = "black") +
  #geom_point(data = fit_coord_cat.path, aes(x = Dim1, y = Dim2), 
             #shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  #geom_text(data = fit_coord_cat.path, aes(x = Dim1, y = Dim2+0.04), 
  #label = row.names(fit_coord_cat.path), colour = "navy", fontface = "bold") + 
  geom_text(data = fit_coord_cont.path, aes(x = Dim1, y = Dim2), colour = "black", 
            size = 3.5, label = row.names(fit_coord_cont.path)) + 
  theme(axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))

library(cowplot)
ggdraw() +
  draw_plot(gg.total, x =0, y = 1/2, width = 1/2, height = 1/2) +
  draw_plot(gg.phy, x = 1/2, y = 1/2, width = 1/2, height = 1/2) +
  draw_plot(gg.fun, x = 0, y = 0, width = 1/2, height = 1/2) +
  draw_plot(gg.path, x = 1/2, y = 0, width = 1/2, height = 1/2) +
  draw_plot_label(label = c("(a)", "(b)", '(c)', '(d)'), size = 11,
                  x = c(0, 1/2, 0, 1/2), y = c(1, 1, 1/2, 1/2))

mantel.partial(xdis, ydis, zdis, method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))

#fig2
fig2 <- plot_grid(p_taxonomic_distance, p_phylogenetic_distance, p_functional_distance, 
          gg.total, gg.phy, gg.fun, labels = 'auto', ncol = 3, nrow = 2, 
          label_x = .01, label_y = 0.99, hjust = 0, label_size=14,align = "v")
ggsave('./figs/ordination_distacne.pdf', fig2, width = 14, height = 7, dpi = 300)

#variance partition
dom_vars <- c('S275_295', 'SUVA254', 'a300', 'DON', 'TOC')
pH_vars <- c('pH', 'Temp', 'Depth', 'DO')
climate_vars <- c('MAP', 'MAT')
conductivity_vars <- c('K', 'Ca', 'Na', 'Mg', 'Conductivity')
env.table <- data.frame(env.table)

total.otu.rel <- data.frame(total.otu.rel)
varpart(total.otu.rel, ~ S275_295+SUVA254+a300+DON+TOC, 
        ~pH+Temp+Depth+DO, ~ K+Ca+Na+Mg+Conductivity, data = env.table, sqrt.dist = T, 999)
varpart(total.otu.dis, env.table[ ,climate_vars], env.table[ ,dom_vars], env.table[ ,pH_vars], env.table[ ,conductivity_vars])
varpart(du, env.table[ ,climate_vars], env.table[ ,dom_vars], env.table[ ,pH_vars], env.table[ ,conductivity_vars])
varpart(fun.dis, env.table[ ,climate_vars], env.table[ ,dom_vars], env.table[ ,pH_vars], env.table[ ,conductivity_vars])

varpart(total.otu.dis, env.table[ ,climate_vars], env.table[ ,dom_vars], env.table[ ,physicochemical_vars])


#==== core microbiomes analysis ========
#core microbiome are difined as the OTUs with relative abundance ≥ 0.01% and prevalence ≥ 10 accross 188 samples
water.rel <- microbiome::transform(water_physeq, "compositional")
core.phy<- core(water.rel, detection = 0.01/100, prevalence = 10/188, include.lowest = TRUE)
core.phy
core.otu.rel <-otu_table(core.phy)


core.otu.table <- t(otu_table(core.phy))
core.otu.dis <- vegdist(core.otu.table)

ord.core <- cmdscale(core.otu.dis)
fit.core <- envfit(ord.core, env.table, perm=999, na.rm = TRUE)

data.scores.core = as.data.frame(scores(ord.core))

fit_coord_cont.core = as.data.frame(scores(fit.core, "vectors")) * ordiArrowMul(fit.core)
fit_coord_cat.core = as.data.frame(scores(fit.core, "factors")) * ordiArrowMul(fit.core)
data.scores.core$pH = env.table$pH

gg.core <- ggplot(data = data.scores.core, aes(x = Dim1, y = Dim2)) + 
  geom_point(data = data.scores.core, colour = 'steelblue', size = 3, alpha = 0.5) + 
  #scale_colour_gradient(low = "orange", high = "steelblue")  + 
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               data = fit_coord_cont.core, size =1, alpha = 0.5, colour = "grey30") +
  geom_point(data = fit_coord_cat.core, aes(x = Dim1, y = Dim2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = fit_coord_cat.core, aes(x = Dim1, y = Dim2+0.04), 
            label = row.names(fit_coord_cat.core), colour = "navy", fontface = "bold") + 
  geom_text(data = fit_coord_cont.core, aes(x = Dim1, y = Dim2), colour = "grey30", 
            fontface = "bold", label = row.names(fit_coord_cont.core)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30"))



env.vars <- as(env.table[ ,c(2:8,12:17)], 'matrix')
sub.vars <- as(env.table[ ,c(9,10,11)], 'matrix')

mantel(env.dis, core.otu.dis,method="pearson", permutations=999, strata = NULL,
       na.rm = FALSE, parallel = getOption("mc.cores"))
env.vars <- as(env.table[ ,c(2:8,12:17)], 'matrix')
sub.vars <- as(env.table[ ,c(9,10,11)], 'matrix')
total.otu.dis

#toxonomic analysis
total.otu.rel <- t(otu_table(water_physeq))
total.otu.dis <- vegdist(sqrt(total.otu.rel))

# create environmental distance matrices
DOC.dist = vegdist(scale(env.table$DOC), 'euclidean')
TC.dist = vegdist(scale(env.table$TC), 'euclidean')
TN.dist = vegdist(scale(env.table$TN), 'euclidean')
NH4_N.dist = vegdist(scale(env.table$NH4_N), 'euclidean')
NO3_N.dist = vegdist(scale(env.table$NO3_N), 'euclidean')
DIN.dist = vegdist(scale(env.table$DIN), 'euclidean')
DON.dist = vegdist(scale(env.table$DON), 'euclidean')
S275_295.dist = vegdist(scale(env.table$S275_295), 'euclidean')
SUVA254.dist = vegdist(scale(env.table$SUVA254), 'euclidean')
bix.dist = vegdist(scale(env.table$bix), 'euclidean')
fi.dist = vegdist(scale(env.table$fi), 'euclidean')
hix.dist = vegdist(scale(env.table$hix), 'euclidean')
a300.dist = vegdist(scale(env.table$a300), 'euclidean')
pH.dist = vegdist(scale(env.table$pH), 'euclidean')
DO.dist = vegdist(scale(env.table$DO), 'euclidean')
Conductivity.dist = vegdist(scale(env.table$Conductivity), 'euclidean')
Salinity.dist = vegdist(scale(env.table$Salinity), 'euclidean')
Temp.dist = vegdist(scale(env.table$Temp), 'euclidean')
Depth.dist = vegdist(scale(env.table$Depth), 'euclidean')
Ca.dist = vegdist(scale(env.table$Ca), 'euclidean')
K.dist = vegdist(scale(env.table$K), 'euclidean')
Mg.dist = vegdist(scale(env.table$Mg), 'euclidean')
Na.dist = vegdist(scale(env.table$Na), 'euclidean')
MAT.dist = vegdist(scale(env.table$MAT), 'euclidean')
MAP.dist = vegdist(scale(env.table$MAP), 'euclidean')


ordination.distmat <- as.data.frame(cbind(bray.distmat, unifac.dismat, fun.dismat, DOC.dist,
                                          TC.dist, TN.dist,NH4_N.dist, NO3_N.dist, DIN.dist,
                                          DON.dist, S275_295.dist, SUVA254.dist, bix.dist,
                                          fi.dist, a300.dist, pH.dist, DO.dist, Conductivity.dist,
                                          Salinity.dist, Temp.dist, Depth.dist, Ca.dist,
                                          K.dist, Mg.dist, Na.dist, MAT.dist,MAP.dist))

long.ordination.distmat <- reshape::melt(ordination.distmat, id.vars = c('bray.distmat', 'unifac.dismat', 'fun.dismat'))
long.ordination.distmat1 <- reshape::melt(long.ordination.distmat, id.vars = c('variable'))

library(ggplot2)
library(ggpubr)
bray.s275_295.plot <- ggplot(ordination.distmat, aes(x = S275_295.dist,y = bray.distmat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="Taxonomic bray distance",x="ΔS275_295")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))


bray.suva254.plot <- ggplot(ordination.distmat, aes(x = SUVA254.dist,y = bray.distmat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="Taxonomic bray distance",x="ΔSUVA254")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

bray.ph.plot <- ggplot(ordination.distmat, aes(x = pH.dist,y = bray.distmat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="Taxonomic bray distance",x="ΔpH")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

bray.map.plot <- ggplot(ordination.distmat, aes(x = MAP.dist,y = bray.distmat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="Taxonomic bray distance",x="ΔMAP")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

bray.k.plot <- ggplot(ordination.distmat, aes(x =K.dist,y = bray.distmat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="Taxonomic bray distance",x="ΔK")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))



unifrac.s275_295.plot <- ggplot(ordination.distmat, aes(x = S275_295.dist,y = unifac.dismat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02, )+
  theme_bw()+
  labs(y="Phylogenetic unweighted unifrac",x="ΔS275_295")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

unifrac.suva254.plot <- ggplot(ordination.distmat, aes(x = SUVA254.dist,y = unifac.dismat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02, )+
  theme_bw()+
  labs(y="Phylogenetic unweighted unifrac",x="ΔSUVA254")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

unifrac.ph.plot <- ggplot(ordination.distmat, aes(x = pH.dist,y = unifac.dismat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02, )+
  theme_bw()+
  labs(y="Phylogenetic unweighted unifrac",x="ΔpH")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

unifrac.map.plot <- ggplot(ordination.distmat, aes(x = MAP.dist,y = unifac.dismat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02, )+
  theme_bw()+
  labs(y="Phylogenetic unweighted unifrac",x="ΔMAP")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

unifrac.k.plot <- ggplot(ordination.distmat, aes(x = K.dist,y = unifac.dismat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02, )+
  theme_bw()+
  labs(y="Phylogenetic unweighted unifrac",x="ΔK")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))


fun.s275_295.plot <- ggplot(ordination.distmat, aes(x = S275_295.dist,y = fun.dismat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="Functional bray distance",x="ΔS275_295")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

fun.suva254.plot <- ggplot(ordination.distmat, aes(x = SUVA254.dist,y = fun.dismat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="Functional bray distance",x="ΔSUVA254")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

fun.ph.plot <- ggplot(ordination.distmat, aes(x = pH.dist,y = fun.dismat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="Functional bray distance",x="ΔpH")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

fun.map.plot <- ggplot(ordination.distmat, aes(x = MAP.dist,y = fun.dismat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="Functional bray distance",x="ΔMAP")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

fun.k.plot <- ggplot(ordination.distmat, aes(x = K.dist,y = fun.dismat))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="Functional bray distance",x="ΔK")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

#
library(cowplot)
linner.distmat.between.env.com <- plot_grid(bray.s275_295.plot, bray.suva254.plot, bray.ph.plot,
                                            bray.map.plot, bray.k.plot, unifrac.s275_295.plot,
                                            unifrac.suva254.plot, unifrac.ph.plot, unifrac.map.plot,
                                            unifrac.k.plot, fun.s275_295.plot, fun.suva254.plot,
                                            fun.ph.plot, fun.map.plot, fun.k.plot,
                                            labels = 'auto', ncol = 5, nrow = 3, 
                  label_x = .01, label_y = 0.99, hjust = 0, label_size=14,align = "v")

ggsave('./figs/linner.distmat.between.env.com.pdf', linner.distmat.between.env.com, width = 15, height = 7, dpi = 300)




cor.test(DOC.dist, total.otu.dis)
cor.test(TOC.dist, total.otu.dis)
cor.test(TN.dist, total.otu.dis)
cor.test(NH4_N.dist, total.otu.dis)
cor.test(NO3_N.dist, total.otu.dis)
cor.test(DIN.dist, total.otu.dis)
cor.test(DON.dist, total.otu.dis)
cor.test(S275_295.dist, total.otu.dis)
cor.test(SUVA254.dist, total.otu.dis)
cor.test(a300.dist, total.otu.dis)
cor.test(bix.dist, total.otu.dis)
cor.test(fi.dist, total.otu.dis)
cor.test(hix.dist, total.otu.dis)
cor.test(pH.dist, total.otu.dis)
cor.test(DO.dist, total.otu.dis)
cor.test(Conductivity.dist, total.otu.dis)
cor.test(Salinity.dist, total.otu.dis)
cor.test(Temp.dist, total.otu.dis)
cor.test(Depth.dist, total.otu.dis)
cor.test(Ca.dist, total.otu.dis)
cor.test(K.dist, total.otu.dis)
cor.test(Mg.dist, total.otu.dis)
cor.test(Na.dist, total.otu.dis)
cor.test(MAT.dist, total.otu.dis)
cor.test(MAP.dist, total.otu.dis)

mantel(DOC.dist, total.otu.dis, "spearman", permutations = 999)
mantel(TC.dist, total.otu.dis, "spearman", permutations = 999)
mantel(TN.dist, total.otu.dis, "spearman", permutations = 999)
mantel(NH4_N.dist, total.otu.dis, "spearman", permutations = 999)
mantel(NO3_N.dist, total.otu.dis, "spearman", permutations = 999)
mantel(DIN.dist, total.otu.dis, "spearman", permutations = 999)
mantel(DON.dist, total.otu.dis, "spearman", permutations = 999)
mantel(S275_295.dist, total.otu.dis, "spearman", permutations = 999)
mantel(SUVA254.dist, total.otu.dis, "spearman", permutations = 999)
mantel(a300.dist, total.otu.dis, "spearman", permutations = 999)
mantel(bix.dist, total.otu.dis, "spearman", permutations = 999)
mantel(fi.dist, total.otu.dis, "spearman", permutations = 999)
mantel(hix.dist, total.otu.dis, "spearman", permutations = 999)
mantel(pH.dist, total.otu.dis, "spearman", permutations = 999)
mantel(DO.dist, total.otu.dis, "spearman", permutations = 999)
mantel(Conductivity.dist, total.otu.dis, "spearman", permutations = 999)
mantel(Salinity.dist, total.otu.dis, "spearman", permutations = 999)
mantel(Temp.dist, total.otu.dis, "spearman", permutations = 999)
mantel(Depth.dist, total.otu.dis, "spearman", permutations = 999)
mantel(Ca.dist, total.otu.dis, "spearman", permutations = 999)
mantel(K.dist, total.otu.dis, "spearman", permutations = 999)
mantel(Mg.dist, total.otu.dis, "spearman", permutations = 999)
mantel(Na.dist, total.otu.dis, "spearman", permutations = 999)

adonis(total.otu.dis~Na+K+Salinity+Conductivity+S275_295+SUVA254+Mg+hix+bix+TOC+Ca+pH+
         Depth+fi+Temp+NO3_N+NH4_N+DOC+a300+DIN+DON+TN+DO,data = env.table, permutations = 999)


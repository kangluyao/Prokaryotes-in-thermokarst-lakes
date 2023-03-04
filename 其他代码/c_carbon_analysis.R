library(reshape)
library(plyr)
# Identify microbial traits
# create object of trans_func
meco_dat <- phyloseq2meco(water_physeq)
t2 <- trans_func$new(meco_dat)

# mapping the taxonomy to the database
# the function can recognize prokaryotes or fungi automatically.
t2$cal_spe_func()
# return t2$res_spe_func, 1 represent function exists, 0 represent no or cannot confirmed.
t2$res_spe_func[1:5, 1:2]
# If you want to change the group list, reset the list t2$func_group_list
t2$func_group_list
# use show_prok_func to see the detailed information of prokaryotic traits
t2$show_prok_func("methanotrophy")
# calculate the percentages for communities
t2$cal_spe_func_perc(use_community = TRUE)
# t2$res_spe_func_perc[1:5, 1:2]
# then we try to correlate the res_spe_func_perc of communities to environmental variables
t3 <- trans_env$new(dataset = meco_dat, add_data = t2$sample_table[, c('pH', 'MAP', 'MAT','C4...', 'C3...', 'S275_295', 'SUVA254', 'DOC',
                                                                        'NH4_N', 'Na', 'K', 'Ca', 'Conductivity', 'Depth')])
t3$cal_cor(add_abund_table = t2$res_spe_func_perc[ ,colnames(t2$res_spe_func_perc) %in% 
                                                     c(t2$func_group_list$`C-cycle`)], cor_method = "spearman")
t3$plot_corr(pheatmap = TRUE)

#calculate the relative abundance of C cycle process
c_cycle_perc <- t2$res_spe_func_perc[ ,colnames(t2$res_spe_func_perc) %in% c(t2$func_group_list$`C-cycle`)]
group <- c(rep(1, 20), rep(2, 124), rep(3, 44))
group <- as.factor(group)
meta_data <- data.frame(group, c_cycle_perc)
melt_df <- melt(meta_data, id.vars = c('group'))
colnames(melt_df) <- c('group', 'cat_carbon','relative_abundance')
mode <- lmer(relative_abundance ~ cat_carbon + (1|group), data = melt_df)
anova(mode)
library(multcomp)
summary(glht(mode, linfct = mcp(cat_carbon = "Tukey")), test = adjusted("holm"))

new_df <- ddply(melt_df, c("cat_carbon"), summarise,
                mean = mean(relative_abundance), sd = sd(relative_abundance),
                se = sd(relative_abundance)/sqrt(length(relative_abundance)))
new_df <- data.frame(new_df, sig = c('bc', 'a', 'd', 'ab', 'c', 'e'))
ggplot(new_df,aes(x = cat_carbon, y = mean, fill = cat_carbon))+
  geom_bar(stat = 'identity', width = 0.7)+
  scale_fill_manual(values=c("#F2B379","#479E9B","#FDDC7B","#4169B2","#ACCDDC","#DD5F60","#B1A4C0"))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se), width=.2)+
  scale_x_discrete(limits = c("chitinolysis", "cellulolysis", "fermentation", 
                              "methanogenesis", "methanotrophy", "methylotrophy"))+
  labs(x = 'Carbon cycle', y = 'Mean relative abundance (%)',
       fill = "Carbon cycle")+
  scale_y_continuous(expand = c(0,0), limits = c(0,5))+
  geom_text(aes(label = sig, y = mean + se+0.02*max(mean)), position = position_dodge(0.9), vjust = 0)+
  theme_bw()+
  theme(legend.position = c(0.85,0.8),
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

#extract asvs involved in carbon cycling
c_cycle <- t2$res_spe_func[ ,c(t2$func_group_list$`C-cycle`)]
c_cycle_asvs <- rownames(c_cycle[rowSums(c_cycle) != 0,])

carbon_phy <- subset_taxa(water_physeq, ASV %in% c_cycle_asvs)
carbon_phy_rel <- subset_taxa(water.rel, ASV %in% c_cycle_asvs)
env.table <- read.csv(file = './edge_analysis/env.csv', header = T, row.names = 1, stringsAsFactors = F)
sample_data(carbon_phy_rel) <- env.table

#phylogenetic tree
library(ggtreeExtra)
library(dplyr)
ps6 <- tax_glom(carbon_phy, taxrank="Genus")
taxa_sums(ps6)
melt_simple <- data.frame(Phylum = data.frame(tax_table(ps6))[2],
                          Genus = data.frame(tax_table(ps6))[6],
                          ASV = taxa_names(ps6),
                          abundance = taxa_sums(ps6))
f<-phy_tree(ps6)
tax = as.data.frame(tax_table(ps6))
# change tip labels to genera / need to gsub
library(ggtree)
mytree = phy_tree(f)$tip.label
foo = data.frame(label=mytree, label2=tax$Genus, label3=tax$Phylum, stringsAsFactors = F)
#foo$label2[foo$label2 == NA] <- 'uncultured'
# Make a less putrid color palette
library(RColorBrewer)
colourCount = length(unique(tax$Phylum))
Mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(colourCount)
names(Mycolors) <- levels(as.factor(tax$Phylum))
colScale <- scale_colour_manual(name = "grp",values = Mycolors)

p<-ggtree(f, branch.length='none', layout='circular', open.angle = 20) %<+% foo + 
  geom_tiplab(size=2, align=TRUE, linesize=0.5, aes(label=label2), offset=5.5, hjust =0) +
  geom_tippoint(size=2, aes(label=label2, group=label3, color=label3) ) + 
  theme_tree2() + 
  scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(colourCount))

p_tree <- p + geom_fruit(
  data = melt_simple,
  geom = geom_bar,
  mapping = aes(x = abundance, y = ASV, fill = 'red4'),
  orientation="y",
  stat="identity",
  size=.2,
  outlier.size=0.5,
  outlier.stroke=0.08,
  outlier.shape=21,
  axis.params=list(
    axis       = "x",
    text.size  = 1.8,
    hjust      = 1,
    vjust      = 0.5,
    nbreak     = 3,
  ),
  grid.params=list()
  ) 
p_tree
#comunity composition
##determine the genus compositions within each family##
ps6 <- tax_glom(carbon_phy_rel, taxrank="Genus")

genus.ra.table <- otu_table(ps6)
MRA <- rowMeans(genus.ra.table)
group <- tax_table(ps6)[,c(5,6)]
genus.mra.table <- data.frame(group,MRA)

#arrange the genuss table
library(tidyr)
genus.mra.table <- genus.mra.table %>% spread(Genus, MRA)
genus.mra.table[is.na(genus.mra.table)] <- 0
rownames(genus.mra.table)<-genus.mra.table$Family
genus.mra.table<-as.matrix(t(genus.mra.table[,-1])*100)
colsum <-apply(genus.mra.table,2,sum)
rowsum<-apply(genus.mra.table,1,sum)
topgenus_table<-(genus.mra.table[order(rowsum,decreasing=TRUE),order(colsum,decreasing=TRUE)])[,1:10]
head(topgenus_table)
topgenus_table<-as.matrix(topgenus_table)
# Get the stacked barplot
# create color palette:
#library(RColorBrewer)
#coul <- brewer.pal(nrow(top10_phy), "Pastel2") 
mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,
          558,43,652,31,610,477,588,99,81,503,562,76,96,495,77,12,90,
          345,255,401,366,276,158,436)
mycol <-colors()[rep(mycol,nrow(topgenus_table))]

#tiff(file="bar.order1.level.tiff",width=750,height=700,pointsize=15)
layout(matrix(1:2,2,1),heights=c(1.5:1))
opar <- par(no.readonly = T)
par(mar=c(2,5,2,2))
barplot(topgenus_table, width = 1.8, space = 0.4, plot = T,las = 2,
        col = mycol[1:nrow(topgenus_table)], cex.axis = 1, cex.names = 1, border = NA,
        xlab = 'Family',ylab = "Mean relative abundance (%)",
        offset = 0, cex.lab = 1.2)
par(mar=c(2,3.5,3.5,1))
plot.new()
legend("topleft",legend=rownames(topgenus_table),
       ncol=3,fill=mycol[1:nrow(topgenus_table)],cex=0.8,bty="n")
par(opar)

#edgeR analysis
# Glom OTUs to genus level for further statistical analysis & reasonable power
ps6 <- tax_glom(carbon_phy_rel, "Genus", NArm = TRUE)
#Use code snippet straight provided by authors of PhyloSeq to export phyloseq object to an EdgeR object
phyloseq_to_edgeR = function(physeq, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

# Make normalized phyloseq object (ps6) into an edgeR object. It needs a grouping factor. We use location.
dge = phyloseq_to_edgeR(ps6)
#The crunching to follow is much easier if metadata is pulled out this way into object "a"
a = sample_data(ps6)

#So many things don't understand what kind of variables they need to be. Make sure they understand
chitinolysis <- as.numeric(a$chitinolysis)
cellulolysis <- as.numeric(a$cellulolysis)
fermentation <- as.numeric(a$fermentation)
methanogenesis <- as.numeric(a$methanogenesis)
methanotrophy <- as.numeric(a$methanotrophy)
methylotrophy <- as.numeric(a$methylotrophy)
DOC <- as.numeric(a$DOC)
TN <- as.numeric(a$TN)
NH4_N <- as.numeric(a$NH4_N) 
NO3_N <- as.numeric(a$NO3_N)
DON <- as.numeric(a$DON)
S275_295 <- as.numeric(a$S275_295)
SUVA254 <- as.numeric(a$SUVA254) 
conductivity <- as.numeric(a$Conductivity)
a300 <- as.numeric(a$a300)
fi <- as.numeric(a$fi)
bix <- as.numeric(a$bix)
hix <- as.numeric(a$hix)
Comp1 <- as.numeric(a$C1...)
Comp2 <- as.numeric(a$C2...)
Comp3 <- as.numeric(a$C3...)
Comp4 <- as.numeric(a$C4...)
pH <- as.numeric(a$pH)
DO <- as.numeric(a$DO)
Salinity <- as.numeric(a$Salinity) 
Depth <- as.numeric(a$Depth)
MAT <- as.numeric(a$MAT)
MAP <- as.numeric(a$MAP)
Conductivity <- as.numeric(a$Conductivity) # transforming this fixes it so it doesn't kill linear model
K <- as.numeric(a$K)
Ca <- as.numeric(a$Ca)
Na <- as.numeric(a$Na)
Mg <- as.numeric(a$Mg)
##Test and forward selection of environmental variables
Env_select <- function (Dist_Matrix,Env,Number_Permutations=999) {
  mod1<-capscale(Dist_Matrix~.,Env,add=T,na.action = na.omit)
  mod0<-capscale(Dist_Matrix~1,Env,add=T,na.action = na.omit)
  mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999,trace = F)
  Env.P<-anova.cca(mod,permutations=Number_Permutations)[[4]][1]
  Env_se<-mod$CCA$biplot
  Env_se<-Env[,rownames(Env_se)]
  return(Env_se)
}
carbon_asvs <- otu_table(carbon_phy_rel)
carbon.dist <- vegdist(sqrt(t(carbon_asvs)))
env.table <- sample_data(carbon_phy_rel)
env_vars <- data.frame(scale(env.table[ ,-c(1:6)]))
env_se_carbon <- Env_select (carbon.dist, env_vars)
colnames(env_se_carbon)
# Design for my linear model
design <-model.matrix(~ chitinolysis + cellulolysis + fermentation + methanogenesis + 
                        methanotrophy + methylotrophy + pH + MAP + MAT +
                        Comp4 + Comp3 + S275_295 + SUVA254 + DOC + NH4_N +
                        Na + K + Ca + Conductivity + Depth) 

# EdgeR needs to calculate dispersion again after you've fed it the design.
#This step projectile vomits if chloride or conductivity aren't sqrt transformed.
x = calcNormFactors(dge, method="RLE")
x = estimateGLMCommonDisp(dge, design)
x = estimateGLMTrendedDisp(dge, design)
x = estimateGLMTagwiseDisp(dge, design)

fit <-glmFit(x, design)

# grab the coefficients I care about
lrt <- glmLRT(fit, coef=2:21)

# lrt to z-scores
table<-lrt$table
table<-apply(table, 2, function(x) scale(x, center = TRUE, scale = TRUE))
rownames(table) = rownames(lrt$table)
table<-table[,1:20]
table<-as.data.frame(table)

q<-lrt$genes
table1<-cbind(table, Kingdom = q$Kingdom, Phylum = q$Phylum, Class = q$Class, Order = q$Order,
             Family = q$Family, Genus = q$Genus, OTU = rownames(table))

# Sometimes want to do this to check what has at least one z score beyond threshold
table2 <- table1[apply(table1[,1:20], 1, function(x) any(abs(x)>1.96)), ]

#write.csv(table2, "./tables/EdgeR-Zscores-Genera.csv")
melted<-melt(table2)
melt2<-subset(melted, abs(melted$value) > 1.96 )

filter<-unique(melt2$OTU)
mini<-prune_taxa((rownames(otu_table(water_physeq)) %in% filter), water_physeq)
f<-phy_tree(mini)
tax = as.data.frame(tax_table(mini))

tax$Genus = gsub("D_5__", "", tax$Genus)
tax$Phylum = gsub("D_1__", "", tax$Phylum)

# change tip labels to genera / need to gsub
library(ggtree)
mytree = phy_tree(f)$tip.label
foo = data.frame(label=mytree, label2=paste(tax$Genus, tax$ASV, sep='_'), label3=tax$Phylum, stringsAsFactors = F)
#foo$label2[foo$label2 == NA] <- 'uncultured'
# Make a less putrid color palette

library(RColorBrewer)

colourCount = length(unique(tax$Phylum))
Mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(colourCount)
names(Mycolors) <- levels(as.factor(tax$Phylum))
colScale <- scale_colour_manual(name = "grp",values = Mycolors)

p<-ggtree(f) %<+% foo + 
  geom_tiplab(size=2, align=TRUE, linesize=0.5, aes(label=label2), hjust =0) +
  geom_tippoint(size=2, aes(label=label2, group=label3, color=label3) ) + 
  theme_tree2() + 
  scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(colourCount))
p + xlim (0, 3) 

test<-as.matrix(table2[1:20])
tested<-apply(test, 2, function(x){cut(x, br=c(-14, -4, -2, 2, 4, 14))})
rownames(tested) = rownames(table2)
colnames(tested) = gsub("logFC.", "", colnames(tested))
colnames(tested) = gsub("logCPM", "Rotorua", colnames(tested))

heatmapcols <-colorRampPalette(brewer.pal(6, "Spectral"))(6)
names(heatmapcols) <- levels(as.factor(tested[1:10000]))
heatmapcols <- c(heatmapcols[5],heatmapcols[4],heatmapcols[2],heatmapcols[3],heatmapcols[1],heatmapcols[6])
# Make sure to check that the order of these is correct if this is rerunhe
heatmapcols[1] = "#881111"
heatmapcols[2] = "#FFCCCF"
heatmapcols[3] = "#F2F2F2"
heatmapcols[4] = "#AABBDD"
heatmapcols[5] = "#112288"
heatmapcols[6] = "#112288"

mycols <- c("#881111","#D9444D","#FFCCCF","#F2F2F2","#AABBDD","#112288")
tested <- tested[,c('chitinolysis', 'cellulolysis', 'fermentation', 'methanogenesis', 
                      'methanotrophy', 'methylotrophy',  'pH', 'MAP', 'MAT',
                      'Comp4', 'Comp3', 'S275_295', 'SUVA254', 'DOC', 'NH4_N',
                      'Na', 'K', 'Ca', 'Conductivity', 'Depth')] 
p1<-gheatmap(p, tested, offset = 0.25, width = 3, font.size=2.5, 
             colnames_angle=-90, colnames_position = "top", hjust=1) + 
  scale_fill_manual(values=heatmapcols)
p1 <-p1 + theme(legend.position="right")
p1
#ggsave(filename="output/edgeR-Genus.pdf", plot=p, width=8, height=10, units="in")


#edgeR for functional gene
meanfun <- colSums(fun.table)/nrow(fun.table)
fun.order.table<-(fun.table[ ,order(meanfun,decreasing=TRUE)])
fun.order.table[1:5,1:5]
class.ra.table <- otu_table(class.phy.rel)
MRA <- rowSums(class.ra.table)/ncol(class.ra.table)
group <- tax_table(class.phy.rel)[,c(2,6)]
class.mra.table <- data.frame(group,MRA)

#arrange the class table
library(tidyr)
class.mra.table$Genus <- paste(class.mra.table$Genus, rownames(class.mra.table), sep = '_')
class.mra.table <- class.mra.table %>% spread(Genus, MRA)
class.mra.table[is.na(class.mra.table)] <- 0
rownames(class.mra.table)<-class.mra.table$Phylum
class.mra.table<-as.matrix(t(class.mra.table[,-1])*100)
colsum <-apply(class.mra.table,2,sum)
rowsum<-apply(class.mra.table,1,sum)
top10phy_table<-(class.mra.table[order(rowsum,decreasing=TRUE),order(colsum,decreasing=TRUE)])[,1:10]
head(top10phy_table)
top10phy_table<-as.matrix(top10phy_table)

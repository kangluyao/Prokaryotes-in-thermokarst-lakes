setwd('E:/thermokast_lakes/water_microbes/tibet_dada2_asv/data/methano_data/')
#load libraries
library(ggplot2)
library(phyloseq)
library(ape)
library(Biostrings)

tax.table <- read.csv('./taxonomy.csv',sep=",",row.names=1)
methy_asvs <- tax.table$ASV
methy_phy <- subset_taxa(water.rel, ASV %in% methy_asvs)
otu.table <- otu_table(methy_phy)
ref.seqs <- refseq(methy_phy)
tree <- phy_tree(methy_phy)
metadata <- read.csv("E:/thermokast_lakes/water_microbes/tibet_dada2_asv/data/metadata.csv", row.names=1, header = T)

tax.table <- tax_table(as.matrix(tax.table))
metadata <- sample_data(metadata)
methano_phylo <- phyloseq(tax.table, otu.table, metadata, tree, ref.seqs)

#comunity composition
##determine the genus compositions within each family##
genus.phy <- tax_glom(methano_phylo, taxrank="Genus")

genus.ra.table <- otu_table(genus.phy)
MRA <- rowMeans(genus.ra.table)
group <- tax_table(genus.phy)[,c(5,6)]
genus.mra.table <- data.frame(group,MRA)

#arrange the genuss table
library(tidyr)
genus.mra.table <- genus.mra.table %>% spread(Genus, MRA)
genus.mra.table[is.na(genus.mra.table)] <- 0
rownames(genus.mra.table)<-genus.mra.table$Family
genus.mra.table<-as.matrix(t(genus.mra.table[,-1])*100)
colsum <-apply(genus.mra.table,2,sum)
rowsum<-apply(genus.mra.table,1,sum)
topgenus_table<-(genus.mra.table[order(rowsum,decreasing=TRUE),order(colsum,decreasing=TRUE)])
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

#relative abundance analysis
methaon_rel_abun <- colSums(genus.ra.table)
type1_methano_phylo <- subset_taxa(methano_phylo, Type == 'Type1') 
type1_methano_rel_abun <- colSums(otu_table(type1_methano_phylo))
type2_methano_phylo <- subset_taxa(methano_phylo, Type == 'Type2')
type2_methano_rel_abun <- colSums(otu_table(type2_methano_phylo))

methano_abundance_table <- cbind(methaon_rel_abun, type1_methano_rel_abun, type2_methano_rel_abun)
env.table <- sample_data(methano_phylo)
cor.methano.env <- cor.mat.cal(methano_abundance_table, env.table[,-(1:5)], 'pearson')

#####corplot  
ggplot(cor.methano.env, aes(y = Taxa, x = Env, colour = Taxa)) +
  geom_point(aes(size = abs(Correlation)))+ #,shape=1
  scale_colour_manual(values=c("#B1A4C0","#479E9B","#FDDC7B","#4169B2",
                               "#ACCDDC","#DD5F60","#F2B379","#7CB6DB",
                               "pink2","skyblue","Gold2","#009E73",
                               "#E69F00","dodgerblue","#E69F00"))+
  theme_bw(base_size = 12)+
  geom_text(aes(label=Significance), color="black", size = 5.5)+
  scale_size(range = c(1, 12),breaks=c(0.2,0.4,0.6,0.8), limits=c(0,0.8))+ 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  #theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = 'Environmental factors', y = 'Phylum',
       colour = "Phylum", size = "Pearson's r")+ #legend title
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


library(lme4)
library(lmerTest)
library(MuMIn)
methan_env_table <- cbind(methano_abundance_table, env.table)
vars <- colnames(env.table)[-(1:5)]
model1 <- lapply(vars, function(x) {
  lmer(substitute(sqrt(methaon_rel_abun) ~ i + (1|Site),list(i = as.name(x))), data = methan_env_table)})
lapply(model1, summary)

model2 <- lapply(vars, function(x) {
  lmer(substitute(sqrt(type1_methano_rel_abun) ~ i + (1|Site),list(i = as.name(x))), data = methan_env_table)})
lapply(model2, summary)

model3 <- lapply(vars, function(x) {
  lmer(substitute(sqrt(type2_methano_rel_abun) ~ i + (1|Site),list(i = as.name(x))), data = methan_env_table)})
lapply(model3, summary)

taxa_sums(Methanotrophy_table)

methan_MRA <- data.frame(OTU = names(methan_MRA), MRA = as.numeric(methan_MRA))

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

# Make normalized phyloseq object (ps2) into an edgeR object. It needs a grouping factor. We use location.
mge = phyloseq_to_edgeR(methano_phylo)
#The crunching to follow is much easier if metadata is pulled out this way into object "a"
a = sample_data(methano_phylo)

#So many things don't understand what kind of variables they need to be. Make sure they understand
DOC <- as.numeric(a$DOC)
#TC <- as.numeric(a$TC)
TN <- as.numeric(a$TN)
NH4_N <- as.numeric(a$NH4_N) 
NO3_N <- as.numeric(a$NO3_N)
#DIN <- as.numeric(a$DIN)
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
#Temp <- as.numeric(a$Temp)
Depth <- as.numeric(a$Depth)
MAT <- as.numeric(a$MAT)
MAP <- as.numeric(a$MAP)
Conductivity <- as.numeric(a$Conductivity) # transforming this fixes it so it doesn't kill linear model
K <- as.numeric(a$K)
Ca <- as.numeric(a$Ca)
Na <- as.numeric(a$Na)
Mg <- as.numeric(a$Mg)
# Design for my linear model
options(na.action='na.pass')
design <-model.matrix(~ DOC + TN + NH4_N + NO3_N + fi + bix + hix + Comp1 + Comp2 + Comp3 + Comp4 + pH +
                        MAP + S275_295 + SUVA254 + pH + K + Conductivity + Depth + 
                        Salinity + Na + DO + MAT + a300 + Ca + Mg) 

# EdgeR needs to calculate dispersion again after you've fed it the design. Why? I don't know.
#This step projectile vomits if chloride or conductivity aren't sqrt transformed.
x = calcNormFactors(mge, method="RLE")
x = estimateGLMCommonDisp(mge, design)
x = estimateGLMTrendedDisp(mge, design)
x = estimateGLMTagwiseDisp(mge, design)

fit <-glmFit(x, design)

# grab the coefficients I care about
lrt <- glmLRT(fit, coef=2:24)

# lrt to z-scores
table<-lrt$table
table<-apply(table, 2, function(x) scale(x, center = TRUE, scale = TRUE))
rownames(table) = rownames(lrt$table)
table<-table[,1:23]
table<-as.data.frame(table)

q<-lrt$genes
table<-cbind(table, Kingdom = q$Kingdom, Phylum = q$Phylum, Class = q$Class, Order = q$Order,
             Family = q$Family, Genus = q$Genus, OTU = rownames(table))

# Sometimes want to do this to check what has at least one z score beyond threshold
#table2<-table[apply(table[1:15], 1, function(x) any(abs(x)>1.96)), ]

#write.csv(table2, "./tables/EdgeR.methaon-Zscores-Genera.csv")
melted<-melt(table)
#melt2<-subset(melted, abs(melted$value) > 1.96 )

filter<-unique(melted$OTU)
mini<-prune_taxa((rownames(otu_table(water_physeq)) %in% filter), water_physeq)
f<-phy_tree(mini)
tax = as.data.frame(tax_table(mini))

tax$Genus = gsub("D_5__", "", tax$Genus)
tax$Phylum = gsub("D_1__", "", tax$Phylum)

# change tip labels to genera / need to gsub
library(ggtree)
mytree = phy_tree(f)$tip.label
foo = data.frame(label=mytree, label2= paste(tax$Genus, tax$OTU, sep = '_'), label3=tax$Family, stringsAsFactors = F)
#foo$label2[foo$label2 == NA] <- 'uncultured'
# Make a less putrid color palette

library(RColorBrewer)

colourCount = length(unique(tax$Family))
Mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(colourCount)
names(Mycolors) <- levels(as.factor(tax$Phylum))
colScale <- scale_colour_manual(name = "grp",values = Mycolors)

p<-ggtree(f) %<+% foo + 
  geom_tiplab(size=2, align=TRUE, linesize=0.5, aes(label=label2), hjust =0) +
  geom_tippoint(size=2, aes(label=label2, group=label3, color=label3) ) + 
  #geom_hilight(node=34, fill="pink", alpha=0.5) +
  #geom_hilight(node=25, fill="yellow", alpha=0.5) +
  theme_tree2() + 
  scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(colourCount))
p + xlim (0, 3) 

test<-as.matrix(table[1:23])
tested<-apply(test, 2, function(x){cut(x, br=c(-8, -6, -4, -1.96, 1.96, 4, 6, 8))})
rownames(tested) = rownames(table)
colnames(tested) = gsub("logFC.", "", colnames(tested))
colnames(tested) = gsub("logCPM", "Rotorua", colnames(tested))

heatmapcols <-colorRampPalette(brewer.pal(7, "Spectral"))(7)
names(heatmapcols) <- levels(as.factor(tested[1:10000]))

# Make sure to check that the order of these is correct if this is rerunhe
heatmapcols[1] = "#F2F2F2" # -1.96 to 1.96
heatmapcols[2] = "#AABBDD" #"#D7E3F4" # -1.96 to -4
heatmapcols[3] = "#112288" # -4 to -6 - dk blue
heatmapcols[4] = "#FFCCCF"
heatmapcols[5] = "#D9444D"
heatmapcols[6] = "#881111" 
heatmapcols[7] = 
  
p1<-gheatmap(p, tested, offset = 0.5, width = 4.5, font.size=2.5, 
             colnames_angle=-90, colnames_position = "top", hjust=1)+
  scale_fill_manual(values=heatmapcols)+
  theme(legend.position="right")
p1


















env.table <- as(sample_data(water_physeq), 'data.frame')
sub.envtable <- cbind(env.table[,1:3], 
                      diver_table[,!colnames(diver_table) %in% c('se.chao1', 'se.ACE', 'Observed', 'ACE', 'InvSimpson','Fisher','SR')])
geo.diver <- aggregate(sub.envtable[,-1], by = list(sub.envtable$Site), mean)

dat <- cbind(env.table, sub.envtable)
meadow.dat <- subset(dat, vegetable_type == 'Meadow')
plot(meadow.dat$MAP, meadow.dat$Chao1)
library(ggplot2)
library(ggpubr)
ggplot(meadow.dat, aes(x = MAP, y = Shannon))+
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


group<-data.frame(c(rep('all',188)))
x1 <- data.frame(t(otu.table))
rownames(group)<-rownames(x1)
library(NST)
tnst.bacteria=tNST(comm=x1, group=group, dist.method="jaccard",
                   abundance.weighted=TRUE, rand=100,output.rand = T,
                   nworker=4, null.model="PF", between.group=F,
                   SES=T, RC=T)

bacteria.nst.bt=nst.boot(nst.result=tnst.bacteria, group=NULL, rand=99,
                         trace=TRUE, two.tail=FALSE, out.detail=T,
                         between.group=FALSE, nworker=1)
NST.bacteria<-bacteria.nst.bt$detail$NST.boot$all
tnst.bacteria$index.pair

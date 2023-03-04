library(phyloseq)
library(biomformat)
library(edgeR)
library(reshape2)
library(ggplot2)
library(ggtree)
library(microbiome)

#prune out samples that I don't want to keep.
x1.genus <- tax_glom(water.rel, taxrank="Genus")
pruned<- microbiome::core(x1.genus, detection = 0.1/100, prevalence = 5/100)
pruned <- prune_taxa(taxa_names(pruned), water_physeq)


# Glom OTUs to genus level for further statistical analysis & reasonable power
ps2 = tax_glom(pruned, "Genus", NArm = TRUE)

# Just for fun - how much variance do we have in the original dataset vs after the OTUs are glommed to genus level?

hist(log10(apply(otu_table(water_physeq), 1, var)),
     xlab="log10(variance)", breaks=50,
     main="Variance within the dataset")

hist(log10(apply(otu_table(ps2), 1, var)),
     xlab="log10(variance)", breaks=50,
     main="Variance within the dataset")

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
dge = phyloseq_to_edgeR(ps2)
#The crunching to follow is much easier if metadata is pulled out this way into object "a"
a = sample_data(ps2)

#So many things don't understand what kind of variables they need to be. Make sure they understand
DOC <- as.numeric(a$DOC)
TC <- as.numeric(a$TC)
TN <- as.numeric(a$TN)
NH4_N <- as.numeric(a$NH4_N) 
NO3_N <- as.numeric(a$NO3_N)
DIN <- as.numeric(a$DIN)
DON <- as.numeric(a$DON)
S275_295 <- as.numeric(a$S275_295)
SUVA254 <- as.numeric(a$SUVA254) 
conductivity <- as.numeric(a$Conductivity)
a300 <- as.numeric(a$a300)
fi <- as.numeric(a$fi)
bix <- as.numeric(a$bix)
hix <- as.numeric(a$hix)
pH <- as.numeric(a$pH)
DO <- as.numeric(a$DO)
Salinity <- as.numeric(a$Salinity) 
Temp <- as.numeric(a$Temp)
Depth <- as.numeric(a$Depth)
MAT <- as.numeric(a$MAT)
MAP <- as.numeric(a$MAP)
Conductivity <- as.numeric(a$Conductivity) # transforming this fixes it so it doesn't kill linear model
K <- as.numeric(a$K)
Ca <- as.numeric(a$Ca)
Na <- as.numeric(a$Na)
Mg <- as.numeric(a$Mg)
# Design for my linear model
design <-model.matrix(~ MAP + S275_295 + SUVA254 + pH + K + Conductivity + TC + Depth
                      +TN + Na + Temp + DON + DO + MAT + a300) 

# EdgeR needs to calculate dispersion again after you've fed it the design. Why? I don't know.
#This step projectile vomits if chloride or conductivity aren't sqrt transformed.
x = calcNormFactors(dge, method="RLE")
x = estimateGLMCommonDisp(dge, design)
x = estimateGLMTrendedDisp(dge, design)
x = estimateGLMTagwiseDisp(dge, design)

fit <-glmFit(x, design)

# grab the coefficients I care about
lrt <- glmLRT(fit, coef=2:16)

# lrt to z-scores
table<-lrt$table
table<-apply(table, 2, function(x) scale(x, center = TRUE, scale = TRUE))
rownames(table) = rownames(lrt$table)
table<-table[,1:15]
table<-as.data.frame(table)

q<-lrt$genes
table<-cbind(table, Kingdom = q$Kingdom, Phylum = q$Phylum, Class = q$Class, Order = q$Order,
             Family = q$Family, Genus = q$Genus, OTU = rownames(table))

# Sometimes want to do this to check what has at least one z score beyond threshold
table2<-table[apply(table[1:15], 1, function(x) any(abs(x)>1.96)), ]

write.csv(table2, "./tables/EdgeR-Zscores-Genera.csv")
melted<-melt(table2)
melt2<-subset(melted, abs(melted$value) > 1.96 )

filter<-unique(melt2$OTU)
mini<-prune_taxa((rownames(otu_table(core.phy)) %in% filter), core.phy)
f<-phy_tree(mini)
tax = as.data.frame(tax_table(mini))

tax$Genus = gsub("D_5__", "", tax$Genus)
tax$Phylum = gsub("D_1__", "", tax$Phylum)

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

p<-ggtree(f) %<+% foo + 
  geom_tiplab(size=2, align=TRUE, linesize=0.5, aes(label=label2), hjust =0) +
  geom_tippoint(size=2, aes(label=label2, group=label3, color=label3) ) + 
  theme_tree2() + 
  scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(colourCount))
p + xlim (0, 3) 

test<-as.matrix(table2[1:15])
tested<-apply(test, 2, function(x){cut(x, br=c(-8, -6, -4, -1.96, 1.96, 4, 6, 8))})
rownames(tested) = rownames(table2)
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


p1<-gheatmap(p, tested, offset = 1.5, width = 3.5, font.size=2.5, 
             colnames_angle=-90, colnames_position = "top", hjust=1) + 
  scale_fill_manual(values=heatmapcols)
p1 <-p1 + theme(legend.position="right")
p1
ggsave(filename="output/edgeR-Genus.pdf", plot=p, width=8, height=10, units="in")


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

library(GUniFrac)
library(phangorn)
library(vegan)
library(picante)
library(phyloseq)
library(ggplot2)
library(microbiome)

#read in otu table
otu.table <- read.table("./module_data/ito/module.anno.2.tsv", header = T, row.names = 1)
otu.table <- as.matrix(t(otu.table))

#read in taxonomy
#seperated by kingdom phylum class order family genus species 
tax.table <- read.table("./module_data/ito/module.anno.1.tsv", header = T, row.names = 1)
tax.table <- as.matrix(tax.table)

#read in metadata
meta.table <- read.csv("./new_data/metadata.csv", row.names=1, header = T)

# read in tree
phy_tree <- read_tree("./module_data/ito/module.ito.tree.tre")

#read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "./new_data/ref.seqs.fasta",
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#import as phyloseq objects
otu.table <- otu_table(otu.table, taxa_are_rows = TRUE)
tax.table <- tax_table(tax.table)
meta.table <- sample_data(meta.table)
#(tree was already imported as a phyloseq object)

#merge into one phyloseq object
module.phy <- phyloseq(otu.table, tax.table, meta.table, phy_tree, ref_seqs)

#prune the phyloseq object after discard the OTUs less than 10
tips.name <- NULL
for (i in module.num) {
  tips <- wtc[[i]]
  if(is.null(tips.name)){
    tips.name <- tips
  } else {
    tips.name = c(tips.name,tips)
  }
}

module.13 <- prune_taxa(tips.name, module.phy)
module.13
module.rel <- microbiome::transform(module.13, "compositional")
tax_table(module.rel)

#Make an PCoA ordination plot based on abundance using unifrac distances
ord_rel <- ordinate(module.rel, method="PCoA", distance="unifrac", weighted = TRUE)
ord_rel_un <- ordinate(module.rel, method="PCoA", distance="unifrac", weighted = F)

#plot
mycol <-c(119,132,147,454,89,404,329,123,463,104,552,28,54)
mycol <-colors()[mycol]

myshape <-c(119,132,147,454,89,404,329,123,463,104,552,28,54)

p <- plot_ordination(module.rel, ord_rel_un, type="taxa", shape = 'Module', color = 'Module', 
                     title="Phyloseq's unWeighted Unifrac")
  #stat_ellipse(aes(group = Module), level=0.95, lwd=1, linetype = 2)
p <- p + geom_point(size=1) +
  scale_color_manual(values = mycol)+
  scale_shape_manual(values=c(1,rep(c(0:2,5:6,9:10,11:12,14), times=6)))+
  theme_bw()
print(p)


















#export the files for ito object
module.ito.13.tree<-phy_tree(module.13)
write.tree(module.ito.13.tree, file = './module_data/ito13/module.ito.13.tree.tre', append = FALSE,
           digits = 10, tree.names = FALSE)
module.13.anno.1<-data.frame(as(tax_table(module.13), "matrix"))
module.13.anno.1<-data.frame(Tip=rownames(module.13.anno.1),module.13.anno.1)

write.table(module.13.anno.1, file = "./module_data/ito13/module.13.anno.1.tsv", row.names=FALSE, sep="\t")


module.13.abund.table<-data.frame(t(otu_table(module.13)))
module.13.abund.table<-data.frame(Tip=rownames(module.13.abund.table),module.13.abund.table)
write.table(module.13.abund.table, file = "./module_data/ito13/module.13.anno.2.tsv", row.names=FALSE, sep="\t")

ref.13.seqs <- refseq(module.13)
writeXStringSet(ref.13.seqs, './module_data/tax4fun2/ref.13.seqs.fasta', append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

df <- NULL
for (i in module.num) {
  tips <- wtc[[i]]
  wp_modul <- prune_taxa(tips, water_physeq)
  abundance<-taxa_sums(wp_modul)
  tmp<-data.frame(OTU_ID = tips, abundance)
  names(tmp)[2] <- paste('Module', i, sep = '')
  if(is.null(df)){
    df <- tmp
  }
  else{
    df <- merge(df, tmp, by = 'OTU_ID', all = TRUE)
  }
  abun.table <- df
}

abun.table[is.na(abun.table)] <- 0
write.table(abun.table, file = "./module_data/tax4fun2/otu_table.txt", row.names=FALSE, sep="\t")


source('C:/Program Files/R/R-3.6.0/library/table2itol/table2itol.R')

#read.table(file = 'drug_info.tsv', sep = '\t', header = TRUE)

create_itol_files(infiles = c("./module_data/ito13/module.13.anno.1.tsv",
                              './module_data/ito13/module.13.anno.2.tsv'),
                  identifier = "Tip", label = "OTU", na.strings = "X")


#The GUniFrac package can also be used to calculate unifrac distances and has additional features. 
#Unifrac distances are traditionally calculated on either presence/absence data, or abundance data. 
#The former can be affected by PCR and sequencing errors leading to a high number of spurious and 
#usually rare OTUs, and the latter can give undue weight to the more abundant OTUs. 
#GUniFrac's methods include use of a parameter alpha that controls the weight given to abundant OTUs 
#and also a means of adjusting variances.


#The function GUniFrac requires a rooted tree, but unlike phyloseq's ordination function 
#will not try to root an unrooted one. We will apply mid-point rooting with the midpoint function 
#from the phangorn package

unifracs <- GUniFrac(t(as.matrix(otu_table(module.rel))), 
                     midpoint(phy_tree(module.rel)), alpha = c(0, 0.5, 1))$unifracs
# We can extract a variety of distance matrices with different weightings.
dw <- unifracs[, , "d_1"]  # Weighted UniFrac
du <- unifracs[, , "d_UW"]  # Unweighted UniFrac

# use vegan's cmdscale function to make a PCoA ordination from a distance matrix.
pcoa_dw <- cmdscale(dw, k = nrow(t(as.matrix(otu_table(module.rel)))) - 1, eig = TRUE, add = TRUE)
pcoa_du <- cmdscale(du, k = nrow(t(as.matrix(otu_table(module.rel)))) - 1, eig = TRUE, add = TRUE)

env.table <- data.frame(sample_data(water_physeq))
fit.phy.core <- envfit(pcoa_du, env.table, perm=999, na.rm = TRUE)


p_dw<-plot_ordination(module.rel, pcoa_dw, type="taxa", color = 'Module', shape = 'Module',
                   title="GUniFrac Weighted Unifrac") +
  stat_ellipse(aes(group = Module), level=0.95, lwd=1, linetype = 2)+
  geom_point(size=2)+
  theme_bw()
p_dw <- p_dw + geom_point(size=2) +
  scale_color_manual(values = mycol)+
  scale_shape_manual(values=c(1,rep(c(0:2,5:6,9:10,11:12,14), times=6)))+
  theme_bw()
print(p_dw)

p_du<-plot_ordination(module.rel, pcoa_du, type="taxa", color = 'Module',
                   title="GUniFrac Unweighted UniFrac") + geom_point(size=2)+ theme_bw()
print(p_du)




#Permanova - Distance based multivariate analysis of variance

species.scores <- as.data.frame(scores(ord_rel, "species"))

adonis(as.dist(d5) ~ as.matrix(tax_table(physeq)[,"Module"]))

PermanovaG(unifracs[, , c("d_0", "d_0.5", "d_1")]  ~ as.matrix(sample_data(physeq)[,"Forest"]))

# > PermanovaG(unifracs[, , c("d_0", "d_0.5", "d_1")]  ~ as.matrix(sample_data(physeq_subset)[,"Country"]))
# $aov.tab
# F.Model p.value
# as.matrix(sample_data(physeq_subset)[, "Country"]) 14.55806   0.001

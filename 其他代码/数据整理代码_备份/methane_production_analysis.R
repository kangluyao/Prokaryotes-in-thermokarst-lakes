#methane
methane_taxs <- c('Methanopyrales', 'Methanococcales', 'Methanobacteriales', 'Methanomicrobiales',
                  'Methanosarcinales', 'Methanocellales', 'Methanomassiliicoccales')
methane_phylo <- subset_taxa(water_physeq, Order %in% methane_taxs)
ref.seqs <- refseq(methane_phylo)
tree <- phy_tree(methane_phylo)
tax <- tax_table(methane_phylo)
otu.table <- otu_table(methane_phylo)

write.csv(tax, './taxonomy.csv')
write.csv(otu.table, './methano_table.csv')
write.tree(tree, './methy_tree.nwk')
writeXStringSet(ref.seqs, './Methy.ref.seqs.fasta', append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")


#comunity composition
##determine the genus compositions within each family##
methane_phylo_rel <- subset_taxa(water.rel, Order %in% methane_taxs)
genus.phy <- tax_glom(methane_phylo_rel, taxrank="Genus")
genus.ra.table <- otu_table(genus.phy)
MRA <- rowMeans(genus.ra.table)
group <- tax_table(genus.phy)[,c(4,6)]
genus.mra.table <- data.frame(group,MRA)

#arrange the genuss table
library(tidyr)
genus.mra.table <- genus.mra.table %>% spread(Genus, MRA)
genus.mra.table[is.na(genus.mra.table)] <- 0
rownames(genus.mra.table)<-genus.mra.table$Order
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


setwd('E:/thermokast_lakes/water_microbes/tibet_dada2_asv/data/methy_data/')
Methanotrophy_table <- subset_taxa(water_physeq, Order == "Methylococcales" | Family == 'Methylacidiphilaceae' | Family == 'Beijerinckiaceae'| Family == 'Methylocystaceae' | Family == 'Methylomirabilaceae')
ref.seqs <- refseq(Methanotrophy_table)
tree <- phy_tree(Methanotrophy_table)
tax <- tax_table(Methanotrophy_table)
otu.table <- otu_table(Methanotrophy_table)

write.csv(tax, './taxonomy.csv')
write.csv(otu.table, './methano_table.csv')
write.tree(tree, './methy_tree.nwk')
writeXStringSet(ref.seqs, './Methy.ref.seqs.fasta', append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

##file for iTol
methy_tax <- read.csv('tax.csv',row.names=1, header = T)
methy_asvs <- methy_tax$ASV
methy_phy <- subset_taxa(water_physeq, ASV %in% methy_asvs)

ref.seqs <- refseq(methy_phy)
tree <- phy_tree(methy_phy)
tax.table <- tax_table(as.matrix(methy_tax))
otu.table <- otu_table(methy_phy)
env.data <- read.csv("E:/thermokast_lakes/water_microbes/tibet_dada2_asv/data/metadata.csv", row.names=1, header = T)
methy_phy <- phyloseq(otu.table, tax.table, env.table, tree, ref.seqs)


#write the tree, otu table and taxonomy table
tree <- phy_tree(methy_phy)
tax <- tax_table(methy_phy)
otu.table <- otu_table(methy_phy)
ref.seqs <- refseq(methy_phy)
sample.data <- sample_data(methy_phy)

setwd('E:/thermokast_lakes/water_microbes/tibet_dada2_asv/data/methano_data/')

write.csv(tax, './taxonomy.csv')
write.csv(otu.table, './otu.table.csv')
write.csv(sample.data, './sample.data.csv')
write.tree(tree, './tree.nwk')
writeXStringSet(ref.seqs, './ref.seqs.fasta', append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

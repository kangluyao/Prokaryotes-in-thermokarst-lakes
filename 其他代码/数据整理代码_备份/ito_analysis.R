
nodes <- names(V(core.net))
module.phy <- prune_taxa(nodes, core.phy)

#group the otus into module
wtc <- cluster_fast_greedy(core.net)
Module <- data.frame(Tip = wtc$names, Module = paste('Module', wtc$membership, sep = ''))

##file for iTol
module.ito.tree<-phy_tree(module.phy)
write.tree(module.ito.tree, file = './module_data/ito/module.ito.tree.tre', append = FALSE,
           digits = 10, tree.names = FALSE)
module.anno.1<-data.frame(as(tax_table(module.phy), "matrix"))
module.anno.1<-data.frame(Tip=rownames(module.anno.1),module.anno.1)
module.anno.1 <- merge(module.anno.1, Module, by = 'Tip', all.x = TRUE)

write.table(module.anno.1, file = "./module_data/ito/module.anno.1.tsv", row.names=FALSE, sep="\t")


module_abund_table<-data.frame(t(otu_table(module.phy)))
module_abund_table<-data.frame(Tip=rownames(module_abund_table),module_abund_table)
write.table(module_abund_table, file = "./module_data/ito/module.anno.2.tsv", row.names=FALSE, sep="\t")























source('C:/Program Files/R/R-3.6.0/library/table2itol/table2itol.R')


#read.table(file = 'drug_info.tsv', sep = '\t', header = TRUE)

create_itol_files(infiles = c("./module_data/ito/module.anno.1.tsv",
                              './module_data/ito/module.anno.2.tsv'),
                  identifier = "Tip", label = "OTU", na.strings = "X")

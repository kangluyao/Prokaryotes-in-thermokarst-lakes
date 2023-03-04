setwd('E:/thermokast_lakes/water_microbes/meta_analysis/')
#=============================================
otu<-read.table(file="./data/data_from_qiime2/otu_table.tsv", sep= '\t', header=TRUE)
tax<-read.table(file="./data/data_from_qiime2/taxonomy.tsv", sep= '\t', header= TRUE,  fill = TRUE)
merged_file<- merge(otu, tax, by.x= c("OTU_ID"), by.y= c("OTU_ID"))
write.table(merged_file, file = './data/data_from_qiime2/combined_otu_tax',sep = '\t',
            col.names = T, row.names = F)
##将combined_otu_tax在excel中整理成OTU_TABLE和Taxonomy的csv文档
#==============================================
#load libraries
library(ggplot2)
library(phyloseq)
library(ape)
library(Biostrings)


#read in otu table
meta.otu.table <- read.csv("./data/data_from_qiime2/otu_table.csv",sep=",", row.names=1)
meta.otu.table <- as.matrix(meta.otu.table)

#read in taxonomy
#seperated by kingdom phylum class order family genus species 
meta.taxonomy <- read.csv("./data/data_from_qiime2/taxonomy.csv",sep=",",row.names=1)
meta.taxonomy <- as.matrix(meta.taxonomy)

#read in metadata
metadata <- read.csv("./data/metadata.csv", row.names=1, header = T)

# read in tree
meta.phy.tree <- read_tree("./data/data_from_qiime2/tree.nwk")

#read in represent dna sequences
meta.ref.seqs <- readDNAStringSet(file = "./data/data_from_qiime2/dna-sequences.fasta",
                                  format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#import as phyloseq objects
meta.otu.table <- otu_table(meta.otu.table, taxa_are_rows = TRUE)
meta.tax.table <- tax_table(meta.taxonomy)
meta.table <- sample_data(metadata)

#merge into one phyloseq object
meta_physeq <- phyloseq(meta.otu.table, meta.tax.table,
                        meta.table, meta.phy.tree, meta.ref.seqs)
meta_physeq



#change HashCode ID to OTUIDs
taxa_names(meta_physeq) <- paste0("OTU", seq(ntaxa(meta_physeq)))

#write the tree, otu table and taxonomy table
meta.tree <- phy_tree(meta_physeq)
meta.taxonomy <- tax_table(meta_physeq)
meta.otu.table <- otu_table(meta_physeq)
meta.ref.seqs <- refseq(meta_physeq)
sample.data <- sample_data(meta_physeq)

write.csv(meta.taxonomy, './data/meta_data/meta_taxonomy.csv')
write.csv(meta.otu.table, './data/meta_data/meta_otu_table.csv')
write.csv(sample.data, './data/meta_data/sample.data.csv')
write.tree(meta.tree, './data/meta_data/meta_tree.nwk')
writeXStringSet(meta.ref.seqs, './data/meta_data/meta_ref_seqs.fasta', append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

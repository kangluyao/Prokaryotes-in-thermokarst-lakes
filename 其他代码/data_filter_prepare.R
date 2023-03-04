#====================
setwd('E:/thermokast_lakes/water_microbes/tibet_dada2_asv/')
otu<-read.table(file="./data/data_from_qiime2/ASV_table.tsv", sep= '\t', header=TRUE)
tax<-read.table(file="./data/data_from_qiime2/taxonomy.tsv", sep= '\t', header= TRUE, fill = TRUE)
merged_file<- merge(otu, tax, by.x= c("OTU_ID"), by.y= c("OTU_ID"))
write.table(merged_file, file = './data/data_from_qiime2/combined_otu_tax',sep = '\t',
            col.names = T, row.names = F)

#==============================================
#load libraries
library(ggplot2)
library(phyloseq)
library(ape)
library(Biostrings)


#read in otu table
otu.table <- read.csv("./data/data_from_qiime2/otu_table.csv",sep=",", row.names=1)
otu.table <- as.matrix(otu.table)

#read in taxonomy
#seperated by kingdom phylum class order family genus species 
taxonomy <- read.csv("./data/data_from_qiime2/taxonomy.csv",sep=",",row.names=1)
taxonomy <- as.matrix(taxonomy)

#read in metadata
metadata <- read.csv("./data/metadata.csv", row.names=1, header = T)

# read in tree
phy_tree <- read_tree("./data/data_from_qiime2/tree.nwk")

#read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "./data/data_from_qiime2/dna-sequences.fasta",
                                   format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#import as phyloseq objects
otu.table <- otu_table(otu.table, taxa_are_rows = TRUE)
tax.table <- tax_table(taxonomy)
meta.table <- sample_data(metadata)

#merge into one phyloseq object
water_physeq <- phyloseq(otu.table, tax.table, meta.table, phy_tree, ref_seqs)
water_physeq



#change HashCode ID to ASVIDs
taxa_names(water_physeq) <- paste0("ASV", seq(ntaxa(water_physeq)))

#write the tree, otu table and taxonomy table
tree <- phy_tree(water_physeq)
tax <- tax_table(water_physeq)
otu.table <- otu_table(water_physeq)
ref.seqs <- refseq(water_physeq)
sample.data <- sample_data(water_physeq)

write.csv(tax, './data/total_taxa_data/taxonomy.csv')
write.csv(otu.table, './data/total_taxa_data/otu.table.csv')
write.csv(sample.data, './data/total_taxa_data/sample.data.csv')
write.tree(tree, './data/total_taxa_data/tree.nwk')
writeXStringSet(ref.seqs, './data/total_taxa_data/ref.seqs.fasta', append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")



setwd('E:/thermokast_lakes/water_microbes/lakeshore/data_from_qiime2/data_unrarefied/')
#=============================================
# you have otu_table.txt and taxonomy.tsv
# open them up in text edit and change #OTUID to OTUID in otu_table.txt
# and change Feature ID to OTUID in taxonomy.tsv

#	read in	OTU	table
otu<-read.table(file="./otu_table.tsv", sep= '\t', header=TRUE)
# read in taxonomy table
tax<-read.table(file="./taxonomy.tsv", sep= '\t', header= TRUE, fill = T)
# merge files
merged_file<- merge(otu, tax, by.x= c("OTU_ID"), by.y= c("OTU_ID"))
# note: number of rows should equal your shortest file length, drops taxonomy for OTUs that don't exist in your OTU table
# output merged .txt file
write.table(merged_file, file = './combined_otu_tax',sep = '\t',
            col.names = T, row.names = F)
# It seems tedious but you need to open the merged .txt file in excel and split into two files: 
# one for taxonomy (containing only the columns OTUID and taxonomic info) and 
# the other for the OTU matrix (containing only OTUID and abundances in each sample).
# Note: for the taxonomy file, you need to use data —> text-to-columns 
# in Excel and separate on semicolon to get columns for 
# kingdom, phylum, class, etc… once you make these two separate files in excel, 
# save each as a .csv

#====================
setwd('E:/thermokast_lakes/water_microbes/OTU_97')
otu<-read.table(file="./data/otu_table.tsv", sep= '\t', header=TRUE)
tax<-read.table(file="./data/taxonomy.tsv", sep= '\t', header= TRUE)
merged_file<- merge(otu, tax, by.x= c("OTU_ID"), by.y= c("OTU_ID"))
write.table(merged_file, file = './data/combined_otu_tax',sep = '\t',
            col.names = T, row.names = F)

#==============================================
#load libraries
library(ggplot2)
library(phyloseq)
library(ape)
library(Biostrings)


#read in otu table
otu.table <- read.csv("./otu.csv",sep=",", row.names=1)
otu.table <- as.matrix(otu.table)

#read in taxonomy
#seperated by kingdom phylum class order family genus species 
taxonomy <- read.csv("./taxonomy.csv",sep=",",row.names=1)
taxonomy <- as.matrix(taxonomy)

#read in metadata
metadata <- read.csv("./data/data_from_qiime2/metadata.csv", row.names=1, header = T)

# read in tree
phy_tree <- read_tree("./tree.nwk")

#read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "./dna-sequences.fasta",
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#import as phyloseq objects
otu.table <- otu_table(otu.table, taxa_are_rows = TRUE)
tax.table <- tax_table(taxonomy)
meta.table <- sample_data(metadata)

#merge into one phyloseq object
lakeshores_phylo <- phyloseq(otu.table, tax.table, phy_tree, ref_seqs)
lakeshores_phylo


#change HashCode ID to ASVIDs
taxa_names(lakeshores_phylo) <- paste0("OTU", seq(ntaxa(lakeshores_phylo)))

#write the tree, otu table and taxonomy table
total.tree <- phy_tree(lakeshores_phylo)
total.taxonomy <- tax_table(lakeshores_phylo)
total.otu.table <- otu_table(lakeshores_phylo)
total.ref.seqs <- refseq(lakeshores_phylo)

write.csv(total.taxonomy, '../taxonomy.csv')
write.csv(total.otu.table, '../otu_table.csv')
write.tree(total.tree, '../tree.nwk')
writeXStringSet(total.ref.seqs, '../ref_seqs.fasta', append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")


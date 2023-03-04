setwd('E:/thermokast_lakes/water_microbes/lakeshore/data/')
#conduct a phyloseq project
#read in metadata
metadata <- read.csv("./sample_data.csv",
                     header = T, row.names = 1)
#read in otu table
otu <- read.csv("./otu_table.csv",sep=",", row.names=1)
otu <- as.matrix(otu)
#read in taxonomy
taxonomy <- read.csv("./taxonomy.csv",sep=",",row.names=1)
taxonomy <- as.matrix(taxonomy)

# read in tree
phy_tree <- read_tree("./tree.nwk")

#read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "./ref_seqs.fasta",
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#import as phyloseq objects
otu <- otu_table(otu, taxa_are_rows = TRUE)
tax <- tax_table(taxonomy)
meta.table <- sample_data(metadata)

#merge into one phyloseq object
lakeshores_phylo <- phyloseq(otu, tax, phy_tree, meta.table, ref_seqs)
lakeshores_phylo
set.seed(1234)
lakeshores_phylo_even = rarefy_even_depth(lakeshores_phylo, sample.size = 33859, replace = TRUE)



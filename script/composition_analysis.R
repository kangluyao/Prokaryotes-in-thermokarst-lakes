# set work dorectory
setwd('E:/thermokast_lakes/water_microbes/')
# loading packages
# ipak function: install and load multiple R packages.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
pkgs <- c('ggplot2', 'ape', 'Biostrings', 'vegan', 'tidyverse',
          'reshape', 'cowplot', 'microbiome', 'RColorBrewer',
          'ggpubr', 'dplyr', 'phyloseq')
ipak(pkgs)
# lapply(PKGs, require, character.only = TRUE, warn.conflicts = FALSE)

#conduct a phyloseq project
#read in metadata
metadata <- read.delim("./meta_analysis/data1/meta_data/sample_data.txt",
                     header = T, row.names = 1)

#read in otu table
meta.otu.table <- read.delim("./meta_analysis/data1/meta_data/meta_otu_table.txt",
                           header = T, row.names = 1, stringsAsFactors = F)
meta.otu.table <- as.matrix(meta.otu.table)

#read in taxonomy
meta.taxonomy <- read.delim("./meta_analysis/data1/meta_data/taxonomy.txt",
                            header = T, row.names = 1, stringsAsFactors = F)
meta.taxonomy <- as.matrix(meta.taxonomy)

# read in tree
meta.phy.tree <- read_tree("./meta_analysis/data1/meta_data/meta_tree.nwk")

#read in represent dna sequences
meta.ref.seqs <- readDNAStringSet(file = "./meta_analysis/data1/meta_data/meta_ref_seqs.fasta",
                                  format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
meta.otu.table <- otu_table(meta.otu.table, taxa_are_rows = TRUE)
meta.tax.table <- tax_table(meta.taxonomy)
meta.table <- sample_data(metadata)
meta_physeq <- phyloseq(meta.tax.table, meta.otu.table, meta.table, meta.phy.tree, meta.ref.seqs)
meta_physeq
#Convert to relative abundance
#meta_physeq_rel = phyloseq::transform_sample_counts(meta_physeq, function(x){x / sum(x)})
#phyloseq::otu_table(meta_physeq)[1:5, 1:5]
meta_physeq_rel <- microbiome::transform(meta_physeq, "compositional")
meta_table <- as(sample_data(meta_physeq), "data.frame")

pa_phylo <- subset_samples(meta_physeq, Region == "Pan-Arctic")
pa_phylo <- prune_taxa(taxa_sums(pa_phylo) > 0, pa_phylo) 
tp_phylo <- subset_samples(meta_physeq, Region == "Tibetan Plateau")
tp_phylo <- prune_taxa(taxa_sums(tp_phylo) > 0, tp_phylo) 

# determine the OTU numbers within each dominant phylum
## data preparation
otu_num_in_phylum <- data.frame(table(tax_table(meta_physeq)[,"Phylum"]))
otu_num_in_phylum <- otu_num_in_phylum %>% arrange(desc(Freq))
otu_num_in_phylum <- rbind(otu_num_in_phylum[1:11, ], data.frame(Var1 = c('Others'), Freq = sum(otu_num_in_phylum[-c(1:11), 2])))

otu_num_in_phylum <- data.frame(Phylum = otu_num_in_phylum$Var1, otu_num = otu_num_in_phylum$Freq ,
                                prop = otu_num_in_phylum$Freq/sum(otu_num_in_phylum$Freq)*100)
otu_count.data <- otu_num_in_phylum %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
otu_count.data$Phylum <- factor(otu_count.data$Phylum, ordered = T, levels = otu_num_in_phylum$Phylum)

## pie plot 
### Define the colors you want
mycols <- c("#89c5da", "#ffc15c", "#74d944", "#CE50CA", "#5e738f", "#C0717C", "#CBD5ec", "#5F7FC7", 
                     "#00718b", "#d3d93e", "#169e77", "#d95f02", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                     "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
                     "#8A7C64", "#599861")

pie_for_otu_num_phylum <- ggplot(otu_count.data, aes(x="", y=prop, fill=reorder(Phylum, -lab.ypos))) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  geom_text(aes(x = 1.35, y = lab.ypos, 
                label = paste0(otu_num, ' (', round(prop, 1), '%', ')', sep= '')),
            color = "black", size=3) +
  scale_fill_manual('Phylum', values = mycols) +
  guides(fill = guide_legend(reverse=T)) +
  theme_void() +
  theme(legend.position = "left")
pie_for_otu_num_phylum

# determine the Order compositions within top 10 phyla
arrange.tab <- function(phylo, N, taxrank, vect) {
  subphylo <- tax_glom(phylo, taxrank)
  subphylo.rel <- microbiome::transform(subphylo, "compositional")
  ra.tab <- otu_table(subphylo.rel)
  MRA <- rowMeans(ra.tab)
  group <- tax_table(subphylo.rel)[,vect]
  mra.tab <- data.frame(group,MRA)
  colnames(mra.tab) <- c('level1', 'level2', 'MRA')
  #arrange the class table
  mra.tab_level1 = mra.tab %>% group_by(level1) %>% 
    summarise(sum_MRA = sum(MRA)) %>% 
    arrange(desc(sum_MRA))
  top_N_level1 = mra.tab_level1[1:N, ]$'level1'
  top_N_tab = mra.tab[mra.tab$'level1' %in% top_N_level1, ]
  mra.tab_level2 = top_N_tab %>% group_by(level2) %>% 
    summarise(sum_MRA = sum(MRA)) %>% 
    arrange(desc(sum_MRA))
  order_level2 = mra.tab_level2$'level2'
  top_N_tab$'level1' = factor(top_N_tab$'level1', ordered = T, levels = top_N_level1)
  top_N_tab$'level2' = factor(top_N_tab$'level2', ordered = T, levels = rev(order_level2))
  top_N_tab$MRA = top_N_tab$MRA*100
  return(top_N_tab)
}
top10phylum_meta <- arrange.tab(meta_physeq, 10, 'Order', c(2,4))
mra.tab_level1 = top10phylum_meta %>% group_by(level1) %>% 
  summarise(sum_MRA = sum(MRA)) %>% 
  arrange(desc(sum_MRA))
mra.tab_level1
mra.tab_level2 = top10phylum_meta %>% group_by(level2) %>% 
  summarise(sum_MRA = sum(MRA)) %>% 
  arrange(desc(sum_MRA))
mra.tab_level2 [1:20, ]
order_level2 = mra.tab_level2$'level2'
#stack bar plot
taxa_barplot <- ggplot(top10phylum_meta, aes(fill=level2, y=MRA, x=level1)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual('Order', breaks = order_level2[1:15], 
                    values = rep(c(rev(mycols[1:12]), mycols[-c(1:12)]), 20)[1:nrow(top10phylum_meta)]) + #only the top 10 phylum and top 10 order are showed
  labs(x = 'Phylum', y = 'Mean relative abundance (%)') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 55)) +
  theme_classic()+
  theme(legend.position = c(0.6,0.6),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
taxa_barplot

# plot the community composition for TP and PA at order level
subphylo <- tax_glom(meta_physeq, taxrank = 'Order')
subphylo.rel  = transform_sample_counts(subphylo, function(x) x / sum(x))
ntaxa(subphylo.rel)
# meta_physeq_rel <- microbiome::transform(meta_physeq, "compositional")
# meta.com.cla <- microbiome::aggregate_rare(meta_physeq_rel, level = "Order", 
#                                            detection = 1/100, prevalence = 10/100)

ra.tab <- otu_table(subphylo.rel)
sum(ra.tab[, 1])
subtaxa_tab <- tax_table(subphylo.rel)[, 4]

# set the plot theme
mycols <- c("#d95f02", "#169e77", "#d3d93e", "#00718b", "#5F7FC7", "#CBD5ec", "#C0717C", "#5e738f", 
            "#CE50CA", "#74d944", "#ffc15c", "#89c5da", "#D7C1B1", "#689030", "#AD6F3B", "grey", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")

main_theme = theme_linedraw() + 
  theme(panel.grid=element_blank(), 
        strip.text = element_text(colour = 'black', size = 12),
        strip.background = element_rect(colour = 'grey', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(colour = 'black', size = 12),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white"))

order_tax_table <- data.frame(subtaxa_tab, ra.tab) %>% 
  mutate(MRA = rowMeans(select(., colnames(ra.tab)))) %>%
  arrange(desc(MRA)) %>% dplyr::top_n(15, MRA) %>%
  select(., -c('MRA')) %>% 
  bind_rows(summarise_all(., ~if(is.numeric(.)) 1-sum(.) else "Others")) %>%
  mutate(Order = factor(Order, levels = Order)) %>%
  pivot_longer(cols = -c(Order), names_to = "Sample_name", values_to = 'rel_abun') %>%
  right_join(data.frame(Sample_name = rownames(metadata), Region = metadata$Region), by = c("Sample_name")) %>%
  select(., -c('Sample_name')) %>% 
  group_by(Region, Order) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))

print(order_tax_table, n = 32)
# write.csv(order_tax_table, "E:/thermokast_lakes/water_microbes/meta_analysis/results/tables/order_tax_table.csv")
# bar plot at the order level
bar_plot <- order_tax_table %>%
  ggplot(aes(x = Region, y = 100*rel_abun, fill = Order))+
  geom_bar(stat = "identity") +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_discrete(limits = c('Pan-Arctic', 'Tibetan Plateau')) +
  scale_fill_manual(values =  mycols) +
  labs(x = 'Region', y = 'Mean relative abundance (%)') +
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  main_theme +
  guides(fill = guide_legend(ncol = 5)) +
  coord_flip()
bar_plot
# Combining Plots
library(cowplot)
compositional_plot <- ggdraw() +
  draw_plot(pie_for_otu_num_phylum, x = 0, y = 1/2, width = 0.5, height = 1/2) +
  draw_plot(taxa_barplot, x = 0.5, y = 1/2, width = 0.5, height = 1/2) +
  draw_plot(bar_plot, x = 0, y = 0, width = 1, height = 1/2) +
  draw_plot_label(label = c("A", "B", "C"), size = 14,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))
compositional_plot
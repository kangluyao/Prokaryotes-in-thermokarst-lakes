# set work dorectory
setwd('E:/thermokast_lakes/water_microbes/')
# loading packages
PKGs <- c('ggplot2', 'ape', 'Biostrings', 'vegan', 'tidyverse',
          'reshape', 'cowplot', 'microbiome', 'RColorBrewer',
          'ggpubr', 'dplyr')

lapply(PKGs, require, character.only = TRUE, warn.conflicts = FALSE)

#conduct a phyloseq project
#read in metadata
metadata <- read.csv("./meta_analysis/data1/meta_data/sample_data.csv",
                     header = T, row.names = 1)

#read in otu table
meta.otu.table <- read.csv("./meta_analysis/data1/meta_data/meta_otu_table.csv",
                           header = T, row.names = 1, stringsAsFactors = F)
meta.otu.table <- as.matrix(meta.otu.table)

#read in taxonomy
meta.taxonomy <- read.csv("./meta_analysis/data1/meta_data/taxonomy.csv", sep=",", row.names=1)
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
  phylo.rel <- microbiome::transform(phylo, "compositional")
  ra.tab <- otu_table(phylo.rel)
  MRA <- rowMeans(ra.tab)
  group <- tax_table(subphylo.rel)[, vect]
  mra.tab <- data.frame(group, MRA)
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

# Combining Plots
library(cowplot)
compositional_plot <- ggdraw() +
  draw_plot(pie_for_otu_num_phylum, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(taxa_barplot, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 14,
                  x = c(0, 0.5), y = c(1, 1))

compositional_plot







# archea analysis
archaea_phylo <- subset_taxa(meta_physeq, Kingdom == 'Archaea', level = "Order")
top10phylum_archaea <- arrange.tab(archaea_phylo, 10, 'Order', c(2,4))
mra.archaea_level1 = top10phylum_archaea %>% group_by(level1) %>% 
  summarise(sum_MRA = sum(MRA)) %>% 
  arrange(desc(sum_MRA))
mra.archaea_level1
mra.archaea_level2 = top10phylum_archaea %>% group_by(level2) %>% 
  summarise(sum_MRA = sum(MRA)) %>% 
  arrange(desc(sum_MRA))
mra.archaea_level2

# plot the community composition for TP and PA at order level

##prune out Order below 1% in each sample and prevalence lower than 10/100 at class level
meta_physeq_rel <- microbiome::transform(meta_physeq, "compositional")
meta.com.ord <- microbiome::aggregate_rare(meta_physeq_rel, level = "Order", 
                                           detection = 1/100, prevalence = 10/100)

plot.composition.relAbun <- microbiome::plot_composition(meta.com.ord, 
                                                         average_by = "Region", 
                                                         otu.sort = "abundance") +
  scale_fill_manual(values = mycols) + 
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_y_continuous(labels = function(x) paste0(x*100), expand = c(0, 0)) + 
  labs(x = NULL, y = 'Mean relative abundance (%)') +
  theme_bw() +
  guides(fill=guide_legend(ncol=1)) +
  theme(legend.position = 'right',
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 14, colour = 'black'),
        legend.text = element_text(size = 12, colour = 'black'))
plot.composition.relAbun 



## determine the composition within gammaproteobacteria
gammaproteobacteria_phylo <- subset_taxa(meta_physeq, Class == 'Gammaproteobacteria', level = "Order")
gammaproteobacteria_phylo <- tax_glom(gammaproteobacteria_phylo, taxrank="Order")
gammaproteobacteria_phylo_rel <- microbiome::transform(gammaproteobacteria_phylo, "compositional")
gammaproteobacteria_phylo_rel <- aggregate_rare(gammaproteobacteria_phylo_rel, level = "Order", detection = 1/100, prevalence = 10/100)

plot.gammaproteobacteria.relAbun <- microbiome::plot_composition(gammaproteobacteria_phylo_rel, 
                                                                 average_by = "Region", 
                                                                 otu.sort = "abundance") +
  scale_fill_brewer("Order", palette = "Paired") + 
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_y_continuous(label = scales::percent, expand = c(0, 0)) + 
  labs(x = NULL, y = 'Relative abundance') +
  theme_bw() +
  theme(legend.position =  'bottom',
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 12, colour = 'black'),
        legend.text = element_text(size = 10, colour = 'black'))

## determine the composition within actinobacteria
actinobacteria_phylo <- subset_taxa(meta_physeq, Class == 'Actinobacteria', level = "Order")
actinobacteria_phylo <- tax_glom(actinobacteria_phylo, taxrank="Order")
actinobacteria_phylo_rel <- microbiome::transform(actinobacteria_phylo, "compositional")
actinobacteria_phylo_rel <- aggregate_rare(actinobacteria_phylo_rel, level = "Order", detection = 1/100, prevalence = 10/100)

plot.actinobacteria.relAbun <- microbiome::plot_composition(actinobacteria_phylo_rel, 
                                                            average_by = "Region", 
                                                            otu.sort = "abundance") +
  scale_fill_brewer("Order", palette = "Set1") + 
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_y_continuous(label = scales::percent, expand = c(0, 0)) + 
  labs(x = NULL, y = 'Relative abundance') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 12, colour = 'black'),
        legend.text = element_text(size = 10, colour = 'black'))

## determine the composition within bacteroidia
bacteroidia_phylo <- subset_taxa(meta_physeq, Class == 'Bacteroidia', level = "Order")
bacteroidia_phylo <- tax_glom(bacteroidia_phylo, taxrank="Order")
bacteroidia_phylo_rel <- microbiome::transform(bacteroidia_phylo, "compositional")
bacteroidia_phylo_rel <- aggregate_rare(bacteroidia_phylo_rel, level = "Order", detection = 1/100, prevalence = 10/100)

plot.bacteroidia.relAbun <- microbiome::plot_composition(bacteroidia_phylo_rel, 
                                                         average_by = "Region", 
                                                         otu.sort = "abundance") +
  scale_fill_brewer("Order", palette = "Accent") + 
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_y_continuous(label = scales::percent, expand = c(0, 0)) + 
  labs(x = NULL, y = 'Relative abundance') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 12, colour = 'black'),
        legend.text = element_text(size = 10, colour = 'black'))
## arrange the plots
domin_class_plot <- plot_grid(plot.gammaproteobacteria.relAbun, plot.actinobacteria.relAbun, plot.bacteroidia.relAbun,
                              labels = c('(a)', '(b)', '(c)'), 
                              ncol = 3, nrow = 1, 
                              label_x = .01, label_y = 1.01, hjust = 0, 
                              label_size=14, align = "v")
domin_class_plot


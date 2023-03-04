# loading packages
PKGs <- c('ggplot2', 'ape', 'Biostrings', 'vegan', 
              'reshape', 'cowplot', 'microbiome', 'RColorBrewer',
              'ggpubr', 'dplyr')

lapply(PKGs, require, character.only = TRUE, warn.conflicts = FALSE)

# set work dorectory
setwd('E:/thermokast_lakes/water_microbes/')
#conduct a phyloseq project
#read in metadata
metadata <- read.csv("./meta_analysis/data/meta_data/sample_data.csv",
                     header = T, row.names = 1)

#read in otu table
meta.otu.table <- read.csv("./meta_analysis/data/meta_data/meta_otu_table.csv",
                           header = T, row.names = 1, stringsAsFactors = F)
meta.otu.table <- as.matrix(meta.otu.table)

#read in taxonomy
meta.taxonomy <- read.csv("./meta_analysis/data/meta_data/meta_taxonomy.csv",sep=",",row.names=1)
meta.taxonomy <- as.matrix(meta.taxonomy)

# read in tree
meta.phy.tree <- read_tree("./meta_analysis/data/meta_data/meta_tree.nwk")

#read in represent dna sequences
meta.ref.seqs <- readDNAStringSet(file = "./meta_analysis/data/meta_data/meta_ref_seqs.fasta",
                                  format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
meta.otu.table <- otu_table(meta.otu.table, taxa_are_rows = TRUE)
meta.tax.table <- tax_table(meta.taxonomy)
meta.table <- sample_data(metadata)
meta_physeq <- phyloseq(meta.tax.table, meta.otu.table, meta.table, meta.phy.tree, meta.ref.seqs)
meta_physeq
meta_physeq_rel <- microbiome::transform(meta_physeq, "compositional")
meta_table <- as(sample_data(meta_physeq), "data.frame")

pa_phylo <- subset_samples(meta_physeq, Region == "Pan-Arctic")
pa_phylo <- prune_taxa(taxa_sums(pa_phylo) > 0, pa_phylo) 
tp_phylo <- subset_samples(meta_physeq, Region == "Tibetan Plateau")
tp_phylo <- prune_taxa(taxa_sums(tp_phylo) > 0, tp_phylo) 


# melt to long format (for ggploting) 
# prune out class below 1% in each sample and prevalence lower than 10/100
meta.com.cla <- aggregate_rare(meta_physeq_rel, level = "Class", detection = 1/100, prevalence = 10/100)
plot.composition.relAbun <- microbiome::plot_composition(meta.com.cla, 
                                                         average_by = "Region", 
                                                         otu.sort = "abundance") +
  scale_fill_brewer("Class", palette = "Paired") + 
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_y_continuous(label = scales::percent, expand = c(0, 0)) + 
  labs(x = NULL, y = 'Relative abundance') +
  theme_bw() +
  theme(legend.position = 'left',
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 14, colour = 'black'),
        legend.text = element_text(size = 12, colour = 'black'))

# determine the composition within gammaproteobacteria
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

# determine the composition within actinobacteria
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

# determine the composition within bacteroidia
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
#arrange the plots
domin_class_plot <- plot_grid(plot.gammaproteobacteria.relAbun, plot.actinobacteria.relAbun, plot.bacteroidia.relAbun,
                              labels = c('(a)', '(b)', '(c)'), 
                              ncol = 3, nrow = 1, 
                              label_x = .01, label_y = 1.01, hjust = 0, 
                              label_size=14, align = "v")
domin_class_plot


# alpha diversity
diversity <- estimate_richness(meta_physeq, measures = c("Chao1", 'Shannon', 'Simpson'))
meta_diversity <- cbind(diversity, Region = meta_table$Region, 
                        Site = meta_table$Site)
# test the difference
library(lme4)
library(lmerTest)
library(multcomp)
mode1 <- lmer(Chao1 ~ Region + (1|Site), meta_diversity)
summary(mode1)
mode2 <- lmer(Shannon ~ Region + (1|Site), meta_diversity)
summary(mode2)
mode3 <- lmer(Simpson ~ Region + (1|Site), meta_diversity)
summary(mode3)

# diversity plot
melted <- melt(meta_diversity[,c("Chao1", "Shannon", 'Simpson', "Region")], id.vars = c("Region"))
alpha_region_plot <- ggplot(melted, aes(x = Region, y = value, fill = Region)) +
  geom_violin(trim=T, width=0.5, aes(fill = Region), colour = "#000000") +
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_fill_manual(values= c('#d95f02', '#1b9e77')) +
  geom_boxplot(width=0.1, fill="white", colour = "#000000") +
  labs(x = NULL, y = 'Alpha diversity') +
  facet_wrap(~variable, scales = 'free') +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12, 
                                   angle = 45, hjust = 1),
        panel.grid = element_blank(), legend.position = 'none')
alpha_region_plot

alpha_site_plot <- plot_richness(meta_physeq, x="Sitegroup", color = 'Sitegroup', 
                                   measures=c("Chao1", "Shannon", 'Simpson'))+
  geom_boxplot(aes(fill = Sitegroup), alpha=0.2)+
  theme(strip.text = element_text(size = 14),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 90, hjust = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
alpha_site_plot

#NMDS PLOT AND PERMANOVA ANALYSIS
meta.ord <- ordinate(meta_physeq, "NMDS", "bray")
meta.ord
NMDS_plot <- plot_ordination(meta_physeq, meta.ord, type="samples", color="Region") +
  geom_point(size = 2.5) + 
  scale_color_manual(values = c('#d95f02', '#1b9e77')) +
  #stat_ellipse(type = "norm", linetype = 1) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill = Region)) +
  theme_bw() +
  theme(legend.position = c(0.85,0.88),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white"))
metadata <- as(sample_data(meta_physeq), "data.frame")
adonis2(phyloseq::distance(meta_physeq, method = "bray") ~ Region,
       data = metadata)
#nmds analysia for tp
meta.tp.ord <- ordinate(tp_phylo, "NMDS", "bray")
meta.tp.ord
NMDS_tp_plot <- plot_ordination(tp_phylo, meta.tp.ord, type="samples", color="Sitegroup1") +
  geom_point(size = 2.5) + 
  #scale_color_manual(values = c('#d95f02', '#1b9e77')) +
  #stat_ellipse(type = "norm", linetype = 1) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill = Sitegroup1)) +
  theme_bw() +
  theme(legend.position = c(0.85,0.88),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white"))
metadata_tp <- as(sample_data(tp_phylo), "data.frame")
adonis2(phyloseq::distance(tp_phylo, method = "bray") ~ Sitegroup1,
       data = metadata_tp)

alpha_NMDS_plot <- ggdraw() +
  draw_plot(plot.composition.relAbun, x = 0, y = 0, width = 2/5, height = 1) +
  draw_plot(NMDS_plot, x = 2.3/5, y = 2/5, width = 2.7/5, height = 3/5) +
  draw_plot(alpha_region_plot, x = 2.3/5, y = 0, width = 2.7/5, height = 2/5) +
  draw_plot_label(label = c("(a)", "(b)", '(c)'), size = 14,
                  x = c(0, 2.3/5, 2.3/5), y = c(1, 1, 2/5))
alpha_NMDS_plot
#diff class using microeco package
library(microeco)
meco_df <- phyloseq2meco(meta_physeq)
#calculate the abundance table
m1 <- meco_df$cal_abund()

#rf: which depends on the random forest[14, 15] and the non-parametric test. 
# use Genus level for parameter rf_taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
m1 <- trans_diff$new(dataset = meco_df, method = "rf", 
                     group = "Region", rf_taxa_level = "Order")

# m1$res_rf is the result stored in the object
# plot the result
m1 <- m1$plot_diff_abund(use_number = 1:20, only_abund_plot = FALSE)
gridExtra::grid.arrange(m1$p1, m1$p2, ncol=2, nrow = 1, widths = c(2,2))
# the middle asterisk represent the significances


#Lefse
ps_genus <- phyloseq::tax_glom(meta_physeq, taxrank = 'Genus')
meco_genus_df <- phyloseq2meco(ps_genus)
#calculate the abundance table
m1_genus <- meco_genus_df$cal_abund()

m1_genus <- trans_diff$new(dataset = meco_genus_df, method = "lefse", 
                           group = "Region", alpha = 0.01, 
                           lefse_subgroup = NULL)
# t1$res_lefse is the LEfSe result
# t1$res_abund is the abundance information
m1_genus$plot_lefse_bar(LDA_score = 4, color_values = c('#d95f02', '#1b9e77'), width = 0.5)
m1_genus$plot_diff_abund(use_number = 1:30, color_values = c('#d95f02', '#1b9e77'))

# we can format the data for Lefse analysis (http://huttenhower.sph.harvard.edu/galaxy)
ps_genus_sel <- subset_taxa(ps_genus, Family == "Sporichthyaceae" |
                              Genus == "Candidatus Aquiluna" |
                              Genus == "Candidatus Limnoluna" |
                              Genus == "Candidatus Planktoluna" |
                              Genus == "Sediminibacterium" |
                              Genus == "Algoriphagus" |
                              Genus == "Emticicia" |
                              Genus == "Flavobacterium" |
                              Genus == "Pedobacter" |
                              Family == "NS11_12marinegroup" |
                              Class == "Firmicutes" |
                              Family == "Caulobacteraceae" |
                              Genus == "Rhodoblastus" |
                              Genus == "Roseiarcus" |
                              Genus == "Sphingorhabdus" |
                              Genus == "Polynucleobacter" |
                              Genus == "Variovorax" |
                              Family == "Rubritaleaceae" |
                              Family == "Verrucomicrobiaceae" |
                              Genus == "uncultured Opitutaebacterium")

phyloseq2lefse <- function(
  ps,
  covars,
  file.name = "lefse_data.txt",
  taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  transpose.otus = TRUE
) {
  
  # grab the taxa table from the phyloseq object and coerce into a matrix
  tax.tbl <- as(phyloseq::tax_table(ps), "matrix")
  tax.tbl <- tax.tbl[, taxa.levels]
  tax.tbl.join <- as.data.frame(do.call(paste, c(as.data.frame(tax.tbl), sep = "|")))
  
  row.names(tax.tbl.join) <- row.names(tax.tbl)
  colnames(tax.tbl.join) <- NULL
  
  # grab the otu table from the phyloseq object
  otu.tbl <- as(phyloseq::otu_table(ps), "matrix")
  colnames(otu.tbl) <- NULL
  tax.otu.tbl <- cbind(tax.tbl.join, otu.tbl)
  
  # grab the sample table from the phyloseq object
  smpl.data <- as(phyloseq::sample_data(ps), "data.frame")
  smpl.data$Sample <- row.names(smpl.data)
  t.smpl.data <- t(smpl.data)
  t.smpl.data <- as.matrix(t.smpl.data[c("Sample", covars), ])
  t.smpl.data <- cbind(rownames(t.smpl.data), t.smpl.data)
  colnames(t.smpl.data) <- colnames(tax.otu.tbl)
  
  final.data <- rbind(t.smpl.data, tax.otu.tbl)
  write.table(final.data, file = file.name, row.names = F, col.names = F, sep = "\t", quote = FALSE)
}

phyloseq2lefse(
  ps = ps_genus_sel,
  covars = 'Region',
  file.name = "./meta_analysis/results/tables/lefse_data_sel_genus.txt",
  taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  transpose.otus = TRUE
)

# Functional analysis
# create object of trans_func
m2 <- trans_func$new(meco_df)
# mapping the taxonomy to the database
# the function can recognize prokaryotes or fungi automatically.
m2$cal_spe_func()
# return m2$res_spe_func, 1 represent function exists, 0 represent no or cannot confirmed.
m2$res_spe_func[1:5, 1:2]
m2$func_group_list
#carbon cycling taxa analysis
# If you want to change the group list, reset the list t2$func_group_list
C_cycle_process <- c("methanotrophy", 
                     "hydrogenotrophic_methanogenesis", "methanogenesis", 
                     "methanol_oxidation", "methylotrophy", "chitinolysis",
                     "cellulolysis", "xylanolysis", "ligninolysis", "fermentation",
                     'aromatic_hydrocarbon_degradation','aromatic_compound_degradation',
                     'hydrocarbon_degradation')
m2$func_group_list$`C-cycle` <- C_cycle_process
# use show_prok_func to see the detailed information of prokaryotic traits
m2$show_prok_func("methanotrophy")
# calculate the percentages for communities
m2$cal_spe_func_perc(use_community = TRUE)
m2$res_spe_func_perc[1:5, 1:2]

meta_c_cycle <- data.frame(m2$res_spe_func_perc[ ,colnames(m2$res_spe_func_perc) %in% c(m2$func_group_list$`C-cycle`)], 
                           m2$sample_table[, c('Site', 'Region')])

vars <- colnames(m2$res_spe_func_perc[ ,colnames(m2$res_spe_func_perc) %in% c(m2$func_group_list$`C-cycle`)])
# test the significance of the functional taxa among tregion with LMM 
library(lme4)
library(lmerTest)
mode <- lapply(vars, function(x) {
  lmer(substitute(i ~ Region + (1|Site), list(i = as.name(x))), data = meta_c_cycle)})
summary.model <- function(model){
  F.value <- anova(model)$`F value`
  p.value <- anova(model)$'Pr(>F)'
  p.stars <- function(p.values) {
    unclass(symnum(p.values, corr = FALSE, 
                   na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")))}
  sig <- p.stars(p.value)
  results<-data.frame(F.value, p.value, sig)
  return(results)
}
df <- NULL
for(i in 1:length(vars)) {
  tmp <- summary.model(mode[[i]])
  if (is.null(df)){
    df <- tmp
  } else {
    df <- rbind(df, tmp)
  }
}
result_carbon <-data.frame(vars, df)
result_carbon

### arrange data for plot
library(plyr)
library(reshape)
melt_df <- melt(meta_c_cycle[, !colnames(meta_c_cycle) %in% 'Site'], id.vars = c('Region'))
new_df <- ddply(melt_df, c('Region','variable'), summarise,
                mean = mean(value), sd = sd(value),
                sem = sd(value)/sqrt(length(value)))
sig <- c('a', 'a', 'a', 'a', 'a', 'a', 'b', 'b', 'b', 'a', 'b', 'a',
         'b', 'b', 'b', 'a', 'b', 'a', 'a', 'a', 'a', 'a', 'a', 'b')
new_df <- cbind(new_df ,sig)
new_df$Region <- factor(new_df$Region, levels = c('Tibetan Plateau', 'Pan-Arctic'))
#plot
carbon_fun_plot <- ggplot(new_df,aes(x = variable, y = mean, fill = Region))+
  geom_bar(position = 'dodge', stat = 'identity', colour = 'black', width = 0.7)+
  scale_fill_manual(values=c('#1b9e77', '#d95f02'))+
  geom_errorbar(aes(ymin = mean, ymax = mean + sem), width=.2,  position = position_dodge(0.7))+
  scale_x_discrete(limits = c("cellulolysis", "chitinolysis", "xylanolysis", 
                              'aromatic_hydrocarbon_degradation','aromatic_compound_degradation',
                              "fermentation", 'hydrocarbon_degradation',
                              "methanotrophy", "methanogenesis", 
                              "methanol_oxidation", "methylotrophy"))+
  #scale_x_discrete(limits = c("methanotrophy", "methanogenesis"))+
  labs(x = 'Carbon cycle', y = 'Mean relative abundance (%)', fill = "Region")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6))+
  geom_text(aes(label = sig, y = mean + sem + 0.02*max(mean)), 
            position = position_dodge(0.7), vjust = 0)+
  theme_bw()+
  theme(legend.position = c(0.7,0.88),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))


#extract asvs involved in carbon cycling
meta_c_cycle_otus <- m2$res_spe_func[ ,c(m2$func_group_list$`C-cycle`)]
meta_c_cycle_otus <- rownames(meta_c_cycle_otus[rowSums(meta_c_cycle_otus) != 0,])

meta_carbon_phy <- subset_taxa(meta_physeq, OTU %in% meta_c_cycle_otus)
meta_carbon_phy_rel <- subset_taxa(meta_physeq_rel, OTU %in% meta_c_cycle_otus)

meta_meco_df <- phyloseq2meco(meta_carbon_phy)
m3 <- meta_meco_df$cal_abund()
#trans_diff class, The third approach is rf, which depends on the random forest[14, 15] and the non-parametric test. 
# use Genus level for parameter rf_taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
m3 <- trans_diff$new(dataset = meta_meco_df, method = "rf", 
                     group = "Region", rf_taxa_level = "Genus")

# m3$res_rf is the result stored in the object
# plot the result
m4 <- m3$plot_diff_abund(use_number = 1:15, only_abund_plot = F)
m4$p1 <- m4$p1 + theme(axis.text = element_text(size = 12, colour = 'black')) +
  coord_flip()
m4$p2 <- m4$p2 + scale_fill_manual(values = c('#d95f02', '#1b9e77')) + 
  scale_color_manual(values = c('#d95f02', '#1b9e77')) +
  theme(legend.position = c(0.8, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12, colour = 'black'))
carbon_rf_plot <- gridExtra::grid.arrange(m4$p1, m4$p2, ncol=2, nrow = 1, widths = c(2,2)) # the middle asterisk represent the significances

carbon_plot <- ggdraw() +
  draw_plot(carbon_fun_plot, x = 0, y = 0.45, width = 1, height = 0.55) +
  draw_plot(carbon_rf_plot, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("(a)", "(b)"), size = 14,
                  x = c(0, 0), y = c(1, 0.45))

carbon_plot

#N cycling taxa analysis
# If you want to change the group list, reset the list t2$func_group_list
N_cycle_process <- c('aerobic_ammonia_oxidation', 'aerobic_nitrite_oxidation', 'nitrification',
                     'denitrification', 'nitrogen_fixation',
                     'nitrite_respiration',
                     'nitrate_respiration', 'nitrate_reduction', 'nitrogen_respiration')
m2$func_group_list$`N-cycle` <- N_cycle_process
# use show_prok_func to see the detailed information of prokaryotic traits
m2$show_prok_func("nitrification")

meta_n_cycle <- data.frame(m2$res_spe_func_perc[ ,colnames(m2$res_spe_func_perc) %in% c(m2$func_group_list$`N-cycle`)], 
                           m2$sample_table[, c('Site', 'Region')])

vars <- colnames(m2$res_spe_func_perc[ ,colnames(m2$res_spe_func_perc) %in% c(m2$func_group_list$`N-cycle`)])
library(lme4)
library(lmerTest)
mode <- lapply(vars, function(x) {
  lmer(substitute(i ~ Region + (1|Site), list(i = as.name(x))), data = meta_n_cycle)})
lapply(mode, summary)

library(plyr)
melt_df <- melt(meta_n_cycle[, !colnames(meta_n_cycle) %in% 'Site'], id.vars = c('Region'))
new_df <- ddply(melt_df, c('Region','variable'), summarise,
                mean = mean(value), sd = sd(value),
                sem = sd(value)/sqrt(length(value)))
sig <- c('b', 'a', 'b', 'b', 'b', 'b', 'a', 'a', 'a', 
         'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a')
new_df <- cbind(new_df ,sig)

#[1] "aerobic_ammonia_oxidation"    **                                      
#[2] "aerobic_nitrite_oxidation"                           
#[3] "nitrification"  **
#[4] "nitrate_denitrification"   **                        
#[5] "nitrite_denitrification"   **             
#[6] "nitrous_oxide_denitrification" **
#[7] "denitrification"   **                     
#[8] "nitrogen_fixation"   **                                     
#[9] "nitrate_ammonification"                                     
#[10] "nitrite_ammonification"                                         
#[11] "nitrite_respiration"  **                                          
#[12] "nitrate_respiration"                                        
#[13] "nitrate_reduction"                                         
#[14] "nitrogen_respiration"

#plot
nitrogen_fun_plot <- ggplot(new_df,aes(x = variable, y = mean, fill = Region))+
  geom_bar(position = 'dodge', stat = 'identity', width = 0.7)+
  scale_fill_manual(values=c('#d95f02', '#1b9e77'))+
  geom_errorbar(aes(ymin = mean-sem, ymax = mean + sem), width=.2,  position = position_dodge(0.7))+
  scale_x_discrete(limits = c(vars))+
  labs(x = 'Nitrogen cycle', y = 'Mean relative abundance (%)', fill = "Region")+
  scale_y_continuous(expand = c(0,0), limits = c(0,8))+
  geom_text(aes(label = sig, y = mean + sem + 0.02*max(mean)), 
            position = position_dodge(0.7), vjust = 0)+
  theme_bw()+
  theme(legend.position = c(0.7,0.88),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
nitrogen_fun_plot
#extract asvs involved in N cycling
meta_n_cycle_otus <- m2$res_spe_func[ ,c(m2$func_group_list$`N-cycle`)]
meta_n_cycle_otus <- rownames(meta_n_cycle_otus[rowSums(meta_n_cycle_otus) != 0,])

meta_N_phy <- subset_taxa(meta_physeq, OTU %in% meta_n_cycle_otus)
meta_N_phy <- prune_samples(sample_sums(meta_N_phy) > 0, meta_N_phy)
meta_N_phy_rel <- subset_taxa(meta_physeq_rel, OTU %in% meta_n_cycle_otus)
meta_N_phy_rel <- prune_samples(sample_sums(meta_N_phy_rel) > 0, meta_N_phy_rel)

meta_meco_df <- phyloseq2meco(meta_N_phy)
m5 <- meta_meco_df$cal_abund()
#trans_diff class, The third approach is rf, which depends on the random forest[14, 15] and the non-parametric test. 
# use Genus level for parameter rf_taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
m5 <- trans_diff$new(dataset = meta_meco_df, method = "rf", 
                     group = "Region", rf_taxa_level = "Genus")

# m3$res_rf is the result stored in the object
# plot the result
m6 <- m5$plot_diff_abund(use_number = 1:20, only_abund_plot = FALSE)
m6$p1 <- m6$p1 + theme(axis.text = element_text(size = 12, colour = 'black'))
m6$p2 <- m6$p2 + theme(legend.position = c(0.8, 0.2),
                       legend.title = element_text(size = 14),
                       legend.text = element_text(size = 12),
                       axis.text = element_text(size = 12, colour = 'black'))
nitrogen_rf_plot <- gridExtra::grid.arrange(m6$p1, m6$p2, ncol=2, nrow = 1, widths = c(2,2)) 
# the middle asterisk represent the significances

nitrogen_plot <- ggdraw() +
  draw_plot(nitrogen_fun_plot, x = 0, y = 0.45, width = 1, height = 0.55) +
  draw_plot(nitrogen_rf_plot, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("(a)", "(b)"), size = 14,
                  x = c(0, 0), y = c(1, 0.45))

# multivariables analysis for tibet plateau dataset
## read tibet plateau dataset
### read in otu table
otu.table <- read.csv('./tibet_dada2_asv/data/total_taxa_data/otu.table.csv',sep=",", row.names=1)
otu.table <- as.matrix(otu.table)

### read in taxonomy
taxonomy <- read.csv('./tibet_dada2_asv/data/total_taxa_data/taxonomy.csv',sep=",",row.names=1)
taxonomy <- as.matrix(taxonomy)

### read in metadata
metadata <- read.csv("./tibet_dada2_asv/data/metadata.csv", row.names=1, header = T)

### read in tree
total.tree <- read_tree('./tibet_dada2_asv/data/total_taxa_data/tree.nwk')

### read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "./tibet_dada2_asv/data/total_taxa_data/ref.seqs.fasta",
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

### import as phyloseq objects
otu.table <- otu_table(otu.table, taxa_are_rows = TRUE)
taxonomy <- tax_table(taxonomy)
metadata <- sample_data(metadata)

#### merge into one phyloseq object
tp_physeq <- phyloseq(otu.table, taxonomy, metadata, total.tree, ref_seqs)
tp_physeq
tp_physeq_rel <- microbiome::transform(tp_physeq, "compositional")

#diff class using microeco package
library(microeco)
meco_df <- phyloseq2meco(tp_physeq)
#calculate the abundance table
t1 <- meco_df$cal_abund()
# create object of trans_func
t2 <- trans_func$new(meco_df)
# mapping the taxonomy to the database
# the function can recognize prokaryotes or fungi automatically.
t2$cal_spe_func()
# return m2$res_spe_func, 1 represent function exists, 0 represent no or cannot confirmed.
t2$res_spe_func[1:5, 1:2]
t2$func_group_list
#carbon cycling taxa analysis
# If you want to change the group list, reset the list t2$func_group_list
# C_cycle_process <- c("methanotrophy",  "methanogenesis", 
#                     "methanol_oxidation", "methylotrophy", "chitinolysis",
#                     "cellulolysis", "xylanolysis", "ligninolysis", "fermentation",
#                     'aromatic_hydrocarbon_degradation','aromatic_compound_degradation',
#                     'hydrocarbon_degradation')
# t2$func_group_list$`C-cycle` <- C_cycle_process
# use show_prok_func to see the detailed information of prokaryotic traits
t2$show_prok_func("methanotrophy")
# calculate the percentages for communities
t2$cal_spe_func_perc(use_community = TRUE)
t2$res_spe_func_perc[1:5, 1:2]
tp_fun_table <- t2$res_spe_func_perc

#extract asvs involved in carbon cycling
tp_C_cycle <- t2$res_spe_func[ ,c(t2$func_group_list$`C-cycle`)]
tp_C_cycle_otus <- rownames(tp_C_cycle[rowSums(tp_C_cycle) != 0,])

tp_C_phylo <- subset_taxa(tp_physeq, ASV %in% tp_C_cycle_otus)
tp_C_phylo <- prune_samples(sample_sums(tp_C_phylo) > 0, tp_C_phylo)

# alpha diversity analysis
tp_total_alpha_div <- estimate_richness(tp_physeq, measures = c("Chao1", 'Shannon', 'Simpson'))
tp_carbon_alpha_div <- estimate_richness(tp_C_phylo, measures = c("Chao1", 'Shannon', 'Simpson'))

tp_alpha_div <- data.frame(group = c(rep('Total comunity', nrow(sample_data(tp_physeq))), 
                                     rep('Carbon cycling comunity', nrow(sample_data(tp_C_phylo)))), 
                           rbind(cbind(tp_total_alpha_div, sample_data(tp_physeq)[, c('Site', 'Sitegroup1', 'latitude', 'longitude')]),
                                 cbind(tp_carbon_alpha_div, sample_data(tp_C_phylo)[, c('Site', 'Sitegroup1', 'latitude', 'longitude')])))
## remove the se of Chao1 index
tp_alpha_div <- tp_alpha_div[,-3]
## rename the colnames
colnames(tp_alpha_div) <- c('group', "Chao1", 'Shannon', 'Simpson', 'Site', 'Sitegroup1', 'latitude', 'longitude')
tp_alpha_div[1:5, 1:5]

## calculate the average alpha diversity for each site
library(dplyr)
tp_alpha_div_plot_data <-  tp_alpha_div %>% 
  group_by(group, Site, Sitegroup1) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE)) %>%
  tidyr::gather(diversity, value, -c('group', 'Site', 'Sitegroup1', 'latitude', 'longitude'))
tp_alpha_div_plot_data$group <- factor(tp_alpha_div_plot_data$group, 
                                       levels = c('Total comunity', 'Carbon cycling comunity'))
### plot for Chao1
p1 <-  ggplot(tp_alpha_div_plot_data[tp_alpha_div_plot_data$diversity == 'Chao1',], aes(x = longitude, y = latitude, size = value)) +
  geom_point(alpha = 0.3, colour = '#1b9e77') +
  scale_size(range = c(.005, 10)) +
  coord_map('polyconic') +
  facet_grid(group ~ .) +
  xlab('Longitude (E°)')+ylab('Latitude (N°)')+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.title = element_text(colour = 'black',size=14),
        axis.ticks.length = unit(0.4,'lines'), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = 'black'), 
        axis.text = element_text(colour='black',size=12),
        strip.text = element_text(size = 12),
        legend.position = 'top',
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))+
  guides(size = guide_legend(nrow = 1))
### plot for Shannon
p2 <- p1 %+% tp_alpha_div_plot_data[tp_alpha_div_plot_data$diversity == 'Shannon',]
### plot for Simpson
p3 <- p1 %+% tp_alpha_div_plot_data[tp_alpha_div_plot_data$diversity == 'Simpson',]
library(gridExtra)
grid.arrange(p1, p2, p3, ncol = 3)

# beta diversity analysis
## local contribution to beta diversity (LCBD) analysis for total community
library(adespatial)
beta_tax_div <- beta.div(t(as.matrix(otu_table(tp_physeq))), 
                         method = "hellinger", sqrt.D = FALSE, samp = TRUE, 
                         nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)
env_div <- data.frame(LCBD = beta_tax_div$LCBD, sample_data(tp_physeq))
### calculate the average LCBD for each site
library(dplyr)
env_div_agg <-  env_div %>% 
  dplyr::select(-c(2:6, 9, 12)) %>%
  group_by(Site, Sitegroup1) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))
#write.csv(env_div_agg, file = './tibet_dada2_asv/results/tables/env_div_agg.csv')

### determine the relationships between LCBD and envs using linear regression modes
vars <- c("MAP", "MAT", "DOC", "S275_295", "SUVA254", "a300", "FI", "FrI", "HIX",
          "TN", "NH4_N", "DO", "pH", "Conductivity", "Salinity",
          "K", "Ca", "Na", "Mg")
mode <- lapply(vars, function(x) {
  lm(substitute(LCBD ~ i, list(i = as.name(x))), data = env_div_agg)})
sum.mode <- lapply(mode, broom::glance)
### normality test using Shapiro-Wilk test 
res <- lapply(mode, residuals)
norm_test <- lapply(res, shapiro.test)
norm_results <- data.frame(
  variables = vars, 
  w = sapply(norm_test, "[[", "statistic"), 
  pvalue = sapply(norm_test, "[[", "p.value")
)
norm_results

### extract the standardized regression coefficients
sd.coeff <- lapply(mode, QuantPsyc::lm.beta)
### arrange the table for plot
LCBD <- c(rep('LCBD', length(vars)))
sd.coeff <- sapply(sd.coeff, function(x){as.numeric(x[1])})
r.squared <- sapply(sum.mode, "[[", "r.squared")
adj.r.squared <- sapply(sum.mode, "[[", "adj.r.squared")
pvalue <- sapply(sum.mode, "[[", "p.value")
sig <- cut(pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
results <- data.frame(vars, LCBD, sd.coeff, r.squared, adj.r.squared, pvalue, sig)
results

#model selection
library(MASS)
library(glmulti)
A1 <- glmulti(LCBD ~ MAP + MAT + S275_295 + SUVA254 + a300 + FrI + HIX +
              TN + Conductivity + Salinity + Mg + K + Na, data=env_div_agg,
               level=1, fitfunction=lm, crit="aicc", confsetsize= 2^13, plotty = F, trace = 0)
top <- weightable(A1)
###  models with values more than 2 units away are considered substantially 
### less plausible than those with AICc values closer to that of the best model. 
### refrence:Anderson, D. R. (2007). Model based inference in the life sciences: A primer on evidence. New York: Springer. 
top_1 <- top[top$aicc <= min(top$aicc) + 2,] # 
top_1

modes_inf <- NULL
for(i in 1:nrow(top_1)){
  rse_sum <- summary(A1@objects[[i]])
  adj.r.squared <- rse_sum$adj.r.squared # obtain the adjust r squared
  multicollinearity <- any(sqrt(car::vif(A1@objects[[i]])) > 2) # check the multicollinearity
  tmp <- data.frame(adj.r.squared, multicollinearity)
  if(is.null(modes_inf)){
    modes_inf<-tmp
    } else {
      modes_inf <- rbind(modes_inf,tmp)
    } 
}
modes_inf <- cbind(top_1, modes_inf)
modes_inf

vpa.mod <- varpart(env_div_agg$LCBD, ~ env_div_agg$HIX,
                   ~ env_div_agg$MAP)
plot(vpa.mod)

## local contribution to beta diversity (LCBD) analysis for carbon cycling community
beta_c_div <- beta.div(t(as.matrix(otu_table(tp_C_phylo))), 
                       method = "hellinger", sqrt.D = FALSE, samp = TRUE, 
                       nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)

env_div <- data.frame(LCBD = beta_c_div$LCBD, sample_data(tp_physeq))
LCBD <- data.frame(total_LCBD = beta_tax_div$LCBD, carbon_LCBD = beta_c_div$LCBD)
# write.csv(LCBD, file = 'E:/thermokast_lakes/water_microbes/tibet_dada2_asv/results/tables/LCBD.CSV')

### calculate the average LCBD for each site
library(dplyr)
env_div_agg_carbon <-  env_div %>% 
  dplyr::select(-c(2:6, 9, 12)) %>%
  group_by(Site, Sitegroup1) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))
#write.csv(env_div_agg_carbon, file = './tibet_dada2_asv/results/tables/env_div_agg_carbon.csv')

### determine the relationships between LCBD and envs using linear regression modes
vars <- c("MAP", "MAT", "DOC", "S275_295", "SUVA254", "a300", "FI", "FrI", "HIX",
          "TN", "NH4_N", "DO", "pH", "Conductivity", "Salinity",
          "K", "Ca", "Na", "Mg")
mode <- lapply(vars, function(x) {
  lm(substitute(LCBD ~ i, list(i = as.name(x))), data = env_div_agg_carbon)})
sum.mode <- lapply(mode, broom::glance)
### normality test using Shapiro-Wilk test 
res <- lapply(mode, residuals)
norm_test <- lapply(res, shapiro.test)
norm_results_carbon <- data.frame(
  variables = vars, 
  w = sapply(norm_test, "[[", "statistic"), 
  pvalue = sapply(norm_test, "[[", "p.value")
)
norm_results_carbon
env_div_agg_carbon$DOC <- log(env_div_agg_carbon$DOC)
env_div_agg_carbon$S275_295 <- log(env_div_agg_carbon$S275_295)
env_div_agg_carbon$a300 <- log(env_div_agg_carbon$a300)
env_div_agg_carbon$FI <- log(env_div_agg_carbon$FI)
env_div_agg_carbon$TN <- log(env_div_agg_carbon$TN)
env_div_agg_carbon$NH4_N <- log(env_div_agg_carbon$NH4_N)
env_div_agg_carbon$Depth <- log(env_div_agg_carbon$Depth)
env_div_agg_carbon$DO <- log(env_div_agg_carbon$DO)
#env_div_agg_carbon$pH <- log(env_div_agg_carbon$pH)
env_div_agg_carbon$SUVA254 <- log(env_div_agg_carbon$SUVA254)
env_div_agg_carbon$Conductivity <- log(env_div_agg_carbon$Conductivity)
env_div_agg_carbon$Salinity <- log(env_div_agg_carbon$Salinity)
env_div_agg_carbon$Ca <- log(env_div_agg_carbon$Ca)
env_div_agg_carbon$Mg <- log(env_div_agg_carbon$Mg)
env_div_agg_carbon$K <- log(env_div_agg_carbon$K)
env_div_agg_carbon$Na <- log(env_div_agg_carbon$Na)
env_div_agg_carbon$ALT <- log(env_div_agg_carbon$ALT)

### normality test again using Shapiro-Wilk test 
mode <- lapply(vars, function(x) {
  lm(substitute(LCBD ~ i, list(i = as.name(x))), data = env_div_agg_carbon)})
sum.mode <- lapply(mode, broom::glance)
### normality test using Shapiro-Wilk test 
res <- lapply(mode, residuals)
norm_test <- lapply(res, shapiro.test)
norm_results_carbon <- data.frame(
  variables = vars, 
  w = sapply(norm_test, "[[", "statistic"), 
  pvalue = sapply(norm_test, "[[", "p.value")
)
norm_results_carbon
### extract the standardized regression coefficients
sd.coeff <- lapply(mode, QuantPsyc::lm.beta)
### arrange the table for plot
LCBD <- c(rep('LCBD', length(vars)))
sd.coeff <- sapply(sd.coeff, function(x){as.numeric(x[1])})
r.squared <- sapply(sum.mode, "[[", "r.squared")
adj.r.squared <- sapply(sum.mode, "[[", "adj.r.squared")
pvalue <- sapply(sum.mode, "[[", "p.value")
sig <- cut(pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
results_carbon <- data.frame(vars, LCBD, sd.coeff, r.squared, adj.r.squared, pvalue, sig)
results_carbon

## models selection
library(MASS)
library(glmulti)
A2 <- glmulti(LCBD ~ MAP + MAT + S275_295 + SUVA254 + a300 + FrI + HIX +
                TN + NH4_N + pH + Conductivity + Salinity + Mg + K + Na, data=env_div_agg_carbon,
              level=1, fitfunction=lm, crit="aicc", confsetsize= 2^15, plotty = F, trace = 0)
top <- weightable(A2)
###  models with values more than 2 units away are considered substantially 
### less plausible than those with AICc values closer to that of the best model. 
### refrence:Anderson, D. R. (2007). Model based inference in the life sciences: A primer on evidence. New York: Springer. 
### top_2 <- top[top$aicc <= min(top$aicc) + 2,] #
top_2 <- top[1:30,]
top_2

modes_inf_carbon <- NULL
for(i in 1:30){
  rse_sum <- summary(A2@objects[[i]])
  adj.r.squared <- rse_sum$adj.r.squared # obtain the adjust r squared
  multicollinearity <- any(sqrt(car::vif(A2@objects[[i]])) > 2) # check the multicollinearity
  tmp <- data.frame(adj.r.squared, multicollinearity)
  if(is.null(modes_inf_carbon)){
    modes_inf_carbon <- tmp
  } else {
    modes_inf_carbon <- rbind(modes_inf_carbon,tmp)
  } 
}
modes_inf_carbon <- cbind(top_2, modes_inf_carbon)
modes_inf_carbon

vpa.mod_carbon <- varpart(env_div_agg_carbon$LCBD, ~ env_div_agg_carbon$SUVA254,
                   ~ env_div_agg_carbon$MAP)
plot(vpa.mod_carbon)

## plot using ggplot package
### maping the LCBD among the study area
#### total community
geo_div <- ggplot(data = env_div_agg, aes(x = longitude, y = latitude, size = LCBD)) +
  geom_point(aes(alpha = 0.2), colour = '#1b9e77') +
  scale_size(range = c(.005, 10)) +
  coord_map('polyconic') +
  xlab('Longitude')+ylab('Latitude')+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.title = element_text(colour = 'black',size=14),
        axis.ticks.length = unit(0.4,'lines'), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = 'black'), 
        axis.text = element_text(colour='black',size=12),
        legend.position = 'top',
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

#### carbon cycling community
geo_div_carbon <- ggplot(data = env_div_agg_carbon, aes(x = longitude, y = latitude, size = LCBD)) +
  geom_point(aes(alpha = 0.2), colour = '#d95f02') +
  scale_size(range = c(.005, 10)) +
  coord_map("polyconic") +
  xlab('Longitude')+ylab('Latitude')+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.title = element_text(colour = 'black',size=14),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.position = 'top',
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

geo_LCBD_plot <- plot_grid(geo_div, geo_div_carbon, 
                  labels = c('(a)', '(b)'), ncol = 2, 
                  label_x = .01, label_y = 1, 
                  hjust = 0, label_size = 14, align = "v")
geo_LCBD_plot

### heatmap using standardized regression coefficients to explore the relationship between the LCBD and environment factors
results_plot_data <- data.frame(group = c(rep('Total community', nrow(results)), rep('Carbon cycling', nrow(results_carbon))),
                                rbind(results, results_carbon))

results_plot_data$vars <- factor(results_plot_data$vars,levels = rev(vars))
results_plot_data$group <- factor(results_plot_data$group,
                                  levels=c('Total community', 'Carbon cycling'))
p_env_div <- ggplot(aes(x=LCBD, y=vars, fill=sd.coeff), data=results_plot_data) +
  geom_tile() +
  scale_fill_gradient2(low='#1b9e77', mid='white', high='#d95f02') +
  geom_text(aes(label=sig), color="black", size=6) +
  labs(y=NULL, x=NULL, fill='Standardized regression coefficients') +
  facet_wrap( .~ group, ncol = 2) +
  theme_bw()+
  theme(legend.position="bottom", 
        panel.border = element_blank(),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
###  plot the linear regression relationships between the LCBD and best explained variables
#### total community
p_linear <- env_div_agg %>%
  dplyr::select(LCBD, MAP, HIX) %>%
  tidyr::gather(varibales, value, MAP:HIX, factor_key=TRUE) %>%
  ggplot(aes(value, LCBD)) +
  geom_point(size=3.5, alpha=0.8, aes(colour = as.factor(varibales))) +
  geom_smooth(method="lm", size=1, se=T, colour='black') +
  scale_color_manual(values = c('#1b9e77', '#d95f02')) +
  scale_y_continuous(limits = c(0, 0.01)) +
  facet_wrap( .~ varibales, scales="free_x", ncol = 2) +
  ylab('LCBD')+xlab('Values') +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        strip.text = element_text(size = 14),
        legend.position='none')
#### carbon cycling community
plot_carbon_data <- env_div_agg_carbon %>%
  dplyr::select(LCBD, MAP, SUVA254) %>%
  tidyr::gather(varibales, value, MAP:SUVA254, factor_key=TRUE)
p_carbon_linear <- ggplot(plot_carbon_data, aes(value, LCBD)) +
  geom_point(size=3.5, alpha=0.8, aes(colour = as.factor(varibales))) +
  geom_smooth(method="lm", size=1, se=T, colour='black') +
  scale_color_manual(values = c( '#1b9e77', '#d95f02')) +
  scale_y_continuous(limits = c(0, 0.01)) +
  facet_wrap( .~ varibales, scales="free_x", ncol = 2) +
  ylab('LCBD')+xlab('Values') +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        strip.text = element_text(size = 14),
        legend.position='none')

## venn plot
library("VennDiagram")
venn.plot <- draw.pairwise.venn(area1 = 0.27, area2 = 0.31, cross.area = 0.21,
                                category=c('DOM proporties', 'Climate elements'),
                                fill = c( '#d95f02', '#1b9e77'), scaled = 0, 
                                ind = FALSE, cat.col = c(rep('black',2)), 
                                cat.cex = 1.2, cat.dist =  c(0.02, 0.02),
                                cat.pos = c(30, 330), margin = 0.05, lty = 'blank')
venn.plot_carbon <- draw.pairwise.venn(area1 = 0.51, area2 = 0.35, cross.area = 0.31,
                                       category=c('DOM proporties', 'Climate elements'),
                                       fill = c( '#d95f02', '#1b9e77'), scaled = 0, 
                                       ind = FALSE, rotation.degree = 180, cat.col = c(rep('black',2)), 
                                       cat.cex = 1.2, cat.dist =  c(0.02, 0.02),
                                       cat.pos = c(30, 330), margin = 0.05, lty = 'blank')
venn_LCBD_plot <- plot_grid(venn.plot, venn.plot_carbon, 
                           labels = c('(a)', '(b)'), ncol = 2, 
                           label_x = .01, label_y = 1, 
                           hjust = 0, label_size = 14, align = "v")
venn_LCBD_plot


vpa_dat <- data.frame(c(rep('Total community', 3), rep('Carbon cycling community', 3)),
                      c(rep(c('Climate factors', 'Overlap', 'DOM properties'), 2)), c(10, 21, 6, 4, 31, 20))
colnames(vpa_dat) <- c('group', 'fraction', 'value')
vpa_dat$group <- factor(vpa_dat$group, levels=c('Total community', 'Carbon cycling community'))
library(ggplot2)                                                                   
p_varpart <- ggplot(vpa_dat, aes(fill=fraction, y=value, x=group)) + 
  geom_bar(position="stack", stat="identity", alpha=0.8) +
  scale_fill_manual(values = c('#1b9e77', '#d95f02', "grey")) +
  scale_y_continuous(expand = c(0, 0))+
  ylab("Relative contribution (%)")+
  theme_bw()+
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 13),
        axis.text = element_text(size=12), 
        panel.grid = element_blank()) +
  guides(fill = guide_legend(ncol = 1))
ggdraw() +
  draw_plot(p_env_div, x = 0, y = 0, width = 0.25, height = 1) +
  draw_plot(p_linear, x = 0.25, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(p_carbon_linear, x = 0.25, y = 0, width = 0.5, height = 0.5) +
  draw_plot(p_varpart, x = 0.75, y = 0, width = 0.25, height = 1) +
  draw_plot_label(label = c("(a)", "(b)", "(c)", '(d)'), size = 12,
                  x = c(0, 0.25, 0.25, 0.75), y = c(1, 1, 0.5, 1))

print(p_varpart)

## (partial) mantel test for each variable
library(ggcor)
## mantel test function
mantel.fun <- function(phylo) {
  env.table <- data.frame(sample_data(phylo))
  otu_table <- as.matrix(t(otu_table(phylo)))
  otu_table_hel <- decostand(otu_table, 'hellinger')
  otu_table_hel_dist <- vegdist(otu_table_hel, 'bray',upper=F)
  df <- NULL
  vars <- c("DOC", "TN", "NH4_N", "S275_295", "SUVA254",
            "a300", "MAP", "MAT", "DO", "pH", "Conductivity", 
            "Salinity", "Ca", "Mg", "K", "Na", "FI", "FrI", "BIX", "HIX")
  for (x in vars) {
    x.dist <- vegdist(scale(env.table[,x]), 'euclidean')
    mode <- mantel(x.dist, otu_table_hel_dist, 
                   method = "pearson", permutations = 999)
    r <- mode$statistic
    p <- mode$signif
    tmp <- data.frame(variable = x, r = r, p.value = p)
    if(is.null(df))
      df <- tmp
    else
      df <- rbind(df ,tmp)
  }
  return(df)
}
## partial mantel test function
partial.mantel.fun <- function(phylo) {
  env.table <- data.frame(sample_data(phylo))
  otu_table <- as.matrix(t(otu_table(phylo)))
  otu_table_hel <- decostand(otu_table, 'hellinger')
  otu_table_hel_dist <- vegdist(otu_table_hel, 'bray',upper=F)
  df <- NULL
  vars <- c("DOC", "TN", "NH4_N", "S275_295", "SUVA254",
            "a300", "MAP", "MAT", "DO", "pH", "Conductivity", 
            "Salinity", "Ca", "Mg", "K", "Na","FI", "FrI", "HIX")
  for (x in vars) {
    x.dist <- vegdist(scale(env.table[,x]), 'euclidean')
    z.dist <- vegdist(scale(env.table[ , setdiff(vars, x)]), 'euclidean')
    mode <- mantel.partial(x.dist, otu_table_hel_dist, z.dist, 
                           method = "pearson", permutations = 999)
    r <- mode$statistic
    p <- mode$signif
    tmp <- data.frame(env = x, r = r, p.value = p)
    if(is.null(df))
      df <- tmp
    else
      df <- rbind(df ,tmp)
  }
  return(df)
}

envfit.fun <- function(phylo){
  env.table <- data.frame(sample_data(phylo))
  vars <- c("DOC", "TN", "NH4_N", "S275_295", "SUVA254",
            "a300", "MAP", "MAT", "DO", "pH", "Conductivity", 
            "Salinity", "Ca", "Mg", "K", "Na","FI", "FrI", "HIX")
  env.table <- env.table[ , vars]
  otu_table <- as.matrix(t(otu_table(phylo)))
  otu_table_hel <- decostand(otu_table, 'hellinger')
  otu_table_hel_dist <- vegdist(otu_table_hel, 'bray',upper=F)
  ord <-  cmdscale(otu_table_hel_dist,  k = 2, eig = T, add = T)
  fit <- envfit(ord, env.table, perm=999, na.rm = TRUE)
  fit
}

total.commun.par.mant <- partial.mantel.fun(tp_physeq)
carbon.commun.par.mant <- partial.mantel.fun(tp_C_phylo)
## PLOT
## devtools::install_github('hannet91/ggcor')
library(ggcor)
set.seed(123456)

par.man.tibble <- tibble(spec = c(rep('Total community composition', nrow(total.commun.par.mant)), 
                                  rep('Carbon community composition', nrow(carbon.commun.par.mant))), 
                         rbind(total.commun.par.mant, carbon.commun.par.mant))
#par.man.tibble <- tibble(spec = c(rep('Total community composition', nrow(total.commun.par.mant))), 
#                         total.commun.par.mant)
vars <- c("MAT", "MAP", "DOC", "S275_295", "SUVA254", "a300", "FI", "FrI", "HIX",
          "TN", "NH4_N", "pH", "DO", "Conductivity", "Salinity",
          "K", "Ca", "Na", "Mg")
env.table <- sample_data(tp_physeq)[ , vars]

mantel02 <- par.man.tibble %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.2, 0.5, Inf), 
                 labels = c("<0.20", "0.20-0.5", ">0.50"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"),
                       right = T))
quickcor(env.table, type = "upper") + geom_square() + 
  add_link(mantel02, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE) +
  scale_color_manual(values = c('#d95f02', '#1b9e77', '#3C5488FF', 'grey')) +
  scale_size_manual(values = c(0.5, 2)) +
  geom_diag_label() + remove_axis("x")

# lakeshore analysis with DDRs
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
lakeshores_phylo_even = rarefy_even_depth(lakeshores_phylo, sample.size = 33859, rngseed = 1234, replace = TRUE)


## DDRs for water and lakeshores
#colnames(geo_dis)<-c('site','lon','lat')
deg2rad <- function(deg) return(deg*pi/180)

# Calculates the geodesic distance between two points specified by 
# radian latitude/longitude using the Haversine formula
# Ouputs distance between sites 1 and 2 as meters
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = (R * c)*1000
  return(d) # Distance in meters
}


# Fxn to calculate matrix of distances between each two sites
# INPUT: a data frame in which longs are in first column and lats in second column
# OUTPUT: a distance matrix (class dist) between all pairwise sites
# Output distances are in meters
CalcDists <- function(longlats) {
  name <- list(rownames(longlats), rownames(longlats))
  n <- nrow(longlats)
  z <- matrix(0, n, n, dimnames = name)
  for (i in 1:n) {
    for (j in 1:n) z[i, j] <- gcd.hf(long1 = deg2rad(longlats[i, 1]), 
                                     lat1 = deg2rad(longlats[i, 2]), long2 = deg2rad(longlats[j, 1]), 
                                     lat2 = deg2rad(longlats[j, 2]))
  }
  z <- as.dist(z)
  return(z)
}

summary.model <- function(model){
  r.squre <- round(summary(model)$r.squared,3)
  p.value <- round(anova(model)$'Pr(>F)'[1],3)
  p.stars <- function(p.values) {
    unclass(symnum(p.values, corr = FALSE, 
                   na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")))}
  sig <- p.stars(p.value)
  slope <- round(model$coefficients[[2]],3)
  results<-data.frame(slope, r.squre, p.value, sig)
  return(results)
}

dist_table_cal <- function(phylo) {
  env <- sample_data(phylo)
  #geofraphical distance
  geo_dist <- data.frame(env[,c('longitude','latitude')])
  geo.dist <- CalcDists(geo_dist)/100000
  Geography.dist <- as.numeric(geo.dist)
  
  ##calculate the disimilarity matrix of bacterial community based on bray-kurties distance
  otu_table <- as.matrix(t(otu_table(phylo)))
  otu_table_hel <- decostand(otu_table, 'hellinger')
  otu_table_hel_dist <- vegdist(otu_table_hel, 'bray',upper=F)
  
  com_dissimilarty <- as.numeric(otu_table_hel_dist)*100
  com_similarty <- (1 - as.numeric(otu_table_hel_dist))*100
  #arrange the distance table
  dist.table <- as.data.frame(cbind(com_similarty, com_dissimilarty, Geography.dist))
  colnames(dist.table)<-c('com_similarty', 'com_dissimilarty', 'Geography')
  return(list(geo.dist, otu_table_hel_dist, dist.table))
}

dist.table <- dist_table_cal(tp_phylo)[[3]]
dist.table_shores <- dist_table_cal(lakeshores_phylo_even)[[3]]
dist.table_tp <- dist_table_cal(tp_physeq)[[3]]

fit <- lm(com_similarty ~ Geography, data = dist.table)
fit
summary(fit)

fit.shores <- lm(com_similarty ~ Geography, data = dist.table_shores)
fit.shores
summary(fit.shores)

fit.tp <- lm(com_similarty ~ Geography, data = dist.table_tp)
fit.tp
summary(fit.tp)

#partial mantel test
library(vegan)
mantel(dist_table_cal(tp_phylo)[[1]], dist_table_cal(tp_phylo)[[2]], method="pearson")
mantel(dist_table_cal(tp_physeq)[[1]], dist_table_cal(tp_physeq)[[2]], method="pearson")
mantel(dist_table_cal(lakeshores_phylo_even)[[1]], 
       dist_table_cal(lakeshores_phylo_even)[[2]], method="pearson")

##plot
mytheme <- theme(panel.grid = element_blank(),
                 panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA, colour = "black"), 
                 axis.title = element_text(colour = 'black',size=14),
                 axis.ticks.length = unit(0.4,"lines"), 
                 axis.ticks = element_line(color='black'),
                 axis.line = element_line(colour = "black"), 
                 axis.text = element_text(colour='black',size=12),
                 legend.position = c(0.85,0.82),
                 legend.title=element_text(size = 12),
                 legend.text=element_text(size=9),
                 legend.key=element_blank(),
                 legend.background = element_rect(colour = "white"))

p_distance <- ggplot(dist.table_tp, aes(Geography, com_similarty)) + 
  geom_point(shape = 19, alpha = 0.1, colour = '#1b9e77')+
  geom_smooth(method = "lm", size = 1, se = T, colour = 'black') +
  ylab('Community similarity (%)')+xlab('Geographic distance (100km)') +
  mytheme

p_distance_shores <- ggplot(dist.table_shores, aes(Geography, com_similarty)) + 
  geom_point(shape = 19, alpha = 0.1, colour = '#d95f02')+
  geom_smooth(method = "lm", size = 1, se = T, colour = 'black') +
  ylab('Community similarity (%)')+xlab('Geographic distance (100km)') +
  mytheme


#determine the slopes of DDRs among three region:QL, QZR, and MD
tp_ql <- subset_samples(tp_physeq, Sitegroup1 == 'QL')
tp_qzr <- subset_samples(tp_physeq, Sitegroup1 == 'QZR')
tp_md <- subset_samples(tp_physeq, Sitegroup1 == 'MD')

dist.table.ql <- dist_table_cal(tp_ql)[[3]]
dist.table.qzr <- dist_table_cal(tp_qzr)[[3]]
dist.table.md <- dist_table_cal(tp_md)[[3]]

fit1 <- lm(com_similarty ~ Geography, data = dist.table.ql)
fit1
summary(fit1)
ln.fit1 <- lm(com_similarty/100 ~ log(Geography*100+1), data = dist.table.ql)
ln.fit1
summary(ln.fit1)

fit2 <- lm(com_similarty ~ Geography, data = dist.table.qzr)
fit2
summary(fit2)
ln.fit2 <- lm(com_similarty/100 ~ log(Geography*100+1), data = dist.table.qzr)
ln.fit2
summary(ln.fit2)

fit3 <- lm(com_similarty ~ Geography, data = dist.table.md)
fit3
summary(fit3)
ln.fit3 <- lm(com_similarty/100 ~ log(Geography*100+1), data = dist.table.md)
ln.fit3
summary(ln.fit3)


#partial mantel test
library(vegan)
###geo dist matrix
mantel(geo.dist, otu_table_hel_dist, method = "pearson")

env_tp <- sample_data(tp_phylo)
extract_data <- env_tp[,c('Sitegroup1', 'MAP')]
melted <- melt(extract_data, id.vars = c("Sitegroup1"))
library(plyr)
MAP_CAL <- ddply(melted, c("Sitegroup1"), summarise,
                 mean = mean(value), sd = sd(value),
                 sem = sd(value)/sqrt(length(value)))

model <- lm(MAP ~ Sitegroup1, data = data.frame(env_tp))
summary.model(model)
summary(model)
TukeyHSD(aov(model))
sig <- c('a', 'ab', 'b')
MAP_CAL <- cbind(MAP_CAL, sig)


#PLOT
p_map <- ggplot(MAP_CAL, aes(x = Sitegroup1, y = mean, fill = Sitegroup1)) +
  geom_bar(stat = 'identity', color = 'black', width = 0.45) +
  geom_errorbar(aes(ymin = mean, ymax = mean + sem), width = .2) +
  geom_text(aes(label = sig, y = (mean + sem)*1.02), 
            position = position_dodge(0.9),vjust = 0) +
  scale_x_discrete(limits = c('MD', 'QL', 'QZR')) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 600)) +
  xlab('Region') + ylab('MAP (mm)') +
  mytheme

dist.table.3 <- rbind(dist.table.ql, dist.table.qzr, dist.table.md)
dist.table.3 <- cbind(region = c(rep('QL', nrow(dist.table.ql)), rep('QZR', nrow(dist.table.qzr)),
                                 rep('MD', nrow(dist.table.md))) ,dist.table.3)

p_3_DDRs <- ggplot(dist.table.3, aes(x = Geography, y = com_similarty, fill = region)) + 
  geom_point(size = 1.5, alpha = 0.3, aes(colour = region))+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method = lm, level = 0.95, aes(colour = region))+
  #stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  #theme_bw()+
  ylab('Community similarity (%)')+xlab('Geographic distance (100km)')+
  mytheme

#fig3
library(cowplot)
DDRs <- plot_grid(p_distance, p_distance_shores, p_map, p_3_DDRs, 
                  labels = c('(a)', '(b)', '(c)', '(d)'), ncol = 2, 
                  label_x = .01, label_y = 1, 
                  hjust = 0, label_size = 14, align = "v")
DDRs

##comparison of DDRs slope
m.interaction <- lm(com_similarty ~ Geography*taxa, data = dat)
anova(m.interaction)

# Obtain slopes
m.interaction$coefficients
geo.lst <- emtrends(m.interaction, "taxa", var="Geography")

# Compare slopes
pairs(geo.lst)

#log transformation
dat_log<-data.frame(taxa=dat$taxa,com_similarty=dat$com_similarty/100,Geography=log(dat$Geography +1))
m_log.interaction <- lm(com_similarty ~ Geography*taxa, data = dat_log)
anova(m_log.interaction)

m_log.interaction$coefficients
geo.log.lst <- emtrends(m_log.interaction, "taxa", var="Geography")

pairs(geo.log.lst)

























### 三因素
draw.triple.venn(area1 = 0.26, 
                 area2 = 0.30, 
                 area3 = 0.23, 
                 n12 = 0.21,
                 n13 = 0.19, 
                 n23 = 0.18, 
                 n123 = 0.17, 
                 fill = c("#56B4E9", "#009E73", '#D55E00'),
                 cat.col = c(rep('black',3)), cat.cex = 1.2,
                 margin = 0.05, lty = 'blank',
                 category = c("DOM proporties", "Climate elements",
                              "physio-chemical factors"))
plot_data <- env_div_agg %>%
  dplyr::select(LCBD, MAP, K, HIX, S275_295) %>%
  tidyr::gather(varibales, value, MAP:S275_295, factor_key=TRUE)

p_linear <- ggplot(plot_data, aes(value, LCBD)) + 
  geom_point(colour = "#1b9e77", size=3.5, alpha=0.8)+
  geom_smooth(method="lm", size=1, se=T, colour='black') +
  scale_y_continuous(limits = c(0, 0.01)) +
  facet_wrap( .~ varibales , scales="free_x", ncol = 2) +
  ylab('LCBD')+xlab('Varibales') +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        strip.text = element_text(size = 14),
        legend.position='none')
p_linear

p_map_lcbd <- ggplot(env_div_agg, aes(MAP, LCBD)) + 
  geom_point(shape = 19, colour = "#1b9e77", size=3.5, alpha=0.8)+
  geom_smooth(method = "lm", size = 1.5, se = T, colour = 'black') +
  scale_y_continuous(limits = c(0, 0.01)) +
  ylab('LCBD')+xlab('MAP (mm)') +
  mytheme
p_suva254_lcbd <- ggplot(env_div_agg, aes(SUVA254, LCBD)) + 
  geom_point(shape = 19, colour = "#1b9e77", size=3.5, alpha=0.8)+
  geom_smooth(method = "lm", size = 1.5, se = T, colour = 'black') +
  scale_y_continuous(limits = c(0, 0.01)) +
  ylab('LCBD')+xlab('SUVA254') +
  mytheme
p_hix_lcbd <- ggplot(env_div_agg, aes(HIX, LCBD)) + 
  geom_point(shape = 19, colour = "#1b9e77", size=3.5, alpha=0.8)+
  geom_smooth(method = "lm", size = 1.5, se = T, colour = 'black') +
  scale_y_continuous(limits = c(0, 0.01)) +
  ylab('LCBD')+xlab('HIX') +
  mytheme
p_s275_295_lcbd <- ggplot(env_div_agg, aes(S275_295, LCBD)) + 
  geom_point(shape = 19, colour = "#1b9e77", size=3.5, alpha=0.8)+
  geom_smooth(method = "lm", size = 1.5, se = T, colour = 'black') +
  scale_y_continuous(limits = c(0, 0.01)) +
  ylab('LCBD')+xlab('S275_295') +
  mytheme

ggdraw() +
  draw_plot(p_env_div, x = 0, y = 0, width = 0.2, height = 1) +
  draw_plot(p_linear, x = 0.2, y = 0, width = 0.80, height = 1) +
  draw_plot_label(label = c("(a)", "(b)"), size = 12,
                  x = c(0, 0.2), y = c(1, 1))


# SEM
## load the package
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(nlme)

climate_vars <- env_div_agg[ , c('MAP')]
DOM_vars <- env_div_agg[ , c('FrI')]
physio_chemical_vars <- env_div_agg[ , c('Salinity')]
ALT <- env_div_agg[ , c('ALT')]

library(piecewiseSEM)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI + Salinity +ALT, data = env_div_agg),
  lm(ALT ~ MAP, data = env_div_agg),
  lm(FrI ~ MAP + ALT, data = env_div_agg),
  lm(Salinity ~ MAP + ALT + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI + Salinity +ALT, data = env_div_agg),
  lm(ALT ~ MAP, data = env_div_agg),
  lm(FrI ~ MAP + ALT, data = env_div_agg),
  lm(Salinity ~ MAP + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI + Salinity +ALT, data = env_div_agg),
  lm(FrI ~ MAP + ALT, data = env_div_agg),
  lm(Salinity ~ MAP + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI + Salinity, data = env_div_agg),
  lm(FrI ~ MAP + ALT, data = env_div_agg),
  lm(Salinity ~ MAP + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI, data = env_div_agg),
  lm(FrI ~ MAP + ALT, data = env_div_agg),
  lm(Salinity ~ MAP + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI, data = env_div_agg),
  lm(FrI ~ MAP, data = env_div_agg),
  lm(Salinity ~ MAP + FrI, data = env_div_agg))
summary(sem_mod)
sem_mod <- psem(
  lm(LCBD ~ MAP + FrI, data = env_div_agg),
  lm(FrI ~ MAP, data = env_div_agg),
  lm(Salinity ~ MAP, data = env_div_agg))
summary(sem_mod)
fisherC(sem_mod)
rsquared(sem_mod)

coeffs_sem <- coefs(sem_mod, standardize = "scale")[,c(1,2,8)]
coeffs_sem
## Direct effect of climate variables on total community composition.
clim_dir_effect <- coeffs_sem[1,3]
## Ind effect of climate variables on total community composition.
clim_ind_effect <- coeffs_sem[3,3]*coeffs_sem[2,3]

total_clim_effect <- clim_dir_effect + clim_ind_effect

## Direct effect of ALT on total community composition.
ALT_dir_effect <- 0
## Ind effect of ALT on total community composition.
ALT_ind_effect <- 0
total_ALT_effect <- ALT_dir_effect + ALT_ind_effect

## Direct effect of DOM variables on total community composition.
DOM_dir_effect <- coeffs_sem[2,3]
## Ind effect of DOM variables on total community composition.
DOM_ind_effect <- 0
total_DOM_effect <- DOM_dir_effect + DOM_ind_effect

## Direct effect of physio_chemical variables on total community composition.
physio_chemical_dir_effect <- 0
## Ind effect of climate variables on total community composition.
physio_chemical_ind_effect <- 0
total_physio_chemical_effect <- physio_chemical_dir_effect + physio_chemical_ind_effect
## bar plot
sem_effect <- data.frame(variable = c(rep('Climate elements',3),
                                            rep('ALT', 3),
                                            rep('DOM properties', 3),
                                            rep('physio-chemical factors', 3)),
                               type_of_effects = c(rep(c('Total effects', 'Direct effects', 
                                                         'Indirect effects'), 4)),
                               value = c(total_clim_effect, clim_dir_effect, clim_ind_effect,
                                         total_ALT_effect, ALT_dir_effect, ALT_ind_effect,
                                         total_DOM_effect, DOM_dir_effect, DOM_ind_effect,
                                         total_physio_chemical_effect, physio_chemical_dir_effect,
                                         physio_chemical_ind_effect))
sem_effect$type_of_effects <- factor(sem_effect$type_of_effects, ordered = T,
                                           levels = c('Total effects', 'Direct effects', 'Indirect effects'))
sem_effect$variable <- factor(sem_effect$variable, ordered = T,
                                    levels = c('Climate elements', 'DOM properties', 'physio-chemical factors', 'ALT'))
sem_effect_plot <-  ggplot(data = sem_effect, 
                                 aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat="identity", width = 0.6) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-0.8, 0.5)) +
  labs(x = 'Variables', y = 'Standardized effects from SEM') +
  facet_wrap(~ type_of_effects) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(color='black',size=14),
        strip.text = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

sem_effect_plot


























# VPA
## total community
climate_vars <- env.table[ , c('MAP', 'MAT')]
DOM_vars <- env.table[ , c('SUVA254', 'S275_295')]
physio_chemical_vars <- env.table[ , c("pH", "K")]
ALT <- env.table[ , c('ALT')]
Region1 <- c(rep(1, 30), rep(2, 100), rep(3, 58))

## pca for climate
climate_PC1 <- princomp(climate_vars,cor=TRUE)
summary(climate_PC1, loadings=TRUE)
climate_PC1 <- predict(climate_PC1)
climate_PC1 <- climate_PC1[,1]

## pca for DOM
DOM_PC1 <- princomp(DOM_vars,cor=TRUE)
summary(DOM_PC1, loadings=TRUE)
DOM_PC1 <- predict(DOM_PC1)
DOM_PC1 <- DOM_PC1[,1]

## pca for physio_chemical factors
physio_chemical_PC1 <- princomp(physio_chemical_vars,cor=TRUE)
summary(physio_chemical_PC1, loadings=TRUE)
physio_chemical_PC1 <- predict(physio_chemical_PC1)
physio_chemical_PC1 <- physio_chemical_PC1[,1]

## pca for ALT
ALT_PC1 <- princomp(ALT,cor=TRUE)
summary(ALT_PC1, loadings=TRUE)
ALT_PC1 <- predict(ALT_PC1)
ALT_PC1 <- ALT_PC1[,1]

## site
Site <- env_tp$Site
Region1 <- env_tp$Region1
#ordination for community
otu_table <- as.matrix(t(otu_table(tp_phylo)))
otu_table_hel <- decostand(otu_table, method = "hellinger")
otu_table_hel_dist<-vegdist(otu_table_hel, 'bray',upper=F)
ord <-  metaMDS(otu_table_hel_dist,  k = 2, eig = T, add = T)
ord
## round(ord$eig*100/sum(ord$eig),1)[c(1,2)]
commun_NMDS1 <- ord$points[,1]
commun_NMDS2 <- ord$points[,2]

tp_diversity <- diversity[133:320,]
sem_dat <- data.frame(Site, commun_NMDS1, commun_NMDS2, tp_diversity,
                      climate_PC1, DOM_PC1, physio_chemical_PC1, ALT_PC1, Region1)

mod <- varpart(commun_NMDS1, ~ climate_PC1, ~ DOM_PC1, 
               ~ physio_chemical_PC1, ~ ALT_PC1,  data = sem_dat)
plot(mod)
mod <- varpart(otu_table_hel_dist, ~ climate_PC1, ~ DOM_PC1, 
               ~ physio_chemical_PC1, ~ ALT_PC1,  data = sem_dat)
plot(mod)

mod <- varpart(commun_NMDS1, ~ MAP + MAT, ~ SUVA254 + S275_295, 
               ~ pH + K, ~ ALT, data = env.table)
plot(mod)
part_fract <- mod$part$fract
indfract <- mod$part$indfract
contr1 <-mod$part$contr1
contr2 <-mod$part$contr1
## Add name to each set
install.packages("VennDiagram") 
library("VennDiagram")
grid.newpage()
draw.quad.venn(area1 = 0.6, 
               area2 = 0.36, 
               area3 = 0.32, 
               area4 = 0.01, 
               n12 = 0.33,
               n13 = 0.26, 
               n14 = 0.01, 
               n23 = 0.27, 
               n24 = 0.01,
               n34 = 0, 
               n123 = 0.24, 
               n124 = 0.01, 
               n134 = 0, 
               n234 = 0, 
               n1234 = 0,
               direct.area = T,
               fill = c("#56B4E9", "#009E73", "#0072B2", '#D55E00'),
               cat.col = c(rep('black',4)),cat.cex = 1.2,
               margin = 0.05,lty = 'blank',
               category = c("Climate elements", "DOM proporties", 
                            "physio chemical factors", "ALT"))



mod <- varpart(otu_table_hel_dist,  ~ MAP + MAT, ~ SUVA254 + S275_295, 
               ~ pH + K, ~ ALT, data = env.table)
plot(mod)
draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24,
               n34, n123, n124, n134, n234, n1234,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Group 1", "Group 2", "Group 3"))



## check collinearities among all variables， remove the varibles with vif > 10
ord <- capscale(otu_table_hel_dist ~., env.table)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(otu_table_hel_dist ~ MAP + MAT + DOC + 
                  S275_295 + SUVA254 + a300  + FluI + BIX + HIX +
                  DON + NH4_N + NO3_N + Depth + DO + 
                  pH + Conductivity + Ca + 
                  Mg + K + Na + ALT, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(otu_table_hel_dist ~ MAP + MAT + DOC + 
                  S275_295 + SUVA254 + a300  + FluI + BIX + HIX +
                  DON + NH4_N + NO3_N + Depth + DO + 
                  pH + Conductivity + Ca + 
                  Mg + K + ALT, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(otu_table_hel_dist ~ MAP + MAT + DOC + 
                  S275_295 + SUVA254 + a300  + FluI + BIX + HIX +
                  DON + NH4_N + NO3_N + Depth + DO + 
                  pH + Ca + 
                  Mg + K + ALT, env_df)
ord.vif <- vif.cca(ord)
ord.vif
names(ord.vif)
sel.env <- env.table[,names(ord.vif)]
## forward selection of environmental variables
mod1 <- capscale(otu_table_hel_dist ~. , sel.env, add=T)
mod0 <- capscale(otu_table_hel_dist ~1 , sel.env, add=T)
mod <- ordiR2step(mod0, scope = formula(mod1), perm.max = 999, trace = F)
Env.P <- anova.cca(mod,permutations = 999)[[4]][1]
Env.P
Env_sel <- mod$CCA$biplot
Env_sel <- sel.env[ , rownames(Env_sel)]

climate_vars <- Env_sel[ , c('MAP', 'MAT')]
DOM_vars <- Env_sel[ , c('SUVA254', 'S275_295', 'BIX', 'DOC', 'a300')]
phyche_vars <- Env_sel[ , c('pH', 'Depth', 'K', 'Ca', 'DON', 'NO3_N', 'DO', 'Mg')]
ALT <- Env_sel[ , c('ALT')]


mod <- varpart(commun_NMDS1, ~ MAP + MAT, ~ SUVA254 + S275_295 + BIX + DOC + a300,
               ~pH + Depth + K + Ca + DON + NO3_N + DO + Mg, ~ ALT, data = Env_sel)
plot(mod)

mod <- varpart(otu_table_hel_dist, ~ MAP + MAT, ~ SUVA254 + S275_295 + BIX + DOC + a300,
                 ~pH + Depth + K + Ca + DON + NO3_N + DO + Mg, ~ ALT, data = Env_sel)
plot(mod)


# SEM
## load the package
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(nlme)
##total community
env.table <- data.frame(sample_data(tp_phylo))
climate_vars <- env.table[ , c('MAP', 'MAT')]
DOM_vars <- env.table[ , c('SUVA254', 'S275_295')]
physio_chemical_vars <- env.table[ , c("pH", "K", "Na")]
ALT <- env.table[ , c('ALT')]
Region1 <- c(rep(1, 30), rep(2, 100), rep(3, 58))

## pca for climate
climate_PC1 <- princomp(climate_vars,cor=TRUE)
summary(climate_PC1, loadings=TRUE)
climate_PC1 <- predict(climate_PC1)
climate_PC1 <- climate_PC1[,1]

## pca for DOM
DOM_PC1 <- princomp(DOM_vars,cor=TRUE)
summary(DOM_PC1, loadings=TRUE)
DOM_PC1 <- predict(DOM_PC1)
DOM_PC1 <- DOM_PC1[,1]

## pca for physio_chemical factors
physio_chemical_PC1 <- princomp(physio_chemical_vars,cor=TRUE)
summary(physio_chemical_PC1, loadings=TRUE)
physio_chemical_PC1 <- predict(physio_chemical_PC1)
physio_chemical_PC1 <- physio_chemical_PC1[,1]

## pca for ALT
ALT_PC1 <- princomp(ALT,cor=TRUE)
summary(ALT_PC1, loadings=TRUE)
ALT_PC1 <- predict(ALT_PC1)
ALT_PC1 <- ALT_PC1[,1]

## site
Site <- env.table$Site

## ordination for community
otu_table <- as.matrix(t(otu_table(tp_phylo)))
otu_table_dist <- vegdist(otu_table, 'bray',upper=F)
ord <-  cmdscale(otu_table_dist,  k = 2, eig = T, add = T)
round(ord$eig*100/sum(ord$eig),1)[c(1,2)]
commun_NMDS1 <- ord$points[,1]
commun_NMDS2 <- ord$points[,2]

sem_dat <- data.frame(Site, commun_NMDS1, commun_NMDS2, 
                      climate_PC1, DOM_PC1, physio_chemical_PC1, ALT_PC1, Region1)
sem_total <- psem(
  lm(commun_NMDS1 ~ climate_PC1 + physio_chemical_PC1 + DOM_PC1, data = sem_dat),
  lm(ALT_PC1 ~ climate_PC1, data = sem_dat),
  lm(DOM_PC1 ~ climate_PC1 + ALT_PC1, data = sem_dat),
  lm(physio_chemical_PC1 ~ climate_PC1 + ALT_PC1 + DOM_PC1, data = sem_dat))
summary(sem_total)
sem_total <- psem(
  lm(commun_NMDS1 ~ climate_PC1 + physio_chemical_PC1 + DOM_PC1, data = sem_dat),
  lm(ALT_PC1 ~ climate_PC1, data = sem_dat),
  lm(DOM_PC1 ~ climate_PC1 + ALT_PC1, data = sem_dat),
  lm(physio_chemical_PC1 ~ ALT_PC1 + DOM_PC1, data = sem_dat))
summary(sem_total)
sem_total <- psem(
  lm(commun_NMDS1 ~ climate_PC1 + physio_chemical_PC1, data = sem_dat),
  lm(ALT_PC1 ~ climate_PC1, data = sem_dat),
  lm(DOM_PC1 ~ climate_PC1 + ALT_PC1, data = sem_dat),
  lm(physio_chemical_PC1 ~ ALT_PC1 + DOM_PC1, data = sem_dat))
summary(sem_total)
sem_total <- psem(
  lm(commun_NMDS1 ~ climate_PC1 + physio_chemical_PC1, data = sem_dat),
  lm(ALT_PC1 ~ climate_PC1, data = sem_dat),
  lm(DOM_PC1 ~ climate_PC1, data = sem_dat),
  lm(physio_chemical_PC1 ~ ALT_PC1 + DOM_PC1, data = sem_dat))
summary(sem_total)
fisherC(sem_total)
rsquared(sem_total)

coeffs_total <- coefs(sem_total, standardize = "scale")[,c(1,2,8)]
## Direct effect of climate variables on total community composition.
clim_dir_effect <- coeffs_total[1,3]
## Ind effect of climate variables on total community composition.
clim_ind_effect <- coeffs_total[4,3]*coeffs_total[6,3]*coeffs_total[2,3] +
  coeffs_total[3,3]*coeffs_total[6,3]*coeffs_total[2,3]
  
total_clim_effect <- clim_dir_effect + clim_ind_effect

## Direct effect of ALT on total community composition.
ALT_dir_effect <- 0
## Ind effect of ALT on total community composition.
ALT_ind_effect <- coeffs_total[5,3]*coeffs_total[2,3]
total_ALT_effect <- ALT_dir_effect + ALT_ind_effect

## Direct effect of DOM variables on total community composition.
DOM_dir_effect <- 0
## Ind effect of DOM variables on total community composition.
DOM_ind_effect <- coeffs_total[6,3]*coeffs_total[2,3]
total_DOM_effect <- DOM_dir_effect + DOM_ind_effect

## Direct effect of physio_chemical variables on total community composition.
physio_chemical_dir_effect <- coeffs_total[2,3]
## Ind effect of climate variables on total community composition.
physio_chemical_ind_effect <- 0
total_physio_chemical_effect <- physio_chemical_dir_effect + physio_chemical_ind_effect
## bar plot
sem_total_effect <- data.frame(variable = c(rep('Climate elements',3),
                                            rep('ALT', 3),
                                            rep('DOM properties', 3),
                                            rep('physio-chemical factors', 3)),
                                type_of_effects = c(rep(c('Total effects', 'Direct effects', 
                                                          'Indirect effects'), 4)),
                                value = c(total_clim_effect, clim_dir_effect, clim_ind_effect,
                                          total_ALT_effect, ALT_dir_effect, ALT_ind_effect,
                                          total_DOM_effect, DOM_dir_effect, DOM_ind_effect,
                                          total_physio_chemical_effect, physio_chemical_dir_effect,
                                          physio_chemical_ind_effect))
sem_total_effect$type_of_effects <- factor(sem_total_effect$type_of_effects, ordered = T,
                                            levels = c('Total effects', 'Direct effects', 'Indirect effects'))
sem_total_effect$variable <- factor(sem_total_effect$variable, ordered = T,
                                           levels = c('Climate elements', 'DOM properties', 'physio-chemical factors', 'ALT'))
sem_total_effect_plot <-  ggplot(data = sem_total_effect, 
                                  aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat="identity", width = 0.6) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-0.8, 0.5)) +
  labs(x = 'Variables', y = 'Standardized effects from SEM') +
  facet_wrap(~ type_of_effects) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(color='black',size=14),
        strip.text = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

## carbon community
climate_vars <- env.table[ , c('MAP', 'MAT')]
DOM_vars <- env.table[ , c('SUVA254', 'S275_295')]
physio_chemical_vars <- env.table[ , c("pH", "K", "Na")]

## pca for climate
climate_PC1 <- princomp(climate_vars,cor=TRUE)
summary(climate_PC1, loadings=TRUE)
climate_PC1 <- predict(climate_PC1)
climate_PC1 <- climate_PC1[,1]

## pca for DOM
DOM_PC1 <- princomp(DOM_vars,cor=TRUE)
summary(DOM_PC1, loadings=TRUE)
DOM_PC1 <- predict(DOM_PC1)
DOM_PC1 <- DOM_PC1[,1]

## pca for physio_chemical factors
physio_chemical_PC1 <- princomp(physio_chemical_vars,cor=TRUE)
summary(physio_chemical_PC1, loadings=TRUE)
physio_chemical_PC1 <- predict(physio_chemical_PC1)
physio_chemical_PC1 <- physio_chemical_PC1[,1]

## ordination for community
otu_table <- as.matrix(t(otu_table(tp_carbon_phy)))
otu_table_dist <- vegdist(otu_table, 'bray',upper=F)
ord <-  cmdscale(otu_table_dist,  k = 2, eig = T, add = T)
round(ord$eig*100/sum(ord$eig),1)[c(1,2)]
commun_NMDS1 <- ord$points[,1]
commun_NMDS2 <- ord$points[,2]


sem_carbon_dat <- data.frame(Site, commun_NMDS1, commun_NMDS2, climate_PC1, 
                             DOM_PC1, physio_chemical_PC1,  ALT_PC1)
sem_carbon <- psem(
  lm(commun_NMDS1 ~ climate_PC1 + physio_chemical_PC1 + DOM_PC1, data = sem_carbon_dat),
  lm(ALT_PC1 ~ climate_PC1, data = sem_carbon_dat),
  lm(DOM_PC1 ~ climate_PC1 + ALT_PC1, data = sem_carbon_dat),
  lm(physio_chemical_PC1 ~ climate_PC1 + ALT_PC1 + DOM_PC1, data = sem_carbon_dat))
summary(sem_carbon)
sem_carbon <- psem(
  lm(commun_NMDS1 ~ climate_PC1 + physio_chemical_PC1 + DOM_PC1, data = sem_carbon_dat),
  lm(ALT_PC1 ~ climate_PC1, data = sem_carbon_dat),
  lm(DOM_PC1 ~ climate_PC1 + ALT_PC1, data = sem_carbon_dat),
  lm(physio_chemical_PC1 ~ ALT_PC1 + DOM_PC1, data = sem_carbon_dat))
summary(sem_carbon)
sem_carbon <- psem(
  lm(commun_NMDS1 ~ climate_PC1 + physio_chemical_PC1 + DOM_PC1, data = sem_carbon_dat),
  lm(ALT_PC1 ~ climate_PC1, data = sem_carbon_dat),
  lm(DOM_PC1 ~ climate_PC1, data = sem_carbon_dat),
  lm(physio_chemical_PC1 ~ ALT_PC1 + DOM_PC1, data = sem_carbon_dat))
summary(sem_carbon)
fisherC(sem_carbon)
rsquared(sem_carbon)
## Direct effects (all)
coeffs_carbon <- coefs(sem_carbon, standardize = "scale")[,c(1,2,8)]

## Direct effect of climate variables on total community composition.
clim_dir_effect <- coeffs_carbon[1,3]
## Ind effect of climate variables on total community composition.
clim_ind_effect <- coeffs_carbon[5,3]*coeffs_carbon[3,3]+
  coeffs_carbon[5,3]*coeffs_carbon[6,3]*coeffs_carbon[2,3]+
  coeffs_carbon[4,3]*coeffs_carbon[6,3]*coeffs_carbon[2,3]
total_clim_effect <- clim_dir_effect + clim_ind_effect

## Direct effect of ALT on total community composition.
ALT_dir_effect <- 0
## Ind effect of ALT on total community composition.
ALT_ind_effect <- coeffs_total[6,3]*coeffs_total[2,3]
total_ALT_effect <- ALT_dir_effect + ALT_ind_effect

## Direct effect of DOM variables on total community composition.
DOM_dir_effect <- coeffs_carbon[3,3]
## Ind effect of DOM variables on total community composition.
DOM_ind_effect <- coeffs_carbon[7,3]*coeffs_carbon[2,3]
total_DOM_effect <- DOM_dir_effect + DOM_ind_effect

## Direct effect of physio_chemical variables on total community composition.
physio_chemical_dir_effect <- coeffs_carbon[2,3]
# Ind effect of climate variables on total community composition.
physio_chemical_ind_effect <- 0
total_physio_chemical_effect <- physio_chemical_dir_effect + physio_chemical_ind_effect

## bar plot
sem_carbon_effect <- data.frame(variable = c(rep('Climate elements',3),
                                            rep('ALT', 3),
                                            rep('DOM properties', 3),
                                            rep('physio-chemical factors', 3)),
                               type_of_effects = c(rep(c('Total effects', 'Direct effects', 
                                                         'Indirect effects'), 4)),
                               value = c(total_clim_effect, clim_dir_effect, clim_ind_effect,
                                         total_ALT_effect, ALT_dir_effect, ALT_ind_effect,
                                         total_DOM_effect, DOM_dir_effect, DOM_ind_effect,
                                         total_physio_chemical_effect, physio_chemical_dir_effect,
                                         physio_chemical_ind_effect))
sem_carbon_effect$type_of_effects <- factor(sem_carbon_effect$type_of_effects, ordered = T,
                                           levels = c('Total effects', 'Direct effects', 'Indirect effects'))
sem_carbon_effect$variable <- factor(sem_carbon_effect$variable, ordered = T,
                                    levels = c('Climate elements', 'DOM properties', 'physio-chemical factors', 'ALT'))
sem_carbon_effect_plot <-  ggplot(data = sem_carbon_effect, 
                                 aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat="identity", width = 0.6) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-0.8, 0.5)) +
  labs(x = 'Variables', y = 'Standardized effects from SEM') +
  facet_wrap(~ type_of_effects) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(color='black',size=14),
        strip.text = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

library(cowplot)
p_sem_effects <- plot_grid(sem_total_effect_plot, sem_carbon_effect_plot, 
                           labels = 'auto', ncol = 1, nrow = 2, 
                           label_x = .01, label_y = 1.05, hjust = 0, 
                           label_size=14, align = "v")
p_sem_effects


res_fun_sample_df <- cbind(m2$res_spe_func_perc[c(133:320), ], m2$sample_table[c(133:320), ])
res_carbon_fun_sample_df <- res_fun_sample_df[,c("cellulolysis", "chitinolysis", "xylanolysis", 
                                                 'aromatic_hydrocarbon_degradation','aromatic_compound_degradation',
                                                 "fermentation", 'hydrocarbon_degradation',
                                                 "methanotrophy", "methanogenesis", 
                                                 "methanol_oxidation", "methylotrophy", "Site",
                                                 "DOC", "DON", "NH4_N", "NO3_N", "S275_295", "SUVA254",
                                                 "a300", "MAP", "MAT", "Depth", "DO", "pH", "Conductivity", 
                                                 "Salinity", "Ca", "Mg", "K", "Na","FluI", 
                                                 "FrI", "BIX", "HIX", "ALT")]

y <- res_carbon_fun_sample_df[,c("cellulolysis", "chitinolysis", "xylanolysis",
                          'aromatic_hydrocarbon_degradation','aromatic_compound_degradation',
                          "fermentation", 'hydrocarbon_degradation',
                          "methanotrophy", "methanogenesis",
                          "methanol_oxidation", "methylotrophy")]
x <- res_carbon_fun_sample_df[,c("DOC", "DON", "NH4_N", "NO3_N", "S275_295", "SUVA254",
                                 "a300", "MAP", "MAT", "Depth", "DO", "pH", "Conductivity", 
                                 "Salinity", "Ca", "Mg", "K", "Na","FluI", 
                                 "FrI", "BIX", "HIX", "ALT")]

## linner mixed model for matrix to matrix
lmm.mat.cal <- function(y, x, total_df){
  require(lme4)
  require(lmerTest)
  require(MuMIn)
  y <- as.matrix(y)
  x <- as.matrix(x)
  df<-NULL
  for(i in colnames(y)){
    for(j in colnames(x)){
      a <- y[, i, drop = F]
      b <- x[, j, drop = F]
      mode <- lmer(a ~ b + (1|Site), data = total_df, na.action=na.omit)
      coeff <- summary(mode)$coefficients[2,1]
      r.square <- r.squaredGLMM(mode)[1]
      if (coeff>0) r = sqrt(r.square)
      else r = (-1) * sqrt(r.square)
      tmp <- c(i, j, r, r.square, anova(mode)$Pr)
      if(is.null(df)){
        df <- tmp  
      }
      else{
        df <- rbind(df, tmp)
      }    
    }
  }
  df<-data.frame(row.names=NULL,df)
  colnames(df)<-c("carbon_funs","Env","Correlation","r.square", "Pvalue")
  df$Pvalue<-as.numeric(as.character(df$Pvalue))
  df$AdjPvalue<-rep(0,dim(df)[1])
  df$Correlation<-as.numeric(as.character(df$Correlation))
  #You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
  # 1 -> donot adjust
  # 2 -> adjust Env + Type (column on the correlation plot)
  # 3 -> adjust carbon_funs + Type (row on the correlation plot for each type)
  # 4 -> adjust carbon_funs (row on the correlation plot)
  # 5 -> adjust Env (panel on the correlation plot)
  adjustment_label<-c("NoAdj","AdjEnvAndType","Adjcarbon_funsAndType","Adjcarbon_funs","AdjEnv")
  adjustment<-5
  if(adjustment==1){
    df$AdjPvalue<-df$Pvalue
  } else if (adjustment==2){
    for(i in unique(df$Env)){
      for(j in unique(df$Type)){
        sel<-df$Env==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
  } else if (adjustment==3){
    for(i in unique(df$carbon_funs)){
      for(j in unique(df$Type)){
        sel<-df$carbon_funs==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
  } else if (adjustment==4){
    for(i in unique(df$carbon_funs)){
      sel<-df$carbon_funs==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  } else if (adjustment==5){
    for(i in unique(df$Env)){
      sel<-df$Env==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
  #Now we generate the labels for signifant values
  df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  df$carbon_funs <-factor(df$carbon_funs, ordered = T, levels = rev(colnames(y)))
  df$Env <-factor(df$Env, ordered = T, levels = colnames(x))
  return(df)
}

lmm.matrix <- lmm.mat.cal(y, x, res_carbon_fun_sample_df)
lmm.matrix$Correlation[lmm.matrix$AdjPvalue >= 0.05] <- 0
## plot
lmm_carbon_fun_env <- ggplot(aes(x = Env, y = carbon_funs, fill = Correlation), data = lmm.matrix)+
  geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C")+
  geom_text(aes(label=Significance), color="black", size = 5.5)+
  labs(y = 'Carbon cycling', x = 'Environmental factors', fill= 'Correlation') +
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
lmm_carbon_fun_env























multi.plot <- function(a, b){
  tmp_plot <- ggplot(avg_diversity, aes(x = a, y = b))+
    geom_point(shape=19, size=2,colour='tomato3',alpha=0.8)+
    geom_smooth(method="lm", size=1.5, se=T,colour='black') +
    #scale_y_continuous(expand = c(0,0))+
    #scale_x_continuous(expand = c(0,0))+
    theme(axis.line=element_line(colour = 'black'))+
    theme_classic()
  return(tmp_plot)
}




res_carbon_fun_avg_sample_df <- res_carbon_fun_sample_df %>% group_by(Site) %>% summarise_all(mean)

#test correlation between env and diversity
library(lme4)
library(lmerTest)
library(MuMIn)
#Identify, describe, plot, and remove the outliers from the dataset
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
  if(response == "y" | response == "yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "n")
    return(invisible(dt))
  } else{
    cat("Nothing changed", "n")
    return(invisible(var_name))
  }
}
colnames(res_carbon_fun_avg_sample_df)[-1]
outlierKD(res_carbon_fun_avg_sample_df, ALT)
yes






M1 <- lm(chitinolysis ~ MAP, res_carbon_fun_sample_df)
summary(M1)
performance::r2(M1)
parameters::p_value(M1)

# the reationship between biomarker with environmental factors
biomark_c_species <- c('Acinetobacter', 'Rhodoferax', 'Candidatus Methylopumilus', 'Anoxybacillus',
                       'Halomonas', 'Methanobacterium', 'Methylobacter', 'Paracoccus', 'Methylotenera',
                       'Desulfotignum', 'Nocardioides', 'Methanoregula', 'Geothrix', 'Lysobacter',
                       'Aeromonas', 'Dyadobacter', 'Methylophaga', 'Methylobacterium', 'Methylocapsa',
                       'Methylomonas')
biomarker_phylo <- subset_taxa(tp_phylo, Genus %in% biomark_c_species)

#edgeR analysis
# Glom OTUs to genus level for further statistical analysis & reasonable power
ps1 <- tax_glom(tp_carbon_phy, "Genus", NArm = TRUE)
#Use code snippet straight provided by authors of PhyloSeq to export phyloseq object to an EdgeR object
phyloseq_to_edgeR = function(physeq, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

# Make normalized phyloseq object (ps6) into an edgeR object. It needs a grouping factor. We use location.
dge = phyloseq_to_edgeR(ps1)
#The crunching to follow is much easier if metadata is pulled out this way into object "a"
a = sample_data(ps1)

#So many things don't understand what kind of variables they need to be. Make sure they understand
MAT <- as.numeric(a$MAT)
MAP <- as.numeric(a$MAP)
DOC <- as.numeric(a$DOC)
S275_295 <- as.numeric(a$S275_295)
SUVA254 <- as.numeric(a$SUVA254) 
a300 <- as.numeric(a$a300)
Flul <- as.numeric(a$Flul)
Frl <- as.numeric(a$Frl)
Bix <- as.numeric(a$Bix)
Hix <- as.numeric(a$Hix)
DON <- as.numeric(a$DON)
NH4_N <- as.numeric(a$NH4_N) 
NO3_N <- as.numeric(a$NO3_N)
Depth <- as.numeric(a$Depth)
DO <- as.numeric(a$DO)
pH <- as.numeric(a$pH)
conductivity <- as.numeric(a$Conductivity)
Salinity <- as.numeric(a$Salinity) 
K <- as.numeric(a$K)
Ca <- as.numeric(a$Ca)
Na <- as.numeric(a$Na)
Mg <- as.numeric(a$Mg)
ALT <- as.numeric(a$ALT)

# Design for my linear model
design <-model.matrix(~ MAP + MAT + S275_295 + SUVA254 + pH + Na + K + ALT)

# EdgeR needs to calculate dispersion again after you've fed it the design.
#This step projectile vomits if chloride or conductivity aren't sqrt transformed.
x = calcNormFactors(dge, method="RLE")
x = estimateGLMCommonDisp(dge, design)
x = estimateGLMTrendedDisp(dge, design)
x = estimateGLMTagwiseDisp(dge, design)

fit <-glmFit(x, design)

# grab the coefficients I care about
lrt <- glmLRT(fit, coef=2:9)

# lrt to z-scores
table <- lrt$table
table <- apply(table, 2, function(x) scale(x, center = TRUE, scale = TRUE))
rownames(table) = rownames(lrt$table)
table <- table[,1:8]
table <- as.data.frame(table)

q <- lrt$genes
table1 <- cbind(table, Kingdom = q$Kingdom, Phylum = q$Phylum, Class = q$Class, Order = q$Order,
              Family = q$Family, Genus = q$Genus, OTU = rownames(table))

# Sometimes want to do this to check what has at least one z score beyond threshold
table2 <- table1[apply(table1[,1:8], 1, function(x) any(abs(x)>1.96)), ]

#write.csv(table2, "./tables/EdgeR-Zscores-Genera.csv")
melted <- melt(table2)
melt2 <- subset(melted, abs(melted$value) > 1.96 )

filter <- unique(melt2$OTU)
mini <- prune_taxa((rownames(otu_table(tp_phylo)) %in% filter), tp_phylo)
f <- phy_tree(mini)
tax <- as.data.frame(tax_table(mini))

tax$Genus <- gsub("D_5__", "", tax$Genus)
tax$Phylum <- gsub("D_1__", "", tax$Phylum)

# change tip labels to genera / need to gsub
library(ggtree)
mytree <- phy_tree(f)$tip.label
foo <- data.frame(label = mytree, label2 = paste(tax$Genus, tax$OTU, sep='_'), 
                  label3 = tax$Phylum, stringsAsFactors = F)
#foo$label2[foo$label2 == NA] <- 'uncultured'
# Make a less putrid color palette

library(RColorBrewer)

colourCount <- length(unique(tax$Phylum))
Mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(colourCount)
names(Mycolors) <- levels(as.factor(tax$Phylum))
colScale <- scale_colour_manual(name = "grp",values = Mycolors)

p<-ggtree(f) %<+% foo + 
  geom_tiplab(size = 3, align = TRUE, linesize = 0.5, aes(label = label2), hjust = 0) +
  geom_tippoint(size = 2, aes(label = label2, group = label3, color = label3) ) + 
  theme_tree2() + 
  scale_color_manual('Phylum', values = colorRampPalette(brewer.pal(11, "Spectral"))(colourCount))
p + xlim (0, 3) 

test<-as.matrix(table2[1:8])
tested<-apply(test, 2, function(x){cut(x, br=c(-8, -4, -2, 2, 4, 8))})
rownames(tested) = rownames(table2)
colnames(tested) = gsub("logFC.", "", colnames(tested))
colnames(tested) = gsub("logCPM", "Rotorua", colnames(tested))

heatmapcols <-colorRampPalette(brewer.pal(6, "Spectral"))(6)
names(heatmapcols) <- levels(as.factor(tested[1:10000]))
# Make sure to check that the order of these is correct if this is rerunhe
heatmapcols[1] = "#F2F2F2"
heatmapcols[2] = "#AABBDD"
heatmapcols[3] = "#112288"
heatmapcols[4] = "#FFCCCF"
heatmapcols[5] = "#881111"
heatmapcols[6] = "#112288" 

p1 <- gheatmap(p, tested, offset = 0.6, width = 2.5, font.size = 2.5, 
             colnames_angle = -90, colnames_position = "top", hjust = 1) + 
  scale_fill_manual('Z-scores', values = heatmapcols)
p1 <- p1 + theme(legend.position = "right")
p1
#ggsave(filename="output/edgeR-Genus.pdf", plot=p, width=8, height=10, units="in")

ref.seqs <- refseq(methane_phylo)
tree <- phy_tree(methane_phylo)
tax <- tax_table(methane_phylo)
otu.table <- otu_table(methane_phylo)

tax_glom(biomarker_phylo, "Genus", NArm = TRUE)
# the relationship between SUVA254 and MAT, s275_295
mode1 <- lm(SUVA254~MAT, data = env.table)
summary(mode1)

mode2 <- lm(S275_295~MAT, data = env.table)
summary(mode2)

mode3 <- lm(SUVA254~MAP, data = env.table)
summary(mode3)

mode4 <- lm(S275_295~MAP, data = env.table)
summary(mode4)

p1 <- ggplot(env.table, aes(x = MAT, y = SUVA254))+
  geom_point(shape = 19, size = 2, colour = "#DD5F60")+
  geom_smooth(method = "lm", size = 1.5, se = T,colour = 'black') +
  labs(x= 'MAT (℃)', y = 'SUVA254') +
  theme(axis.line = element_line(colour = 'black'))+
  theme_classic()+
  theme(axis.title = element_text(color='black',size=14),
        axis.text = element_text(colour='black',size=12))

p2 <- ggplot(env.table, aes(x = MAT, y = S275_295))+
  geom_point(shape = 19, size = 2, colour = "#DD5F60")+
  geom_smooth(method = "lm", size = 1.5, se = T,colour = 'black') +
  labs(x= 'MAT (℃)', y = 'S275_295') +
  theme(axis.line = element_line(colour = 'black'))+
  theme_classic()+
  theme(axis.title = element_text(color='black',size=14),
        axis.text = element_text(colour='black',size=12))

p3 <- ggplot(env.table, aes(x = MAP, y = SUVA254))+
  geom_point(shape = 19, size = 2, colour = "#DD5F60")+
  geom_smooth(method = "lm", size = 1.5, se = T,colour = 'black') +
  labs(x= 'MAP (mm)', y = 'SUVA254') +
  theme(axis.line = element_line(colour = 'black'))+
  theme_classic()+
  theme(axis.title = element_text(color='black',size=14),
        axis.text = element_text(colour='black',size=12))

p4 <- ggplot(env.table, aes(x = MAP, y = S275_295))+
  geom_point(shape = 19, size = 2, colour = "#DD5F60")+
  geom_smooth(method = "lm", size = 1.5, se = T,colour = 'black') +
  labs(x= 'MAP (mm)', y = 'S275_295') +
  theme(axis.line = element_line(colour = 'black'))+
  theme_classic()+
  theme(axis.title = element_text(color='black',size=14),
        axis.text = element_text(colour='black',size=12))

DOM_climate_plot <- plot_grid(p1, p2, p3, p4,
                              labels = c('(a)', '(b)', '(c)', '(d)'), 
                              ncol = 2, nrow = 2, 
                              label_x = .01, label_y = 1.01, hjust = 0, 
                              label_size=14, align = "v")
DOM_climate_plot






























##nitrogen community
climate_vars <- env.table[ , c('MAP')]
DOM_vars <- env.table[ , c('SUVA254', 'S275_295', "comp2", "comp4")]
physic_vars <- env.table[ , c("pH", "K", "Na")]

library(vegan)
climate_dist <- vegdist(climate_vars, 'euclidean',upper=F)
decorana(climate_dist)
sol <- rda(climate_vars,trace = F)
sol
climate_PC1 <- scores(sol)$sites[,1]

DOM_dist <- vegdist(DOM_vars, 'euclidean',upper=F)
decorana(DOM_dist)
sol <- rda(DOM_vars,trace = F)
sol
DOM_PC1 <- scores(sol)$sites[,1]

physic_dist <- vegdist(physic_vars, 'euclidean',upper=F)
decorana(physic_dist)
sol <- rda(physic_vars,trace = F)
sol
physic_PC1 <- scores(sol)$sites[,1]

otu_table <- as.matrix(t(otu_table(tp_carbon_phy)))
otu_table_hel <- decostand(otu_table, 'hellinger')
MDS <- cmdscale(vegdist(otu_table_hel, method = "bray"), k = 2, eig = T, add = T )
ord_NMDS <- metaMDS(otu_table_hel, k = 2, trymax = 100, trace = F)
round(MDS$eig*100/sum(MDS$eig),1)
commun_NMDS1 <- ord_NMDS$points[,1]

sem_carbon_dat <- data.frame(commun_NMDS1, climate_PC1, DOM_PC1, physic_PC1)
library(piecewiseSEM)
sem_carbon <- psem(
  lm(commun_NMDS1 ~ climate_PC1 + DOM_PC1 + physic_PC1, data = sem_carbon_dat),
  lm(DOM_PC1 ~ climate_PC1, data = sem_carbon_dat),
  lm(physic_PC1 ~ climate_PC1, data = sem_carbon_dat),
  physic_PC1 %~~% DOM_PC1
)
summary(sem_carbon)



# create object of trans_func
meco_tp <- phyloseq2meco(tp_phylo)
t2 <- trans_func$new(meco_tp)

# mapping the taxonomy to the database
# the function can recognize prokaryotes or fungi automatically.
t2$cal_spe_func()
# return t2$res_spe_func, 1 represent function exists, 0 represent no or cannot confirmed.
t2$res_spe_func[1:5, 1:2]
# If you want to change the group list, reset the list t2$func_group_list
C_cycle_process <- c("methanotrophy",  "methanogenesis", 
                     "methanol_oxidation", "methylotrophy", "chitinolysis",
                     "cellulolysis", "xylanolysis", "ligninolysis", "fermentation",
                     'aromatic_hydrocarbon_degradation','aromatic_compound_degradation',
                     'hydrocarbon_degradation')
t2$func_group_list$`C-cycle` <- C_cycle_process
t2$func_group_list
# use show_prok_func to see the detailed information of prokaryotic traits
t2$show_prok_func("methanotrophy")
# calculate the percentages for communities
t2$cal_spe_func_perc(use_community = TRUE)
# t2$res_spe_func_perc[1:5, 1:2]
# then we try to correlate the res_spe_func_perc of communities to environmental variables

t3 <- trans_env$new(dataset = meco_dat, add_data = t2$sample_table[, colnames(Env_C_sel)])
t3$cal_cor(add_abund_table = t2$res_spe_func_perc[ ,colnames(t2$res_spe_func_perc) %in% 
                                                     c(t2$func_group_list$`C-cycle`)], cor_method = "spearman")
t3$plot_corr(pheatmap = TRUE)




#calculate the relative abundance of C cycle process
c_cycle_perc <- t2$res_spe_func_perc[ ,colnames(t2$res_spe_func_perc) %in% c(t2$func_group_list$`C-cycle`)]
group <- c(rep(1, 20), rep(2, 124), rep(3, 44))
group <- as.factor(group)
meta_data <- data.frame(group, c_cycle_perc)
melt_df <- melt(meta_data, id.vars = c('group'))
colnames(melt_df) <- c('group', 'cat_carbon','relative_abundance')
mode <- lmer(relative_abundance ~ cat_carbon + (1|group), data = melt_df)
anova(mode)
library(multcomp)
summary(glht(mode, linfct = mcp(cat_carbon = "Tukey")), test = adjusted("holm"))

new_df <- ddply(melt_df, c("cat_carbon"), summarise,
                mean = mean(relative_abundance), sd = sd(relative_abundance),
                se = sd(relative_abundance)/sqrt(length(relative_abundance)))
new_df <- data.frame(new_df, sig = c('bc', 'a', 'd', 'ab', 'c', 'e'))
ggplot(new_df,aes(x = cat_carbon, y = mean, fill = cat_carbon))+
  geom_bar(stat = 'identity', width = 0.7)+
  #scale_fill_manual(values=c("#F2B379","#479E9B","#FDDC7B","#4169B2","#ACCDDC","#DD5F60","#B1A4C0"))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se), width=.2)+
  scale_x_discrete(limits = c("chitinolysis", "cellulolysis", "fermentation", 
                              "methanogenesis", "methanotrophy", "methylotrophy"))+
  labs(x = 'Carbon cycle', y = 'Mean relative abundance (%)',
       fill = "Carbon cycle")+
  scale_y_continuous(expand = c(0,0), limits = c(0,5))+
  #geom_text(aes(label = sig, y = mean + se+0.02*max(mean)), position = position_dodge(0.9), vjust = 0)+
  theme_bw()+
  theme(legend.position = c(0.85,0.8),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

#extract asvs involved in carbon cycling
c_cycle <- t2$res_spe_func[ ,c(t2$func_group_list$`C-cycle`)]
c_cycle_asvs <- rownames(c_cycle[rowSums(c_cycle) != 0,])

carbon_phy <- subset_taxa(water_physeq, ASV %in% c_cycle_asvs)
carbon_phy_rel <- subset_taxa(water.rel, ASV %in% c_cycle_asvs)
env.table <- read.csv(file = './edge_analysis/env.csv', header = T, row.names = 1, stringsAsFactors = F)
sample_data(carbon_phy_rel) <- env.table

#phylogenetic tree
library(ggtreeExtra)
library(dplyr)
ps6 <- tax_glom(tp_carbon_phy, taxrank="Genus")
taxa_sums(ps6)
melt_simple <- data.frame(Phylum = data.frame(tax_table(ps6))[2],
                          Genus = data.frame(tax_table(ps6))[6],
                          OTU = taxa_names(ps6),
                          abundance = taxa_sums(ps6))
f<-phy_tree(ps6)
tax = as.data.frame(tax_table(ps6))
# change tip labels to genera / need to gsub
library(ggtree)
mytree = phy_tree(f)$tip.label
foo = data.frame(label=mytree, label2=tax$Genus, label3=tax$Phylum, stringsAsFactors = F)
#foo$label2[foo$label2 == NA] <- 'uncultured'
# Make a less putrid color palette
library(RColorBrewer)
colourCount = length(unique(tax$Phylum))
Mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(colourCount)
names(Mycolors) <- levels(as.factor(tax$Phylum))
colScale <- scale_colour_manual(name = "grp",values = Mycolors)

p<-ggtree(f, branch.length='none', layout='circular', open.angle = 20) %<+% foo + 
  geom_tiplab(size=2, align=TRUE, linesize=0.5, aes(label=label2), offset=5.5, hjust =0) +
  geom_tippoint(size=2, aes(label=label2, group=label3, color=label3) ) + 
  theme_tree2() + 
  scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(colourCount))

p_tree <- p + geom_fruit(
  data = melt_simple,
  geom = geom_bar,
  mapping = aes(x = abundance, y = OTU, fill = 'red4'),
  orientation="y",
  stat="identity",
  size=.2,
  outlier.size=0.5,
  outlier.stroke=0.08,
  outlier.shape=21,
  axis.params=list(
    axis       = "x",
    text.size  = 1.8,
    hjust      = 1,
    vjust      = 0.5,
    nbreak     = 3,
  ),
  grid.params=list()
) 
p_tree
#comunity composition
##determine the genus compositions within each family##
ps6 <- tax_glom(carbon_phy_rel, taxrank="Genus")

genus.ra.table <- otu_table(ps6)
MRA <- rowMeans(genus.ra.table)
group <- tax_table(ps6)[,c(5,6)]
genus.mra.table <- data.frame(group,MRA)

#arrange the genuss table
library(tidyr)
genus.mra.table <- genus.mra.table %>% spread(Genus, MRA)
genus.mra.table[is.na(genus.mra.table)] <- 0
rownames(genus.mra.table)<-genus.mra.table$Family
genus.mra.table<-as.matrix(t(genus.mra.table[,-1])*100)
colsum <-apply(genus.mra.table,2,sum)
rowsum<-apply(genus.mra.table,1,sum)
topgenus_table<-(genus.mra.table[order(rowsum,decreasing=TRUE),order(colsum,decreasing=TRUE)])[,1:10]
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

#edgeR analysis
# Glom OTUs to genus level for further statistical analysis & reasonable power
ps6 <- tax_glom(carbon_phy_rel, "Genus", NArm = TRUE)
#Use code snippet straight provided by authors of PhyloSeq to export phyloseq object to an EdgeR object
phyloseq_to_edgeR = function(physeq, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

# Make normalized phyloseq object (ps6) into an edgeR object. It needs a grouping factor. We use location.
dge = phyloseq_to_edgeR(ps6)
#The crunching to follow is much easier if metadata is pulled out this way into object "a"
a = sample_data(ps6)

#So many things don't understand what kind of variables they need to be. Make sure they understand
chitinolysis <- as.numeric(a$chitinolysis)
cellulolysis <- as.numeric(a$cellulolysis)
fermentation <- as.numeric(a$fermentation)
methanogenesis <- as.numeric(a$methanogenesis)
methanotrophy <- as.numeric(a$methanotrophy)
methylotrophy <- as.numeric(a$methylotrophy)
DOC <- as.numeric(a$DOC)
TN <- as.numeric(a$TN)
NH4_N <- as.numeric(a$NH4_N) 
NO3_N <- as.numeric(a$NO3_N)
DON <- as.numeric(a$DON)
S275_295 <- as.numeric(a$S275_295)
SUVA254 <- as.numeric(a$SUVA254) 
conductivity <- as.numeric(a$Conductivity)
a300 <- as.numeric(a$a300)
fi <- as.numeric(a$fi)
bix <- as.numeric(a$bix)
hix <- as.numeric(a$hix)
Comp1 <- as.numeric(a$C1...)
Comp2 <- as.numeric(a$C2...)
Comp3 <- as.numeric(a$C3...)
Comp4 <- as.numeric(a$C4...)
pH <- as.numeric(a$pH)
DO <- as.numeric(a$DO)
Salinity <- as.numeric(a$Salinity) 
Depth <- as.numeric(a$Depth)
MAT <- as.numeric(a$MAT)
MAP <- as.numeric(a$MAP)
Conductivity <- as.numeric(a$Conductivity) # transforming this fixes it so it doesn't kill linear model
K <- as.numeric(a$K)
Ca <- as.numeric(a$Ca)
Na <- as.numeric(a$Na)
Mg <- as.numeric(a$Mg)
##Test and forward selection of environmental variables
Env_select <- function (Dist_Matrix,Env,Number_Permutations=999) {
  mod1<-capscale(Dist_Matrix~.,Env,add=T,na.action = na.omit)
  mod0<-capscale(Dist_Matrix~1,Env,add=T,na.action = na.omit)
  mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999,trace = F)
  Env.P<-anova.cca(mod,permutations=Number_Permutations)[[4]][1]
  Env_se<-mod$CCA$biplot
  Env_se<-Env[,rownames(Env_se)]
  return(Env_se)
}
carbon_asvs <- otu_table(carbon_phy_rel)
carbon.dist <- vegdist(sqrt(t(carbon_asvs)))
env.table <- sample_data(carbon_phy_rel)
env_vars <- data.frame(scale(env.table[ ,-c(1:6)]))
env_se_carbon <- Env_select (carbon.dist, env_vars)
colnames(env_se_carbon)
# Design for my linear model
design <-model.matrix(~ chitinolysis + cellulolysis + fermentation + methanogenesis + 
                        methanotrophy + methylotrophy + pH + MAP + MAT +
                        Comp4 + Comp3 + S275_295 + SUVA254 + DOC + NH4_N +
                        Na + K + Ca + Conductivity + Depth) 

# EdgeR needs to calculate dispersion again after you've fed it the design.
#This step projectile vomits if chloride or conductivity aren't sqrt transformed.
x = calcNormFactors(dge, method="RLE")
x = estimateGLMCommonDisp(dge, design)
x = estimateGLMTrendedDisp(dge, design)
x = estimateGLMTagwiseDisp(dge, design)

fit <-glmFit(x, design)

# grab the coefficients I care about
lrt <- glmLRT(fit, coef=2:21)

# lrt to z-scores
table<-lrt$table
table<-apply(table, 2, function(x) scale(x, center = TRUE, scale = TRUE))
rownames(table) = rownames(lrt$table)
table<-table[,1:20]
table<-as.data.frame(table)

q<-lrt$genes
table1<-cbind(table, Kingdom = q$Kingdom, Phylum = q$Phylum, Class = q$Class, Order = q$Order,
              Family = q$Family, Genus = q$Genus, OTU = rownames(table))

# Sometimes want to do this to check what has at least one z score beyond threshold
table2 <- table1[apply(table1[,1:20], 1, function(x) any(abs(x)>1.96)), ]

#write.csv(table2, "./tables/EdgeR-Zscores-Genera.csv")
melted<-melt(table2)
melt2<-subset(melted, abs(melted$value) > 1.96 )

filter<-unique(melt2$OTU)
mini<-prune_taxa((rownames(otu_table(water_physeq)) %in% filter), water_physeq)
f<-phy_tree(mini)
tax = as.data.frame(tax_table(mini))

tax$Genus = gsub("D_5__", "", tax$Genus)
tax$Phylum = gsub("D_1__", "", tax$Phylum)

# change tip labels to genera / need to gsub
library(ggtree)
mytree = phy_tree(f)$tip.label
foo = data.frame(label=mytree, label2=paste(tax$Genus, tax$ASV, sep='_'), label3=tax$Phylum, stringsAsFactors = F)
#foo$label2[foo$label2 == NA] <- 'uncultured'
# Make a less putrid color palette

library(RColorBrewer)

colourCount <- length(unique(tax$Phylum))
Mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(colourCount)
names(Mycolors) <- levels(as.factor(tax$Phylum))
colScale <- scale_colour_manual(name = "grp",values = Mycolors)

p<-ggtree(f) %<+% foo + 
  geom_tiplab(size=2, align=TRUE, linesize=0.5, aes(label=label2), hjust =0) +
  geom_tippoint(size=2, aes(label=label2, group=label3, color=label3) ) + 
  theme_tree2() + 
  scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(colourCount))
p + xlim (0, 3) 

test<-as.matrix(table2[1:20])
tested<-apply(test, 2, function(x){cut(x, br=c(-14, -4, -2, 2, 4, 14))})
rownames(tested) = rownames(table2)
colnames(tested) = gsub("logFC.", "", colnames(tested))
colnames(tested) = gsub("logCPM", "Rotorua", colnames(tested))

heatmapcols <-colorRampPalette(brewer.pal(6, "Spectral"))(6)
names(heatmapcols) <- levels(as.factor(tested[1:10000]))
heatmapcols <- c(heatmapcols[5],heatmapcols[4],heatmapcols[2],heatmapcols[3],heatmapcols[1],heatmapcols[6])
# Make sure to check that the order of these is correct if this is rerunhe
heatmapcols[1] = "#881111"
heatmapcols[2] = "#FFCCCF"
heatmapcols[3] = "#F2F2F2"
heatmapcols[4] = "#AABBDD"
heatmapcols[5] = "#112288"
heatmapcols[6] = "#112288"

mycols <- c("#881111","#D9444D","#FFCCCF","#F2F2F2","#AABBDD","#112288")
tested <- tested[,c('chitinolysis', 'cellulolysis', 'fermentation', 'methanogenesis', 
                    'methanotrophy', 'methylotrophy',  'pH', 'MAP', 'MAT',
                    'Comp4', 'Comp3', 'S275_295', 'SUVA254', 'DOC', 'NH4_N',
                    'Na', 'K', 'Ca', 'Conductivity', 'Depth')] 
p1<-gheatmap(p, tested, offset = 0.25, width = 3, font.size=2.5, 
             colnames_angle=-90, colnames_position = "top", hjust=1) + 
  scale_fill_manual(values=heatmapcols)
p1 <-p1 + theme(legend.position="right")
p1
#ggsave(filename="output/edgeR-Genus.pdf", plot=p, width=8, height=10, units="in")


#edgeR for functional gene
meanfun <- colSums(fun.table)/nrow(fun.table)
fun.order.table<-(fun.table[ ,order(meanfun,decreasing=TRUE)])
fun.order.table[1:5,1:5]
class.ra.table <- otu_table(class.phy.rel)
MRA <- rowSums(class.ra.table)/ncol(class.ra.table)
group <- tax_table(class.phy.rel)[,c(2,6)]
class.mra.table <- data.frame(group,MRA)

#arrange the class table
library(tidyr)
class.mra.table$Genus <- paste(class.mra.table$Genus, rownames(class.mra.table), sep = '_')
class.mra.table <- class.mra.table %>% spread(Genus, MRA)
class.mra.table[is.na(class.mra.table)] <- 0
rownames(class.mra.table)<-class.mra.table$Phylum
class.mra.table<-as.matrix(t(class.mra.table[,-1])*100)
colsum <-apply(class.mra.table,2,sum)
rowsum<-apply(class.mra.table,1,sum)
top10phy_table<-(class.mra.table[order(rowsum,decreasing=TRUE),order(colsum,decreasing=TRUE)])[,1:10]
head(top10phy_table)
top10phy_table<-as.matrix(top10phy_table)



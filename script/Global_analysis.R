# loading packages
PKGs <- c('ggplot2', 'ape', 'Biostrings', 'vegan', 
          'reshape', 'cowplot', 'microbiome', 'RColorBrewer',
          'ggpubr', 'dplyr', 'phyloseq')

lapply(PKGs, require, character.only = TRUE, warn.conflicts = FALSE)

# set work dorectory
setwd('E:/thermokast_lakes/water_microbes/')
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


# Community composition
#melt to long format (for ggplot) 
table(tax_table(meta_physeq)[,"Phylum"])
# ##prune out class below 1% in each sample and prevalence lower than 10/100 at class level
# meta.com.cla <- microbiome::aggregate_rare(meta_physeq_rel, level = "Class", 
#                                            detection = 1/100, prevalence = 10/100)
# allmean <- rowMeans(otu_table(meta.com.cla))
# 
# regiontype <- as.factor(sample_data(meta.com.cla)$Region)
# table(regiontype)
# mean_region <- sapply(levels(regiontype),function(i){
#   rowMeans(otu_table(meta.com.cla)[,region == i])
# })
# rel_abun_dat_cla <- data.frame(mean_region, All_mean = allmean)
# rel_abun_dat_cla <- dplyr::arrange(rel_abun_dat_cla, desc(All_mean))
# rel_abun_dat_cla

##prune out class below 1% in each sample and prevalence lower than 10/100 at Order level
meta.com.ord <- microbiome::aggregate_rare(meta_physeq_rel, level = "Order", 
                                           detection = 1/100, prevalence = 10/100)
allmean <- rowMeans(otu_table(meta.com.ord))

regiontype <- as.factor(sample_data(meta.com.ord)$Region)
table(regiontype)

mean_region <- sapply(levels(regiontype), function(i){
  rowMeans(otu_table(subset_samples(meta.com.ord, region == i)))
})
mean_region <- data.frame(Pan_Arctic = rowMeans(otu_table(subset_samples(meta.com.ord, region == 'Pan-Arctic'))),Tibetan_Plateau = rowMeans(otu_table(subset_samples(meta.com.ord, region == 'Tibetan Plateau'))))



mean_region <- sapply(levels(regiontype),function(i){
  rowMeans(otu_table(meta.com.ord)[,region == i ])
})
rel_abun_dat_ord <- data.frame(mean_region, All_mean = allmean)
rel_abun_dat_ord <- dplyr::arrange(rel_abun_dat_ord, desc(All_mean))
rel_abun_dat_ord

# pie plot
rel_abun_dat_ord_plot <- data.frame(Order = rownames(rel_abun_dat_ord), prop = rel_abun_dat_ord$All_mean * 100)
rel_abun_dat_ord_plot$Order <- factor(rel_abun_dat_ord_plot$Order, ordered = T,
                                         levels = c(rel_abun_dat_ord_plot$Order[-2], 'Other'))
rel_abun_dat_ord_plot$ypos = cumsum(rel_abun_dat_ord_plot$prop)- 0.5 * rel_abun_dat_ord_plot$prop 

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 9
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
# Create a ggplot with 11 colors 
# Use scale_fill_manual
pie_for_taxa <- ggplot(rel_abun_dat_ord_plot, aes(x="", y=prop, fill=Order)) +
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) +
  scale_fill_manual(values = mycols) +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14)
  ) +
  geom_text(aes(y = ypos, label = scales::percent(prop/100),
            color = "white", size=6))
pie_for_taxa


##bar plot for the community composition for TP and PA
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


# alpha diversity
diversity <- estimate_richness(meta_physeq, measures = c("Chao1", 'Shannon', 'Simpson'))
meta_diversity <- cbind(diversity, Region = meta_table$Region, 
                        Sitegroup = meta_table$Sitegroup,
                        Site = meta_table$Site, MAT = metadata$MAT,
                        MAP = metadata$MAP, DOC = metadata$DOC,
                        SUVA254 = metadata$SUVA254, a320 = metadata$a320,
                        pH = metadata$pH)
meta_diversity %>% dplyr::select(c(1,3,4,5)) %>%
  group_by(Region) %>%
  dplyr::summarise_all(list(mean = mean, sd = sd, se = ~sd(./sqrt(.))))

# test the difference
par(mfrow = c(1, 3))
hist(log(meta_diversity$Chao1))
hist(meta_diversity$Shannon)
hist(log(meta_diversity$Simpson+1))
par(mfrow = c(1, 1))
library(lme4)
library(lmerTest)
library(multcomp)
mode1 <- lmer(log(Chao1) ~ Region + (1|Site), meta_diversity)
summary(mode1)
MuMIn::r.squaredGLMM(mode1)
mode2 <- lmer(Shannon ~ Region + (1|Site), meta_diversity)
summary(mode2)
MuMIn::r.squaredGLMM(mode2)
mode3 <- lmer(log(Simpson) ~ Region + (1|Site), meta_diversity)
summary(mode3)
MuMIn::r.squaredGLMM(mode3)
ggplot(meta_diversity, aes(x=MAT, y=Simpson)) + geom_point(alpha = 0.3)
mode4 <- lmer(log(Chao1) ~ MAT + (1|Site), meta_diversity)
summary(mode4)
MuMIn::r.squaredGLMM(mode4)
mode5 <- lmer(Shannon ~ MAT + (1|Site), meta_diversity)
summary(mode5)
MuMIn::r.squaredGLMM(mode5)
mode6 <- lmer(log(Simpson) ~ MAT + (1|Site), meta_diversity)
summary(mode6)
MuMIn::r.squaredGLMM(mode6)
ggplot(meta_diversity, aes(x=MAP, y=Chao1)) + geom_point(alpha = 0.3)
mode7 <- lmer(log(Chao1) ~ MAP + (1|Site), meta_diversity)
summary(mode7)
MuMIn::r.squaredGLMM(mode7)
mode8 <- lmer(Shannon ~ MAP + (1|Site), meta_diversity)
summary(mode8)
MuMIn::r.squaredGLMM(mode8)
mode9 <- lmer(log(Simpson) ~ MAP + (1|Site), meta_diversity)
summary(mode9)
MuMIn::r.squaredGLMM(mode9)
ggplot(meta_diversity, aes(x=DOC, y=Chao1)) + geom_point(alpha = 0.3)
mode10 <- lmer(log(Chao1) ~ DOC + (1|Site), meta_diversity)
summary(mode10)
MuMIn::r.squaredGLMM(mode10)
mode11 <- lmer(Shannon ~ DOC + (1|Site), meta_diversity)
summary(mode11)
MuMIn::r.squaredGLMM(mode11)
mode12 <- lmer(log(Simpson) ~ DOC + (1|Site), meta_diversity)
summary(mode12)
MuMIn::r.squaredGLMM(mode12)
mode13 <- lmer(log(Chao1) ~ SUVA254 + (1|Site), meta_diversity)
summary(mode13)
MuMIn::r.squaredGLMM(mode13)
mode14 <- lmer(Shannon ~ SUVA254 + (1|Site), meta_diversity)
summary(mode14)
MuMIn::r.squaredGLMM(mode14)
mode15 <- lmer(log(Simpson) ~ SUVA254 + (1|Site), meta_diversity)
summary(mode15)
MuMIn::r.squaredGLMM(mode15)
mode16 <- lmer(log(Chao1) ~ a320 + (1|Site), meta_diversity)
summary(mode16)
MuMIn::r.squaredGLMM(mode16)
mode17 <- lmer(Shannon ~ a320 + (1|Site), meta_diversity)
summary(mode17)
MuMIn::r.squaredGLMM(mode17)
mode18 <- lmer(log(Simpson) ~ a320 + (1|Site), meta_diversity)
summary(mode18)
MuMIn::r.squaredGLMM(mode18)
mode19 <- lmer(log(Chao1) ~ pH + (1|Site), meta_diversity)
summary(mode19)
MuMIn::r.squaredGLMM(mode19)
mode20 <- lmer(Shannon ~ pH + (1|Site), meta_diversity)
summary(mode20)
MuMIn::r.squaredGLMM(mode20)
mode21 <- lmer(log(Simpson) ~ pH + (1|Site), meta_diversity)
summary(mode21)
MuMIn::r.squaredGLMM(mode21)

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


#PERMANOVA analysis and NMDS plot
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

#diff class using microeco package
library(microeco)
meco_df <- phyloseq2meco(meta_physeq)
#calculate the abundance table
m1 <- meco_df$cal_abund()
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
C_cycle_process <- c("methanotrophy",  "methanogenesis", 
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
# test the significance of the functional taxa among the region with LMM 
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
#piechart for all samples
all_df <- melt_df %>% 
  dplyr::select(c(variable, value)) %>%
  group_by(variable) %>% 
  dplyr::summarize(value = sum(value)) %>%
  arrange(desc(value))  %>%
  mutate(prop = value / sum(melt_df$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
  mutate(variable=factor(variable, levels=variable))



## Define the colors you want
mycols <- c("#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", 
            "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", 
            "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64",
            "#599861", "#89C5DA")
# Use scale_fill_manual
pie_for_fun <- ggplot(all_df, aes(x="", y=prop, fill=reorder(variable, -ypos)))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) +
  geom_text(aes(x = 1.35, y = ypos, 
                label = paste0(round(prop, 1), '%',  sep= '')),
            color = "white", size=4) +
  scale_fill_manual('Cabon cycling processes', values = mycols) +
  guides(fill = guide_legend(reverse=T, ncol=1)) +
  theme_void() +
  theme(legend.position = "left")
pie_for_fun
pie_plot <- cowplot::plot_grid(pie_for_taxa, pie_for_fun,
                               labels = c('(a)', '(b)'), ncol = 2, 
                               label_x = .01, label_y = 1, 
                               hjust = 0, label_size = 14, align = "v")

#barplot
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
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5))+
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
carbon_fun_plot

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
res_abund <- data.frame(genus = sapply(strsplit(m3$res_abund$Taxa, "\\|"), "[[", 6),
                        m3$res_abund[, c(2:6)])
res_rf <- data.frame(genus = sapply(strsplit(m3$res_rf$Taxa, "\\|"), "[[", 6),
                     m3$res_rf[, c(2:3)])
top_genus <- res_rf$genus[1:10]

res_abund_plot <- res_abund %>% filter(genus %in% top_genus) %>% 
  mutate(Group = factor(Group, levels = c('Tibetan Plateau', 'Pan-Arctic'))) %>%
  ggplot(aes(x = genus, y = 100*Mean, fill = Group))+
  geom_bar(position = 'dodge', stat = 'identity', colour = 'black', width = 0.7)+
  scale_fill_manual(values=c('#1b9e77', '#d95f02'))+
  geom_errorbar(aes(ymin = 100*Mean, ymax = 100*Mean + 100*SE), width=.2,  position = position_dodge(0.7))+
  scale_x_discrete(limits = top_genus) +
  labs(x = 'Genus', y = 'Relative abundance (%)', fill = "Group") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +
  theme_bw()+
  theme(legend.position = c(0.15,0.9),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
res_abund_plot
library(grid)
vie <- viewport(width=0.5, height=0.45, x=0.68, y=0.75)
res_rf_plot <- res_rf %>% filter(genus %in% top_genus) %>% 
  mutate(genus = factor(genus, levels = top_genus)) %>%
  ggplot(aes(x = genus, y = MeanDecreaseGini)) +
  geom_bar(position="stack", stat="identity", alpha=0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20))+
  labs(x = NULL, y = "MeanDecreaseGini")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size=12, colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank())
print(res_rf_plot, vp = vie)
ggsave(file = "./meta_analysis/results/figs/res_rf_plot.pdf",
       width = 12, height = 10, units = 'in', device='pdf', dpi=300)
carbon_plot <- ggdraw() +
  draw_plot(carbon_fun_plot, x = 0, y = 0.45, width = 1, height = 0.55) +
  draw_plot(carbon_rf_plot, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("(a)", "(b)"), size = 14,
                  x = c(0, 0), y = c(1, 0.45))

carbon_plot


## test the relationship between LCBD and MAP across TP and PA
library(adespatial)
###total community
beta_meta_div <- beta.div(t(as.matrix(otu_table(meta_physeq))), 
                          method = "hellinger", sqrt.D = FALSE, samp = TRUE, 
                          nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)

env_div <- data.frame(LCBD = beta_meta_div$LCBD, sample_data(meta_physeq))

### calculate the average LCBD for each site
library(dplyr)
env_div_agg_meta <-  env_div %>% 
  dplyr::select(c(1, 5, 6, 7, 8, 19)) %>%
  group_by(Sitegroup, Site, Sitegroup1) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))

ggplot(env_div, aes(x=pH, y=LCBD)) + geom_point(alpha = 0.3)
lm1 <- lm(LCBD ~  pH + I(pH^2), data = env_div)
summary(lm1)
ggplot(env_div, aes(x=DOC, y=LCBD)) + geom_point(alpha = 0.3)
lm2 <- lm(LCBD ~  DOC, data = env_div)
summary(lm2)
ggplot(env_div, aes(x=a320, y=LCBD)) + geom_point(alpha = 0.3)
lm3 <- lm(LCBD ~  a320 + I(a320^2), data = env_div)
summary(lm3)
ggplot(env_div, aes(x=SUVA254, y=LCBD)) + geom_point(alpha = 0.3)
lm4 <- lm(LCBD ~  I(SUVA254^2), data = env_div)
summary(lm4)
ggplot(env_div, aes(x=log(MAT+15), y=LCBD)) + geom_point(alpha = 0.3)
lm5 <- lm(LCBD ~ log(MAT+15), data = env_div)
summary(lm5)
ggplot(env_div, aes(x=MAP, y=LCBD)) + geom_point(alpha = 0.3)
lm6 <- lm(LCBD ~  MAP + I(MAP^2), data = env_div)
summary(lm6)

##model selection
library(MASS)
library(glmulti)
A3 <- glmulti(LCBD ~ MAP + MAT + SUVA254 + a320 + DOC + pH, data = env_div,
              level=1, fitfunction=lm, crit="aicc", confsetsize= 2^13, plotty = F, trace = 0)
top_3 <- weightable(A3)
###  models with values more than 2 units away are considered substantially 
### less plausible than those with AICc values closer to that of the best model. 
### refrence:Anderson, D. R. (2007). Model based inference in the life sciences: A primer on evidence. New York: Springer. 
top_3 <- top[top$aicc <= min(top$aicc) + 2,] # 
top_3

modes_inf <- NULL
for(i in 1:nrow(top_3)){
  rse_sum <- summary(A3@objects[[i]])
  adj.r.squared <- rse_sum$adj.r.squared # obtain the adjust r squared
  multicollinearity <- any(car::vif(A3@objects[[i]]) > 2) # check the multicollinearity
  tmp <- data.frame(adj.r.squared, multicollinearity)
  if(is.null(modes_inf)){
    modes_inf<-tmp
  } else {
    modes_inf <- rbind(modes_inf,tmp)
  } 
}
modes_inf <- cbind(top_3, modes_inf)
modes_inf

vpa.mod <- varpart(env_div_agg$LCBD, ~ env_div_agg$HIX,
                   ~ env_div_agg$MAP)
plot(vpa.mod)




## random forest analysis
### Regression:
env_div_rf <- env_div %>%
  select(c('LCBD',  'MAT', 'MAP', 'DOC', 'SUVA254', 'a320', 'pH'))
library(randomForest)
library(rfPermute)
set.seed(131)
lcbd.rf <- rfPermute(LCBD ~ ., data=env_div_rf, ntree=1000, num.rep = 100,
                     importance=TRUE, na.action=na.omit)
impor.dat <-data.frame(variables = rownames(importance(lcbd.rf)), IncMSE=importance(lcbd.rf)[,1])

library(dplyr)
impor_meta_plot <- impor.dat %>% 
  arrange(IncMSE)  %>%
  mutate(variables=factor(variables, levels=variables)) %>%
  ggplot(aes(x = variables, y = IncMSE)) +
  geom_bar(stat="identity", fill="tomato3", colour=NA, width = 0.65) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 45)) +
  ylab('Increases in MSE %') + xlab(NULL) +
  theme_classic() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size =12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.6, 0.75)) +
  coord_flip()

## Show "importance" of variables: higher value mean more important:
mean(lcbd.rf$rf$rsq)
sd(lcbd.rf$rf$rsq)

meta_lcbd_map <- ggplot(env_div, aes(MAP, LCBD, col = Sitegroup)) +
  geom_smooth(method="lm", formula = y ~ poly(x, 2),
              size=1, se=T, colour='black') +
  geom_point(size=1, alpha=0.8) +
  xlab('MAP (mm)') + ylab('LCBD') +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.6, 0.75))

rf_plot <- cowplot::plot_grid(impor_meta_plot, meta_lcbd_map,
                              labels = c('(a)', '(b)'), ncol = 2, 
                              label_x = .01, label_y = 1, 
                              hjust = 0, label_size = 14, align = "v")

#ggsave(rf_plot, file = "./meta_analysis/results/figs/rf_plot.pdf",
#       width = 13, height = 6, units = 'in', device='pdf', dpi=300)

###carbon cycling community
beta_meta_c_div <- beta.div(t(as.matrix(otu_table(meta_carbon_phy))), 
                            method = "hellinger", sqrt.D = FALSE, samp = TRUE, 
                            nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)

env_c_div <- data.frame(LCBD = beta_meta_c_div$LCBD, sample_data(meta_physeq))

## random forest analysis
### Regression:
env_div_c_rf <- env_c_div %>%
  select(c('LCBD',  'MAT', 'MAP', 'DOC', 'SUVA254', 'a320', 'pH'))
library(randomForest)
library(rfPermute)
set.seed(131)
lcbd.c.rf <- rfPermute(LCBD ~ ., data=env_div_c_rf, ntree=1000, num.rep = 100,
                       importance=TRUE, na.action=na.omit)
impor.c.dat <-data.frame(variables = rownames(importance(lcbd.c.rf)), IncMSE=importance(lcbd.c.rf)[,1])

library(dplyr)
impor_meta_c_plot <- impor.c.dat %>% 
  arrange(IncMSE)  %>%
  mutate(variables=factor(variables, levels=variables)) %>%
  ggplot(aes(x = variables, y = IncMSE)) +
  geom_bar(stat="identity", fill="tomato3", colour=NA, width = 0.65) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 45)) +
  ylab('Increases in MSE %') + xlab(NULL) +
  theme_classic() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size =12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.6, 0.75)) +
  coord_flip()

## Show "importance" of variables: higher value mean more important:
mean(lcbd.c.rf$rf$rsq)
sd(lcbd.c.rf$rf$rsq)

meta_c_lcbd_map <- ggplot(env_c_div, aes(MAP, LCBD, col = Sitegroup)) +
  geom_smooth(method="lm", formula = y ~ poly(x, 2),
              size=1, se=T, colour='black') +
  geom_point(size=1, alpha=0.8) +
  xlab('MAP (mm)') + ylab('LCBD') +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.6, 0.75))

rf_c_plot <- cowplot::plot_grid(impor_meta_c_plot, meta_c_lcbd_map,
                                labels = c('(a)', '(b)'), ncol = 2, 
                                label_x = .01, label_y = 1, 
                                hjust = 0, label_size = 14, align = "v")

ggsave(rf_c_plot, file = "./meta_analysis/results/figs/rf_c_plot.pdf",
       width = 13, height = 6, units = 'in', device='pdf', dpi=300)

ggplot(env_c_div, aes(x=pH, y=LCBD)) + geom_point(alpha = 0.3)
lm1 <- lm(LCBD ~  pH + I(pH^2), data = env_c_div)
summary(lm1)
ggplot(env_c_div, aes(x=DOC, y=LCBD)) + geom_point(alpha = 0.3)
lm2 <- lm(LCBD ~  DOC, data = env_c_div)
summary(lm2)
ggplot(env_c_div, aes(x=a320, y=LCBD)) + geom_point(alpha = 0.3)
lm3 <- lm(LCBD ~  a320 + I(a320^2), data = env_c_div)
summary(lm3)
ggplot(env_c_div, aes(x=SUVA254, y=LCBD)) + geom_point(alpha = 0.3)
lm4 <- lm(LCBD ~  SUVA254 + I(SUVA254^2), data = env_c_div)
summary(lm4)
ggplot(env_c_div, aes(x=log(MAT+15), y=LCBD)) + geom_point(alpha = 0.3)
lm5 <- lm(LCBD ~ log(MAT+15), data = env_c_div)
summary(lm5)
ggplot(env_c_div, aes(x=MAP, y=LCBD)) + geom_point(alpha = 0.3)
lm6 <- lm(LCBD ~  MAP + I(MAP^2), data = env_c_div)
summary(lm6)


## partial mantel test function for Tibatan Plateau and Pan Arctic
partial.mantel.fun <- function(phylo) {
  env.table <- data.frame(sample_data(phylo))
  #env.table <- env.table[complete.cases(env.table), ]
  otu_table <- as.matrix(t(otu_table(phylo)))
  otu_table_hel <- decostand(otu_table, 'hellinger')
  otu_table_hel_dist <- vegdist(otu_table_hel, 'bray',upper=F)
  df <- NULL
  vars <- c("DOC", "SUVA254", "a320", "MAP", "MAT", "pH")
  for (x in vars) {
    x.dist <- vegdist(scale(env.table[,x]), 'euclidean', na.rm = T)
    z.dist <- vegdist(scale(env.table[ , setdiff(vars, x)]), 'euclidean', na.rm = T)
    mode <- mantel.partial(x.dist, otu_table_hel_dist, z.dist, 
                           method = "pearson", permutations = 999, na.rm = T)
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

all.total.commun.par.mant <- partial.mantel.fun(meta_physeq)
all.carbon.commun.par.mant <- partial.mantel.fun(meta_carbon_phy)
## PLOT
## devtools::install_github('hannet91/ggcor')
library(ggcor)
set.seed(123456)

all.par.man.tibble <- tibble(spec = c(rep('all.total.commun.par.mant', nrow(all.total.commun.par.mant)), 
                                  rep('all.carbon.commun.par.mant', nrow(all.carbon.commun.par.mant))), 
                         rbind(all.total.commun.par.mant, all.carbon.commun.par.mant))
#par.man.tibble <- tibble(spec = c(rep('Total community composition', nrow(total.commun.par.mant))), 
#                         total.commun.par.mant)
vars <- c("MAT", "MAP", "DOC", "SUVA254", "a320", "pH")
env.table <- sample_data(meta_physeq)[ , vars]

mantel02 <- all.par.man.tibble %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.3, 0.4, Inf), 
                 labels = c("<0.30", "0.30-0.40", ">0.40"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"),
                       right = T))
quickcor(env.table[complete.cases(env.table),], type = "upper") + geom_square() + 
  add_link(mantel02, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE) +
  scale_color_manual(values = c('#d95f02', 'grey', '#1b9e77', '#3C5488FF')) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  geom_diag_label() + remove_axis("x")


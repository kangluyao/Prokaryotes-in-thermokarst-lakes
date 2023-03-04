# test the difference of the environmental factors between TP and PA
##standard error function
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
meta_env_sum <- metadata %>% dplyr::select(c('Region', 'DOC', 'SUVA254', 'a320',
                                             'MAP', 'MAT', 'pH')) %>% group_by(Region) %>%
  dplyr::summarise_all(list(mean = mean, sd = sd, se = stderr), na.rm = TRUE)
# write.csv(meta_env, file = 'E:/thermokast_lakes/water_microbes/meta_analysis/results/tables/meta_env.csv')
# test the difference
library(lme4)
library(lmerTest)
library(multcomp)
par(mfrow = c(2, 3))
hist(metadata$DOC)
hist(metadata$SUVA254)
hist(metadata$a320)
hist(metadata$MAP)
hist(metadata$MAT)
hist(metadata$pH)
par(mfrow = c(1, 1))

par(mfrow = c(2, 3))
hist(log(metadata$DOC))
hist(log(metadata$SUVA254))
hist(log(metadata$a320))
hist(log(metadata$MAP))
hist(log(metadata$MAT+15))
hist(log(metadata$pH))
par(mfrow = c(1, 1))

env_mode1 <- lmer(log(DOC) ~ Region + (1|Site), metadata)
summary(env_mode1)
env_mode2 <- lmer(log(SUVA254) ~ Region + (1|Site), metadata)
summary(env_mode2)
env_mode3 <- lmer(log(a320) ~ Region + (1|Site), metadata)
summary(env_mode3)
env_mode4 <- lmer(log(MAP) ~ Region + (1|Site), metadata)
summary(env_mode4)
env_mode5 <- lmer(MAT ~ Region + (1|Site), metadata)
summary(env_mode5)
env_mode6 <- lmer(log(pH) ~ Region + (1|Site), metadata)
summary(env_mode6)

par(mfrow = c(2, 3))
qqnorm(resid(env_mode1))
qqnorm(resid(env_mode2))
qqnorm(resid(env_mode3))
qqnorm(resid(env_mode4))
qqnorm(resid(env_mode5))
qqnorm(resid(env_mode6))
par(mfrow = c(1, 1))

# test the difference of alpha diversity between PA and TP
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

#test the difference with Wilcoxon Test
par(mfrow = c(1, 3))
hist(meta_diversity$Chao1)
hist(meta_diversity$Shannon)
hist(meta_diversity$Simpson)
par(mfrow = c(1, 1))

#Wilcoxon Test
wilcox.test(Chao1 ~ Region, data = meta_diversity,
            exact = FALSE)
wilcox.test(Shannon ~ Region, data = meta_diversity,
            exact = FALSE)
wilcox.test(Simpson ~ Region, data = meta_diversity,
            exact = FALSE)

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

# test the difference of the community structure between TP and PA
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
NMDS_plot
metadata <- as(sample_data(meta_physeq), "data.frame")
adonis2(phyloseq::distance(meta_physeq, method = "bray") ~ Region,
        data = metadata)

#diff taxa using LEfSe with microeco package
library(microeco)
library(file2meco)
meco_df <- file2meco::phyloseq2meco(meta_physeq)
#calculate the abundance table
m1 <- meco_df$cal_abund()
ps_genus <- phyloseq::tax_glom(meta_physeq, taxrank = 'Genus')
meco_genus_df <- phyloseq2meco(ps_genus)
#calculate the abundance table
m1_genus <- meco_genus_df$cal_abund()

m1_genus <- trans_diff$new(dataset = meco_genus_df, method = "lefse", 
                           group = "Region", alpha = 0.01, 
                           lefse_subgroup = NULL)
#t1$res_lefse is the LEfSe result
#t1$res_abund is the abundance information
lefse_plot <- m1_genus$plot_diff_bar(use_number = 1:37, LDA_score = 4, color_values = c('#d95f02', '#1b9e77'), 
                        plot_vertical = TRUE, width = 0.5)
#m1_genus$plot_diff_abund(use_number = 1:37, color_values = c('#d95f02', '#1b9e77'))

# Combining Plots
library(cowplot)
comparison_plot <- ggdraw() +
  draw_plot(NMDS_plot, x = 0, y = 1/3, width = 0.4, height = 2/3) +
  draw_plot(lefse_plot, x = 0.4, y = 0, width = 0.6, height = 1) +
  draw_plot(alpha_region_plot, x = 0, y = 0, width = 0.4, height = 1/3) +
  draw_plot_label(label = c("A", "B", "C"), size = 14,
                  x = c(0, 0.4, 0), y = c(1, 1, 1/3))

comparison_plot


#we can format the data for Lefse analysis (http://huttenhower.sph.harvard.edu/galaxy)
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


# network analysis
library(phyloseq)
#pruned the OTUs with the mean relative abundance less than 0.005%
physeqr = phyloseq::transform_sample_counts(meta_physeq, function(x) x / sum(x))
physeqrF = filter_taxa(physeqr, function(x) mean(x) < .005/100, TRUE)
rmtaxa = taxa_names(physeqrF)
alltaxa = taxa_names(meta_physeq)
myTaxa = alltaxa[!alltaxa %in% rmtaxa]
physeqaF <- prune_taxa(myTaxa,meta_physeq)
physeqaF

#divided the filtered OTU table into two groups (i.e., high-altitude vs. high-latitude)
pa_phylo_net <- subset_samples(physeqaF, Region == "Pan-Arctic")
pa_phylo_net <- prune_taxa(taxa_sums(pa_phylo_net) > 0, pa_phylo_net) 
tp_phylo_net <- subset_samples(physeqaF, Region == "Tibetan Plateau")
tp_phylo_net <- prune_taxa(taxa_sums(tp_phylo_net) > 0, tp_phylo_net) 

pa_phylo_net
tp_phylo_net

#randomly resample the samples from alpine thermokarst lakes 
set.seed(1234) #To ensure reproducibility, the seed was set as 1234
rand.sample <- sample(sample_data(tp_phylo_net)$Sample_Name, 118, replace = F)
tp_phylo_rand_net <- subset_samples(tp_phylo_net, Sample_Name %in% rand.sample)
tp_phylo_rand_net <- prune_taxa(taxa_sums(tp_phylo_rand_net) > 0, tp_phylo_rand_net)
tp_phylo_rand_net

comm.pa.net <- otu_table(pa_phylo_net)
comm.tp.net <- otu_table(tp_phylo_rand_net)

#prepare the format of input data for Molecular Ecological Network Analyses (MENA) pipeline (http://ieg4.rccc.ou.edu/mena)
comm.pa.net[comm.pa.net == 0] <- ""
comm.tp.net[comm.tp.net == 0] <- ""


write.table(comm.pa.net, file = 'E:/thermokast_lakes/water_microbes/meta_analysis/results/tables/networks_6.30/comm_pa_net.txt')
write.table(comm.tp.net, file = 'E:/thermokast_lakes/water_microbes/meta_analysis/results/tables/networks_6.30/comm_tp_net.txt')
#then import the otu table into MENA pipeline

#To compare the difference in empirical network 
#indices between high-altitude and high-latitude 
#thermokarst lakes, the student t-test was employed 
#using the standard deviations derived from corresponding random networks
library(BSDA)
#Average clustering coefficient
tsum.test(mean.x = 0.643, s.x = 0.006, n.x = 118,
          mean.y = 0.653, s.y = 0.003, n.y = 118)

#Average path distance
tsum.test(mean.x = 2.692, s.x = 0.018, n.x = 118,
          mean.y = 5.080, s.y = 0.015, n.y = 118)

tsum.test(mean.x = 0.544, s.x = 0.003, n.x = 118,
          mean.y = 0.868, s.y = 0.004, n.y = 118)

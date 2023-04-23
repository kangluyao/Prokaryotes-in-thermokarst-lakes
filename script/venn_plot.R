##################################################################
# unique otus profile among two regions
# Load the library
library(limma)
library(phyloseq)
library(VennDiagram)
otu <- data.frame(otu_table(meta_physeq))
# List of items
x <- list(otu %>%
            mutate(rowsum = rowSums(otu[, c(1:108)])) %>%
            filter(rowsum > 0) %>%
            rownames(), 
          otu %>% data.frame() %>%
            mutate(rowsum = rowSums(otu[, c(109:306)])) %>%
            filter(rowsum > 0) %>%
            rownames())
names(x) <- c("Pan-Arctic","Tibetan Plateau")
# plot
# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
# Further customization
display_venn(
  x,
  category.names = c("Pan-Arctic", "Tibetan Plateau"),
  disable.logging = TRUE,
  scaled = F,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#E69F00", "#009E73"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.default.pos = "outer",
  cat.pos = c(-0.5, 0.5),
  cat.dist = c(0.045, 0.045)
)

ggvenn(
  x, c("Pan-Arctic", "Tibetan Plateau"),
  fill_color = c("#E69F00", "#009E73"),
  stroke_size = 0.5, set_name_size = 4
)

#
# Generate example data
pa_venn <- otu %>%
  mutate(rowsum = rowSums(otu[, c(1:108)])) %>%
  filter(rowsum > 0) %>%
  rownames()
tp_venn <- otu %>% data.frame() %>%
  mutate(rowsum = rowSums(otu[, c(109:306)])) %>%
  filter(rowsum > 0) %>%
  rownames()

# What are the possible letters in the universe?
universe <- sort(unique(c(pa_venn, tp_venn)))

# Generate a matrix, with the sets in columns and possible letters on rows
Counts <- matrix(0, nrow = length(universe), ncol = 2)
# Populate the said matrix
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% pa_venn
  Counts[i,2] <- universe[i] %in% tp_venn
}

# Name the columns with the sample names
colnames(Counts) <- c("PA","TP")
rownames(Counts) <- universe

share_otu <- Counts %>% data.frame() %>% filter(PA == 1 & TP == 1) %>% rownames()
pa_otu <- Counts %>% data.frame() %>% filter(PA == 1 & TP == 0) %>% rownames()
tp_otu <- Counts %>% data.frame() %>% filter(PA == 0 & TP == 1) %>% rownames()

share_reads <- sum(otu[row.names(otu) %in% share_otu, ])
pa_reads <- sum(otu[row.names(otu) %in% pa_otu, ])
tp_reads <- sum(otu[row.names(otu) %in% tp_otu, ])
otu_count_venn <- data.frame(group = c("High latitude", "High altitude", "Overlap"),
                             OTU_count = c(length(pa_otu), length(tp_otu), length(share_otu)))
otu_reads_venn <- data.frame(group = c("High latitude", "High altitude", "Overlap"),
                             reads_count = c(pa_reads, tp_reads, share_reads))
otu_count_venn_pro <- otu_count_venn %>%
  mutate(OTUs_prop = otu_count_venn$OTU_count/sum(otu_count_venn$OTU_count)*100) %>%
  mutate(OTUs_lab.ypos = cumsum(OTUs_prop) - 0.5*OTUs_prop)
otu_reads_venn_pro <- otu_reads_venn %>%
  mutate(reads_prop = otu_reads_venn$reads_count/sum(otu_reads_venn$reads_count)*100) %>%
  mutate(reads_lab.ypos = cumsum(reads_prop) - 0.5*reads_prop)

library(ggalluvial)
plot_df <- data.frame(otu_count_venn_pro[, c("group", "OTUs_prop")],
                      reads_prop = otu_reads_venn_pro[, "reads_prop"]) %>%
  dplyr::rename("Propertion of OTUs" = "OTUs_prop", "Relative abundance" = "reads_prop") %>%
  pivot_longer(cols = -c(group), names_to = "otu_reads", values_to = "value") %>%
  mutate(otu_reads = factor(otu_reads, levels = c("Propertion of OTUs", "Relative abundance"))) %>%
  mutate(group = factor(group, levels = c("High altitude", "Overlap", "High latitude"))) %>%
  ggplot(aes(y = value, x = otu_reads)) + 
  geom_flow(aes(alluvium = group), alpha = 0.9, lty = 2, fill = "white",
            color = "black", curve_type = "linear", width = 0.5) + 
  geom_col(aes(fill = group), width = 0.5, color = "black") + 
  labs(x = NULL, y = 'Propertion (%)', fill = 'group') +
  scale_fill_manual(values = c('#1b9e77', "#b3b75d", '#d95f02')) + 
  scale_y_continuous(expand = c(0, 0)) + 
  main_theme

## pie plot 
mycols <- c("#ffc15c", "#74d944", "#CE50CA", "#5e738f", "#C0717C", "#CBD5ec", "#5F7FC7",
"#00718b", "#d3d93e", "#169e77", "#d95f02", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
"#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
"#8A7C64", "#599861")

### Define the colors you want
pie_cols <- c("#b3b75d", '#1b9e77', '#d95f02')
                     
pie_for_otu_num_reads <- ggplot(otu_reads_venn_pro, 
                                aes(x = "", y = reads_prop, fill = reorder(group, -reads_lab.ypos))) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  geom_text(aes(x = 1.35, y = reads_lab.ypos, 
                label = paste0(reads_count, ' (', round(reads_prop, 1), '%', ')', sep= '')),
            color = "black", size = 3) +
  scale_fill_manual('group', values = pie_cols) +
  guides(fill = guide_legend(reverse = T)) +
  theme_void() +
  theme(legend.position = "left")
pie_for_otu_num_reads

# microbial composition of share, unique in TP and PA
share_otu <- Counts %>% data.frame() %>% filter(PA == 1 & TP == 1) %>% rownames()
pa_otu <- Counts %>% data.frame() %>% filter(PA == 1 & TP == 0) %>% rownames()
tp_otu <- Counts %>% data.frame() %>% filter(PA == 0 & TP == 1) %>% rownames()

phylo_uniq_tp <- subset_taxa(tp_phylo, OTU %in% tp_otu) %>% 
  tax_glom(taxrank = 'Order') %>%
  transform_sample_counts(function(x) x / sum(x))
  
phylo_share <- subset_taxa(meta_physeq, OTU %in% share_otu) %>% 
  tax_glom(taxrank = 'Order') %>%
  transform_sample_counts(function(x) x / sum(x))

phylo_uniq_pa <- subset_taxa(pa_phylo, OTU %in% pa_otu) %>% 
  tax_glom(taxrank = 'Order') %>%
  transform_sample_counts(function(x) x / sum(x))
phylo_uniq_pa <- prune_samples(sample_sums(phylo_uniq_pa) > 0, phylo_uniq_pa)
  
pieplot_fun <- function(phylo){
  data.frame(tax_table(phylo)[, 4], otu_table(phylo)) %>% 
    mutate(MRA = rowMeans(dplyr::select(., rownames(sample_data(phylo))))) %>%
    arrange(desc(MRA)) %>% dplyr::top_n(10, MRA) %>%
    dplyr::select(., -c('MRA')) %>% 
    bind_rows(summarise_all(., ~if(is.numeric(.)) 1-sum(.) else "Others")) %>%
    mutate(Order = factor(Order, levels = Order)) %>%
    pivot_longer(cols = -c(Order), names_to = "Sample_name", values_to = 'rel_abun') %>%
    dplyr::select(., -c('Sample_name')) %>% 
    group_by(Order) %>%
    dplyr::summarise(across(, mean, na.rm = TRUE)) %>%
    mutate(rel_abun = rel_abun * 100) %>%
    mutate(lab.ypos = cumsum(rel_abun) - 0.5*rel_abun) %>%
    ggplot(aes(x = "", y = rel_abun, fill = reorder(Order, -lab.ypos))) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start = 0) +
    geom_text(aes(x = 1.35, y = lab.ypos, 
                  label = paste0(round(rel_abun, 1), '%', sep= '')),
              color = "black", size = 3) +
    scale_fill_manual('Order', values = mycols) +
    guides(fill = guide_legend(reverse = T)) +
    theme_void() +
    theme(legend.position = "left")
}
pieplot_share <- pieplot_fun(phylo_share)
pieplot_tp <- pieplot_fun(phylo_uniq_tp)
pieplot_pa <- pieplot_fun(phylo_uniq_pa)

# Combining Plots
library(cowplot)
pie_plot <- ggdraw() +
  draw_plot(pie_for_otu_num_reads, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(pieplot_share, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(pieplot_tp, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(pieplot_pa, x = 0.5, y = 0, width = 0.5, height = 0.5)

pie_plot

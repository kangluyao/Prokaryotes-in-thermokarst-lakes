##################################################################
# unique otus profile among two regions
# Load the library
library(limma)
library(phyloseq)
library(VennDiagram)
otu <- data.frame(otu_table(meta_physeq))
# List of items
x <- list(otu %>% data.frame() %>%
            mutate(rowsum = rowSums(select(., c(1:108)))) %>%
            filter(rowsum > 0) %>%
            rownames(), 
          otu %>% data.frame() %>%
            mutate(rowsum = rowSums(select(., c(109:306)))) %>%
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
pa_venn <- otu %>% data.frame() %>%
  mutate(rowsum = rowSums(select(., c(1:108)))) %>%
  filter(rowsum > 0) %>%
  rownames(),
tp_venn <- otu %>% data.frame() %>%
  mutate(rowsum = rowSums(select(., c(109:306)))) %>%
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
shar_otu <- sum(otu[row.names(otu) %in% share_otu, ])
pa_otu <- sum(otu[row.names(otu) %in% pa_otu, ])
tp_otu <- sum(otu[row.names(otu) %in% tp_otu, ])
otu_count_venn <- data.frame(group = c("Pan-Arctic", "Tibetan Plateau", "Share_otu"),
                             value = c(pa_otu, tp_otu, shar_otu))
otu_count_venn_pro <- otu_count_venn %>%
  mutate(prop = otu_count_venn$value/sum(otu_count_venn$value)*100) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)


## pie plot 
### Define the colors you want
pie_cols <- c("#b3b75d", "#009E73", "#E69F00")
                     
pie_for_otu_num_reads <- ggplot(otu_count_venn_pro, 
                                aes(x = "", y = prop, fill = reorder(group, -lab.ypos))) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  geom_text(aes(x = 1.35, y = lab.ypos, 
                label = paste0(value, ' (', round(prop, 1), '%', ')', sep= '')),
            color = "black", size = 3) +
  scale_fill_manual('group', values = pie_cols) +
  guides(fill = guide_legend(reverse = T)) +
  theme_void() +
  theme(legend.position = "left")
pie_for_otu_num_reads

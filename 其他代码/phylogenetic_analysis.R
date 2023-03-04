library(ape)
library(ggplot2)
library("phytools") # for sims and ASRs
library("ggtree")

#extract the tree at Class rank 
x1.genus <- tax_glom(water.rel, taxrank="Genus")
gen_0.1per<- microbiome::core(x1.genus, detection = 1/100, prevalence = 10/188)

labels<-as.vector(tax_table(gen_0.1per)[,6])
x1 <- gen_0.1per

# add tip labels (making room with xlim first), node labels, background color, 
# branch colors (based on branch legths), and a legend for the branch colors
p3 <- ggtree(x1, branch.length="none") +
  xlim(0, 90) + 
  geom_tiplab(size=2, color="plum1",label=labels) +
  geom_label2(aes(subset=!isTip, label=node), size=2, color="darkred", alpha=0.5) +
  scale_color_continuous(low='white', high='hotpink', name="Branch length (my)") +
  theme(legend.position="bottom")
plot(p3)

# get the correlation matrix with taxa table and env factors
abund_table<-t(otu_table(x1))
cor.table<-cor(abund_table,env.table[ -1], use="complete.obs")
#We can use gheatmap to display the data matrix on the tips of the tree:

# basic plot and highlight clades
p5 <-ggtree(x1,branch.length="none") + 
  xlim(0, 125) +
  geom_tiplab(size=2, offset=86,label=labels) +
  #geom_hilight(node=124, fill="steelblue", alpha=0.5) +
  #geom_hilight(node=113, fill="darkgreen", alpha=0.5) +
  geom_hilight(node=47, fill="gray", alpha=0.5) +
  geom_hilight(node=49, fill="pink", alpha=0.5) +
  geom_hilight(node=75, fill="beige", alpha=0.5) +
  geom_hilight(node=53, fill="yellow", alpha=0.5) 

# add heatmap
p6 <-  gheatmap(p5, cor.table, offset=0.01, width=4, colnames_angle=-90, hjust=0, 
                colnames_position = "top", font.size=2)+
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C")

# plot
plot(p6)


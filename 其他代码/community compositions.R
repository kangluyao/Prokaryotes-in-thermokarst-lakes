setwd('E:/thermokast_lakes/water_microbes/tibet_dada2_asv/')
#load libraries
library(phyloseq)
library(ape)
library(Biostrings)
library(reshape)
library(ggplot2)

#read in otu table
otu.table <- read.csv('./data/total_taxa_data/otu.table.csv',sep=",", row.names=1)
otu.table <- as.matrix(otu.table)

#read in taxonomy
#seperated by kingdom phylum class order family genus species 
taxonomy <- read.csv('./data/total_taxa_data/taxonomy.csv',sep=",",row.names=1)
taxonomy <- as.matrix(taxonomy)

#read in metadata
metadata <- read.csv("E:/thermokast_lakes/water_microbes/tibet_dada2_asv/data/metadata.csv", row.names=1, header = T)

# read in tree
total.tree <- read_tree('./data/total_taxa_data/tree.nwk')

#read in represent dna sequences
ref_seqs <- readDNAStringSet(file = "./data/total_taxa_data/ref.seqs.fasta",
                             format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

#import as phyloseq objects
otu.table <- otu_table(otu.table, taxa_are_rows = TRUE)
taxonomy <- tax_table(taxonomy)
metadata <- sample_data(metadata)

#merge into one phyloseq object
water_physeq <- phyloseq(otu.table, taxonomy, metadata, total.tree, ref_seqs)
water_physeq
water.rel <- microbiome::transform(water_physeq, "compositional")


#phylum.boxplot
phylum_table <- otu_table((tax_glom(water.rel, 'Phylum')))
phylum.names <- as.vector(tax_table(tax_glom(water.rel, 'Phylum'))[,2])
rownames(phylum_table) <- phylum.names

rowmean <-sapply(1:nrow(phylum_table),function(x) mean(phylum_table[x,]))
phylum_table<-phylum_table[order(rowmean,decreasing=TRUE), ]

taxa_list<-rownames(phylum_table)[1:15]#Extract list of top N Taxa
new_df<-phylum_table[rownames(phylum_table) %in% taxa_list,]
new_df <- melt(new_df, id.vars = 'taxa')
colnames(new_df) <- c('taxa','site','relative_abundance')
new_df$taxa <- factor(new_df$taxa, ordered = T, levels = rev(taxa_list))

phylum.boxplot <- ggplot(new_df, aes(x = taxa, y = relative_abundance *100, fill = taxa)) +
  geom_boxplot(outlier.shape = NA) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="steelblue", size=1.4, alpha=0.2) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")  + ylab('Relative abundance (%)')+
  coord_flip()
phylum.boxplot  

#determine the composition within Proteobacteria
Proteobacteria.phy <- subset_taxa(water.rel, Phylum == 'Proteobacteria')
class.Prote.table <- otu_table((tax_glom(Proteobacteria.phy, 'Class')))
class.Prote.names <- as.vector(tax_table(tax_glom(Proteobacteria.phy, 'Class'))[,3])
rownames(class.Prote.table) <- class.Prote.names

rowmean <-sapply(1:nrow(class.Prote.table),function(x) mean(class.Prote.table[x,]))
class.Prote.table<-class.Prote.table[order(rowmean,decreasing=TRUE), ]

taxa_list<-rownames(class.Prote.table)[1:8]#Extract list of top N Taxa
new_df<-class.Prote.table[rownames(class.Prote.table) %in% taxa_list,]
new_df <- melt(new_df, id.vars = 'taxa')
colnames(new_df) <- c('taxa','site','relative_abundance')
new_df$taxa <- factor(new_df$taxa, ordered = T, levels = rev(taxa_list))

class.boxplot <- ggplot(new_df, aes(x = taxa, y = relative_abundance *100, fill = taxa)) +
  geom_boxplot(outlier.shape = NA) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="steelblue", size=1.4, alpha=0.2) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")  + ylab('Relative abundance (%)')+
  coord_flip()
class.boxplot 

##determine the class compositions within top 10 phylums##
class.phy <- tax_glom(water_physeq, taxrank="Class")
class.phy.rel <- microbiome::transform(class.phy, "compositional")

class.ra.table <- otu_table(class.phy.rel)
MRA <- rowMeans(class.ra.table)
group <- tax_table(class.phy.rel)[,c(2,3)]
class.mra.table <- data.frame(group,MRA)

#arrange the class table
library(tidyr)
class.mra.table <- class.mra.table %>% spread(Class, MRA)
class.mra.table[is.na(class.mra.table)] <- 0
rownames(class.mra.table)<-class.mra.table$Phylum
class.mra.table<-as.matrix(t(class.mra.table[,-1])*100)
colsum <-apply(class.mra.table,2,sum)
rowsum<-apply(class.mra.table,1,sum)
top10phy_table<-(class.mra.table[order(rowsum,decreasing=TRUE),order(colsum,decreasing=TRUE)])[,1:10]
head(top10phy_table)
top10phy_table<-as.matrix(top10phy_table)
# Get the stacked barplot
# create color palette:
#library(RColorBrewer)
#coul <- brewer.pal(nrow(top10_phy), "Pastel2") 
mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,
          558,43,652,31,610,477,588,99,81,503,562,76,96,495,77,12,90,
          345,255,401,366,276,158,436)
mycol <-colors()[rep(mycol,nrow(top10phy_table))]

#tiff(file="bar.order1.level.tiff",width=750,height=700,pointsize=15)
layout(matrix(1:2,2,1),heights=c(1.5:1))
opar <- par(no.readonly = T)
par(mar=c(2,5,2,2))
barplot(top10phy_table,width=1.8,space=0.4,plot=T,las=2,
        col=mycol[1:nrow(top10phy_table)],cex.axis=0.8,cex.names=0.7,border=NA,
        xlab = 'Phylum',ylab="Relative abundance(%)",
        offset=0,cex.lab=1)
par(mar=c(2,3.5,3.5,1))
plot.new()
legend("topleft",legend=rownames(top10phy_table),
       ncol=3,fill=mycol[1:nrow(top10phy_table)],cex=0.6,bty="n")
par(opar)









##Correlation.R
# ============================================================
# Tutorial on drawing a correlation map using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================
#Now calculate the correlation between individual Taxa and the environmental data
cor.mat.cal <- function(y, x, method){
  y <- as.matrix(y)
  x <- as.matrix(x)
  df<-NULL
  for(i in colnames(y)){
    for(j in colnames(x)){
      a <- y[, i, drop = F]
      b <- x[, j, drop = F]
      tmp <- c(i, j, cor(a[complete.cases(b), ], b[complete.cases(b), ],
                     use = "everything", method = method),
             cor.test(a[complete.cases(b), ], b[complete.cases(b), ], method = method)$p.value)
      if(is.null(df)){
        df <- tmp  
      }
      else{
        df <- rbind(df, tmp)
      }    
    }
  }
  df<-data.frame(row.names=NULL,df)
  colnames(df)<-c("Taxa","Env","Correlation","Pvalue")
  df$Pvalue<-as.numeric(as.character(df$Pvalue))
  df$AdjPvalue<-rep(0,dim(df)[1])
  df$Correlation<-as.numeric(as.character(df$Correlation))
  #You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
  # 1 -> donot adjust
  # 2 -> adjust Env + Type (column on the correlation plot)
  # 3 -> adjust Taxa + Type (row on the correlation plot for each type)
  # 4 -> adjust Taxa (row on the correlation plot)
  # 5 -> adjust Env (panel on the correlation plot)
  adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
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
    for(i in unique(df$Taxa)){
      for(j in unique(df$Type)){
        sel<-df$Taxa==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
  } else if (adjustment==4){
    for(i in unique(df$Taxa)){
      sel<-df$Taxa==i
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
  df$Taxa <-factor(df$Taxa, ordered = T, levels = rev(colnames(y)))
  df$Env <-factor(df$Env, ordered = T, levels = colnames(x))
  return(df)
}

##phylum Correlation.R
cor.phylum.env <- cor.mat.cal(t(phylum_table), sample_data(water_physeq)[,-(1:5)], 'pearson')
taxa_list <- colnames(top10phy_table)
sub.cor.phylum.env <- subset(cor.phylum.env, Taxa %in% taxa_list)

#####corplot  
ggplot(sub.cor.phylum.env, aes(y = Taxa, x = Env, colour = Taxa)) +
  geom_point(aes(size = abs(Correlation)))+ #,shape=1
  scale_colour_manual(values=c("#B1A4C0","#479E9B","#FDDC7B","#4169B2",
                               "#ACCDDC","#DD5F60","#F2B379","#7CB6DB",
                               "pink2","skyblue","Gold2","#009E73",
                               "#E69F00","dodgerblue","#E69F00"))+
  theme_bw(base_size = 12)+
  geom_text(aes(label=Significance), color="black", size = 5.5)+
  scale_size(range = c(1, 12),breaks=c(0.2,0.4,0.6,0.8), limits=c(0,0.8))+ 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  #theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = 'Environmental factors', y = 'Phylum',
       colour = "Phylum", size = "Pearson's r")+ #legend title
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
  

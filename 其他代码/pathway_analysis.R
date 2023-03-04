#function
fun.module_table <- read.table('./data/tax4fun2/Pred_module13/functional_prediction.txt', 
                        sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
fun.table <- t(fun.table[,-189])
fun.dis <- vegdist(fun.table)
#pathway 
path.module13.table <- read.csv('./data/tax4fun2/Pred_module13/9_categories_metabolism.csv', 
                         header = T, row.names = 1, stringsAsFactors = F)
path.module13.table <- path.module13.table[,!colnames(path.module13.table) %in% c('level1', 'level3')]

ds<-function(x)(c(sum=sum(x)))
dfm<-melt(path.module13.table,id.vars = 'level2')
path.module13.table<-cast(dfm,level2~variable,ds)
rownames(path.module13.table)<-path.module13.table$level2
path.module13.table <- path.module13.table[,!colnames(path.module13.table) %in% 'level2']
#otu <-read.table(file="order.xls",header=T,check.names=FALSE,sep="\t")
#rownames(data) <-sapply(rownames(data),function(x) gsub("o__",'',x,perl = TRUE))
#otu <-otu[,-1]
#al <- which(rownames(data) %in% c("All"))
#if(length(al)) data <-data[-al,]
rowsum <-sapply(1:nrow(path.module13.table),function(x) sum(path.module13.table[x,]))
path.module13.table<-path.module13.table[order(rowsum,decreasing=TRUE),]
dat <-sapply(1:ncol(path.module13.table),function(x) path.module13.table[,x]/sum(path.module13.table[,x])) 
colnames(dat) <-colnames(path.module13.table)
rownames(dat) <-rownames(path.module13.table)
lab <-rownames(dat)
mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,558,43,652,31,610,477,588,99,81,503,562,76,96,495)
mycol <-colors()[rep(mycol,20)]
#tiff(file="bar.order1.level.tiff",width=750,height=700,pointsize=15)
opar <- par(no.readonly = T)
layout(matrix(1:2,2,1),heights=c(1.5:1))
par(mar=c(2,5,2,2))
barplot(dat*100,width=1.8,space=0.4,plot=T,las=2,
        col=mycol[1:nrow(dat)],cex.axis=0.8,cex.names=0.7,
        border=NA,ylab="Proportion (%)",offset=0,cex.lab=1)
par(mar=c(2,4,2,1))
plot.new()
legend("topleft",legend=rownames(dat),ncol=3,fill=mycol[1:nrow(dat)],cex=0.6, bty="n")
par(opar)

###pathway coposition and their correlation with environmental factors
path.table <- read.table('./data/tax4fun2/Pred_total/pathway_prediction.txt', 
                         sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
metabolism.table <- subset(path.table, level3 %in% 'Metabolism')
rownames(metabolism.table) <- metabolism.table$level1
metabolism.table <- metabolism.table[,!colnames(metabolism.table) %in% c('level1', 'level3')]
ds<-function(x)(c(sum=sum(x)))
dfm<-melt(metabolism.table,id.vars = 'level2')
metabolism.table<-cast(dfm,level2~variable,ds)
rownames(metabolism.table)<-metabolism.table$level2
metabolism.table <- metabolism.table[,!colnames(metabolism.table) %in% 'level2']
metabolism.table1 <- t(metabolism.table)
rownames(metabolism.table1) <- colnames(metabolism.table)
colnames(metabolism.table1) <- rownames(metabolism.table)
#pie plot for composition
pie.path.table <- data.frame(pathways = rownames(metabolism.table),
                             proportion = rowSums(metabolism.table)/sum(metabolism.table))

#cluster analysis
library(ggdendro)
betad <- vegdist(metabolism.table, method="bray")
betad <- as.matrix(betad, labels=TRUE)
rownames(betad) <- colnames(betad) <- rownames(metabolism.table)
betad <- as.dist(betad)
#Cluster the samples
hc <- hclust(betad)
hc_d <- dendro_data(as.dendrogram(hc))
hc_d$labels

##phylum Correlation.R
cor.metabolism.env <- cor.mat.cal(metabolism.table1, scale(env.table[,-(1:3)]), 'pearson')
cor.metabolism.env$Correlation[cor.metabolism.env$AdjPvalue >= 0.05] <- NA
cor.metabolism.env$Taxa <- factor(cor.metabolism.env$Taxa, ordered = T, levels = hc_d$labels$label)
#####corplot  
ggplot(cor.metabolism.env, aes(y = Taxa, x = Env, colour = Taxa)) +
  geom_point(aes(size = abs(Correlation)))+ #,shape=1
  scale_colour_manual(values=c("#B1A4C0","#479E9B","#FDDC7B","#4169B2",
                               "#ACCDDC","#DD5F60","#F2B379","#7CB6DB",
                               "pink2","skyblue","Gold2","#009E73",
                               "#E69F00"))+
  theme_bw(base_size = 12)+
  geom_text(aes(label=Significance), color="black", size = 5.5)+
  scale_size(range = c(1, 12),breaks=c(0.1,0.2,0.4,0.6), limits=c(0,0.6))+ 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  #theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = 'Environmental factors', y = 'Pathways',
       colour = "Pathway", size = "Pearson's r")+ #legend title
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













# Install pathview from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("pathview")
library(pathview)
load(url("http://genomedata.org/gen-viz-workshop/pathway_visualization/pathview_Data.RData"))
# View the hsa03430 pathway from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa03430", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa03430")


library(pathview)
metabolite.data <- data.frame(FC = sim.mol.data(mol.type = "cpd", nmol = 3000))
head(metabolite.data)

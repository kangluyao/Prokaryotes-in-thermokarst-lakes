library(phyloseq)
physeqr = phyloseq::transform_sample_counts(meta_physeq, function(x) x / sum(x))
physeqrF = filter_taxa(physeqr, function(x) mean(x) < .005/100, TRUE)

rmtaxa = taxa_names(physeqrF)
alltaxa = taxa_names(meta_physeq)

myTaxa = alltaxa[!alltaxa %in% rmtaxa]

physeqaF <- prune_taxa(myTaxa,meta_physeq)
physeqaF

pa_phylo_net <- subset_samples(physeqaF, Region == "Pan-Arctic")
pa_phylo_net <- prune_taxa(taxa_sums(pa_phylo_net) > 0, pa_phylo_net) 
tp_phylo_net <- subset_samples(physeqaF, Region == "Tibetan Plateau")
tp_phylo_net <- prune_taxa(taxa_sums(tp_phylo_net) > 0, tp_phylo_net) 
pa_phylo_net
tp_phylo_net


set.seed(1234)
rand.sample <- sample(sample_data(tp_phylo_net)$Sample_Name, 118, replace = F)
tp_phylo_rand_net <- subset_samples(tp_phylo_net, Sample_Name %in% rand.sample)
tp_phylo_rand_net <- prune_taxa(taxa_sums(tp_phylo_rand_net) > 0, tp_phylo_rand_net)
tp_phylo_rand_net

comm.pa.net <- otu_table(pa_phylo_net)
comm.tp.net <- otu_table(tp_phylo_rand_net)

write.csv(comm.pa.net, file = 'E:/thermokast_lakes/water_microbes/meta_analysis/results/tables/network/comm_pa_net.csv')
write.csv(comm.tp.net, file = 'E:/thermokast_lakes/water_microbes/meta_analysis/results/tables/network/comm_tp_net.csv')


# test the significance of network properties with random networks

# the folder saving the input files
wd <- 'E:/thermokast_lakes/water_microbes/meta_analysis/data1/networks'

cp.r.pa.file <- 'Galaxy7-[Whole_Network_matrix_0.74_RMT_result].tabular'
cp.r.tp.file <- 'Galaxy21-[Whole_Network_matrix_0.76_RMT_result].tabular'


# the folder to save the output. please change to a new folder even if you are just testing the example data.
save.wd="E:/thermokast_lakes/water_microbes/meta_analysis/results"
# if(!dir.exists(save.wd)){dir.create(save.wd)}

setwd(wd)
cp.r.pa <- t(read.table(cp.r.pa.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))
cp.r.tp <- t(read.table(cp.r.tp.file, header = TRUE, sep = "\t", row.names = 1,
                        as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                        check.names = FALSE))

comm.tp.net <- t(comm.tp[rownames(cp.r.tp),])
comm.meta.net <- t(comm.tp[rownames(cp.r.meta),])


rand_net_properties %>%
  group_by(Region) %>%
  dplyr::summarise_all(list(mean = mean, sd = sd, se = ~sd(./sqrt(.))))


vars <- colnames(rand_net_properties)[c(-1)]
mode <- lapply(vars, function(x) {
  lm(substitute(i~Region,list(i = as.name(x))), data = rand_net_properties)})

summary.model <- function(model){
  F.value <- anova(model)$'F value' [1]
  p.value <- anova(model)$'Pr(>F)' [1]
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
result_net <-data.frame(vars, df)
result_net


rand_net_region_plot <- rand_net_properties %>% 
  dplyr::select(c('Region', 'avgCC', 'GD', 'CB', 'Modularity')) %>%
  tidyr::gather(net_properties, value, -c('Region')) %>%
  ggplot(aes(x = Region, y = value, fill = Region)) +
  geom_violin(trim=T, width=0.5, aes(fill = Region), colour = "#000000") +
  scale_x_discrete(limits = c('TP', 'PA')) +
  scale_fill_manual(values= c('#d95f02', '#1b9e77')) +
  geom_boxplot(width=0.1, fill="white", colour = "#000000") +
  labs(x = NULL, y = 'Net properties') +
  facet_wrap(~net_properties, scales = 'free') +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12, 
                                   angle = 45, hjust = 1),
        panel.grid = element_blank(), legend.position = 'none')
rand_net_region_plot

library(igraph)
ig.pa <- graph_from_adjacency_matrix(cp.r.pa, mode="undirected", weighted=TRUE, diag=FALSE)
ig.tp <- graph_from_adjacency_matrix(cp.r.tp, mode="undirected", weighted=TRUE, diag=FALSE)

# calculate vulnerability of each node
source('https://raw.githubusercontent.com/Mengting-Maggie-Yuan/warming-network-complexity-stability/master/Fig3_and_S7.stability/info.centrality.R')
pa.node.vul<-info.centrality.vertex(ig.pa)
tp.node.vul<-info.centrality.vertex(ig.tp)
max(pa.node.vul)
max(tp.node.vul)

#
setwd(save.wd)
nodes_name <- rownames(cp.r.meta)
physeq_meta_net <-subset_taxa(physeqaF, OTU %in% nodes_name)
table(tax_table(physeq_meta_net)[,2])
write.csv(tax_table(physeq_meta_net)[, c(2, 4)],
          file = './tables/network/meta_nodes_atrr.csv')
nodes_name <- rownames(cp.r.pa)
physeq_pa_net <-subset_taxa(physeqaF, OTU %in% nodes_name)
write.csv(tax_table(physeq_pa_net)[, c(2, 4)],
          file = './tables/network/pa_nodes_atrr.csv')
nodes_name <- rownames(cp.r.tp)
physeq_tp_net <-subset_taxa(physeqaF, OTU %in% nodes_name)
write.csv(tax_table(physeq_tp_net)[, c(2, 4)],
          file = './tables/network/tp_nodes_atrr.csv')

# subnetwork
network.cal <- function(comm, igraph){
  library(igraph)
  library(funrar)
  dat <- matrix_to_stack(dat.net<-as.matrix(comm), "value", "site", "species")
  # Removal of empty rows
  dat <- dat[which(dat$value > 0), ]
  sites <- rownames(comm)
  df <- NULL
  for(i in sites){
    vids <- as.character(subset(dat,site == i)$species)
    tmp <- list(vids)
    if(is.null(df)) {df <- tmp} else { df<-append(df, tmp)} 
  }
  subgraphs <- lapply(df, function(x) induced.subgraph(igraph, as.character(unlist(x))))
  
  # 1计算每一个亚网络（单个样本）的拓扑参数,边数量 The size of the graph (number of edges)
  num.edges <- lapply(subgraphs, function(x) length(E(x)))
  num.edges <- as.vector(unlist(num.edges))
  
  
  # 2顶点数量 Order (number of vertices) of a graph
  degree <- lapply(subgraphs, function(x) length(V(x)))# length(diversity(igraph, weights = NULL, vids = V(igraph)))
  degree <- as.vector(unlist(num.vertices))
  
  
  # 3.connectance (complexity)网络中物种之间实际发生的相互作用数之和（连接数之和）占总的潜在相互作用数（连接数）的比例，
  #等同于网络的密度（edge_density）,可以反映网络的复杂程度
  Complexity <-lapply(subgraphs, function(x) edge_density(x, loops=FALSE))# 同 graph.density;loops如果为TRUE,允许自身环（self loops即A--A或B--B）的存在
  Complexity <- as.vector(unlist(Complexity))
  
  
  # 4平均度(Average degree)
  average.degree = lapply(subgraphs, function(x) mean(igraph::degree(x)))# 或者为2M/N,其中M 和N 分别表示网络的边数和节点数。
  average.degree<-as.vector(unlist(average.degree))
  
  
  # 5平均路径长度(Average path length)
  average.path.length = lapply(subgraphs, function(x) average.path.length(x)) # 同mean_distance(igraph) # mean_distance calculates the average path length in a graph
  average.path.length<-as.vector(unlist(average.path.length))
  
  # 6直径(Diameter)
  diameter = lapply(subgraphs, function(x) diameter(x, directed = FALSE, 
                                                    unconnected = TRUE, weights = NA))
  diameter<-as.vector(unlist(diameter))
  
  # 7群连通度 edge connectivity / group adhesion
  #edge.connectivity = lapply(subgraphs, function(x) edge_connectivity(x))
  #edge.connectivity<-as.vector(unlist(edge.connectivity))
  
  # 8.Clustering coefficient：分局域聚类系数和全局聚集系数，是反映网络中节点的紧密关系的参数，也称为传递性。
  #整个网络的全局聚集系数C表征了整个网络的平均的“成簇性质”。
  Clustering.coefficient <- lapply(subgraphs, function(x) transitivity(x, type = c('global')))
  Clustering.coefficient <- as.vector(unlist(Clustering.coefficient))
  
  # 9介数中心性(Betweenness centralization)
  Betweenness.centralization <- lapply(subgraphs, function(x) centralization.betweenness(x)$centralization)
  Betweenness.centralization <- as.vector(unlist(Betweenness.centralization))
  
  # 10.closeness centralization
  Closeness.centralization <- lapply(subgraphs, function(x) mean(closeness(x, vids = V(x),mode = 'all',normalized = T,weights=NA)))
  Closeness.centralization <- as.vector(unlist(Closeness.centralization))
  
  subgraphs.indexes <- data.frame(num.edges, degree, Complexity, average.degree, 
                                  average.path.length, diameter, Clustering.coefficient, 
                                  Betweenness.centralization, Closeness.centralization)
  rownames(subgraphs.indexes) = rownames(subgraphs.indexes)
  return(subgraphs.indexes)
}
pa_net_indexes <- network.cal(comm.pa.net, ig.pa)
tp_net_indexes <- network.cal(comm.tp.net, ig.tp)

net_indexes <- rbind(pa_net_indexes, tp_net_indexes)
nrow(net_indexes)
net_indexes <- cbind(sample_data(meta_physeq)[c(rownames(comm.pa.net), rownames(comm.tp.net)),c('Region', 'Site')], net_indexes)


# test the difference
library(lme4)
library(lmerTest)
library(multcomp)
mode1 <- lmer(Complexity ~ Region + (1|Site), net_indexes)
summary(mode1)
mode2 <- lmer(Clustering.coefficient ~ Region + (1|Site), net_indexes)
summary(mode2)
mode3 <- lmer(Closeness.centralization ~ Region + (1|Site), net_indexes)
summary(mode3)
mode4 <- lmer(Betweenness.centralization ~ Region + (1|Site), net_indexes)
summary(mode4)

# network indexes violin
melted <- melt(net_indexes[,c("Complexity", "Clustering.coefficient",
                              "Closeness.centralization", "Betweenness.centralization", "Region")], id.vars = c("Region"))
net_indexes_plot <- ggplot(melted, aes(x = Region, y = value, fill = Region)) +
  geom_violin(trim=T, width=0.5, aes(fill = Region), colour = "#000000") +
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_fill_manual(values= c('#1b9e77', '#d95f02')) +
  geom_boxplot(width=0.1, fill="white", colour = "#000000") +
  labs(x = NULL, y = 'Network proporties') +
  facet_wrap(~variable, scales = 'free') +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12, 
                                   angle = 45, hjust = 1),
        panel.grid = element_blank(), legend.position = 'none')
net_indexes_plot




library(BSDA)
tsum.test(mean.x=0.591,   s.x=0.005, n.x=132,
          mean.y=0.615, s.y=0.002, n.y=132)

tsum.test(mean.x=2.51,   s.x=0.018, n.x=132,
          mean.y=5.666, s.y=0.02, n.y=132)

tsum.test(mean.x=0.582,   s.x=0.003, n.x=132,
          mean.y=0.878, s.y=0.005, n.y=132)






















# Network attack analysis
# Libraries
library(SpiecEasi)
require(reshape2)
library(phyloseq)
library(igraph)
# Calculating the natural connectivity from adjacency matrix
ncc <- function(ig) {
  evals <- eigen(ig)$value
  nc <- log(mean(exp(evals)))
}

# Calculating the natural connectivity from adjacency matrix of a graph
natcon <- function(ig) {
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc <- log(mean(exp(evals)))
}

# Targeted attack ordered by betweenness
nc.attackbetweenness <- function(ig) {
  hubord <- order(rank(betweenness(ig)), decreasing=TRUE)
  sapply(1:round(vcount(ig)*.8), function(i) {
    ind <- hubord[1:i]
    tmp <- delete_vertices(ig, V(ig)$name[ind])
    natcon(tmp)
  })
}

# Targeted attack ordered by node degree.
nc.attackdegree <- function(ig) {
  hubord <- order(rank(degree(ig)), decreasing=TRUE)
  sapply(1:round(vcount(ig)*.8), function(i) {
    ind <- hubord[1:i]
    tmp <- delete_vertices(ig, V(ig)$name[ind])
    natcon(tmp)
  })
}


# Node removals
attack<-function (adj.mat, node.sup)
{
  n.nodes <- dim(adj.mat)[1]
  adj.mat[node.sup, ] <- rep(0, n.nodes)
  adj.mat[, node.sup] <- rep(0, n.nodes)
  nc<-ncc(adj.mat)
  list(new.mat = adj.mat, nc=nc)
}

pa.nc.deg.rmt <- nc.attackdegree(ig.pa)
pa.nc.bet.rmt <- nc.attackbetweenness(ig.pa)


tp.nc.deg.rmt <- nc.attackdegree(ig.tp)
tp.nc.bet.rmt <- nc.attackbetweenness(ig.tp)

pa.nc.deg.rmt <- cbind(rm_pro = c(1:length(pa.nc.deg.rmt))/vcount(ig.pa), nat.connet = pa.nc.deg.rmt)
pa.nc.bet.rmt <- cbind(rm_pro = c(1:length(pa.nc.bet.rmt))/vcount(ig.pa), nat.connet = pa.nc.bet.rmt)

tp.nc.deg.rmt <- cbind(rm_pro = c(1:length(tp.nc.deg.rmt))/vcount(ig.tp), nat.connet = tp.nc.deg.rmt)
tp.nc.bet.rmt <- cbind(rm_pro = c(1:length(tp.nc.bet.rmt))/vcount(ig.tp), nat.connet = tp.nc.bet.rmt)

setwd(save.wd)
write.table(pa.nc.deg.rmt, file="pa-deg-att-rmt.txt",
            sep="\t", quote=F, row.names = T, col.names = T)
write.table(pa.nc.bet.rmt, file="pa-bet-att-rmt.txt",
            sep="\t", quote=F, row.names = T, col.names = T)
write.table(tp.nc.deg.rmt, file="tp-deg-att-rmt.txt", 
            sep="\t", quote=F, row.names = T, col.names = T)
write.table(tp.nc.bet.rmt, file="tp-bet-att-rmt.txt",
            sep="\t", quote=F, row.names = T, col.names = T)



robut.df.rmt <- data.frame(rbind(cbind(Region = rep('Pan-Arctic', sum(nrow(pa.nc.deg.rmt), nrow(pa.nc.bet.rmt))),
                                   type = c(rep('Proportion of removed nodes', nrow(pa.nc.deg.rmt)), 
                                            rep('Proportion of removed betweenness', nrow(pa.nc.bet.rmt))),
                                   nat.connet = rbind(pa.nc.deg.rmt, pa.nc.bet.rmt)),
                             cbind(Region = rep('Tibetan Plateau', sum(nrow(tp.nc.deg.rmt), nrow(tp.nc.bet.rmt))),
                                   type = c(rep('Proportion of removed nodes', nrow(tp.nc.deg.rmt)), 
                                            rep('Proportion of removed betweenness', nrow(tp.nc.bet.rmt))),
                                   nat.connet = rbind(tp.nc.deg.rmt, tp.nc.bet.rmt))))
robut.df.rmt$Region <- factor(robut.df.rmt$Region, ordered = T, 
                          levels = c('Tibetan Plateau', 'Pan-Arctic'))
robut.df.rmt$type <- factor(robut.df.rmt$type, ordered = T, 
                        levels = c('Proportion of removed nodes', 'Proportion of removed betweenness'))

robut.df.rmt$rm_pro <- as.numeric(robut.df.rmt$rm_pro)
robut.df.rmt$nat.connet <- as.numeric(robut.df.rmt$nat.connet)

robutness_test_plot_rmt <- ggplot(robut.df.rmt, aes(x = rm_pro, y = nat.connet, colour = Region))+
  geom_point() +
  scale_color_manual(values= c('#1b9e77', '#d95f02')) +
  labs(x = NULL, y = 'Natural connectivity') +
  facet_wrap(~type, scales = 'free') +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12, 
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.9, 0.8),
        panel.grid = element_blank())
robutness_test_plot_rmt




# robustness test
rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

pa.Weighted.simu<-rmsimu(netRaw=cp.r.pa, rm.p.list=seq(0.05,1,by=0.05), sp.ra=comm.pa.net, abundance.weighted=T,nperm=100)
pa.Unweighted.simu<-rmsimu(netRaw=cp.r.pa, rm.p.list=seq(0.05,1,by=0.05), sp.ra=comm.pa.net, abundance.weighted=F,nperm=100)

tp.Weighted.simu<-rmsimu(netRaw=cp.r.tp, rm.p.list=seq(0.05,1,by=0.05), sp.ra=comm.tp.net, abundance.weighted=T,nperm=100)
tp.Unweighted.simu<-rmsimu(netRaw=cp.r.tp, rm.p.list=seq(0.05,1,by=0.05), sp.ra=comm.tp.net, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 year=rep(2014,40),treat=rep("Warmed",40))
dat1<-data.frame(Proportion.removed=seq(0.05,1,by=0.05),tp.Unweighted.simu)
dat2<-data.frame(Proportion.removed=seq(0.05,1,by=0.05),pa.Unweighted.simu)

currentdat<-dat1

write.csv(currentdat,"random_removal_result_Y14_W.csv")


##plot
library(ggplot2)

currentdat$year

ggplot(dat2, aes(x=Proportion.removed, y=remain.mean)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  #scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()
  #facet_wrap(~year, ncol=3)


ggplot(currentdat[currentdat$weighted=="unweighted",], aes(x=Proportion.removed, y=remain.mean, group=treat, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()+
  facet_wrap(~year, ncol=3)










meco_tp <- phyloseq2meco(tp_phylo_rand_net)
meco_pa <- phyloseq2meco(pa_phylo_net)

#calculate the abundance table
tp_1 <- meco_tp$cal_abund()
pa_1 <- meco_pa$cal_abund()

# trans_network class
## require WGCNA package
tp_1 <- trans_network$new(dataset = meco_tp, cal_cor = "WGCNA", taxa_level = "OTU", filter_thres = 0, cor_method = "spearman")
pa_1 <- trans_network$new(dataset = meco_pa, cal_cor = "WGCNA", taxa_level = "OTU", filter_thres = 0, cor_method = "spearman")


#  Optionally
# use arbitrary coefficient threshold to construct network
tp_1$cal_network(p_thres = 0.01, COR_cut = 0.6)
pa_1$cal_network(p_thres = 0.01, COR_cut = 0.6)

# add modules in the network
tp_1$cal_module()
pa_1$cal_module()

# save network
# open the gexf file using Gephi(https://gephi.org/)
# require rgexf package
tp_1$save_network(filepath = "./meta_analysis/results/figs/tp_network_0.6.gexf")
pa_1$save_network(filepath = "./meta_analysis/results/figs/pa_network_0.6.gexf")


pa_igraph <- pa_1$res_network
tp_igraph <- tp_1$res_network



#network analysis
net.cal<-function(phylo){
  library(WGCNA)
  library(impute)
  library(preprocessCore)
  library(igraph)
  x <- data.frame(t(otu_table(phylo)))
  cp<-corAndPvalue(x,
                   use = "pairwise.complete.obs",
                   alternative = c("two.sided")
  )
  #To get a matrix of the corresponding FDR, use:
  cp.r<-cp$cor# 取相关性矩阵R值
  cp.p <- apply(cp$p, 2, p.adjust, method = "fdr")# 取相关性矩阵调整p值
  
  # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
  cp.r[cp.p>0.01|abs(cp.r)<0.6] = 0
  igraph <- graph_from_adjacency_matrix(cp.r,mode="undirected",weighted=TRUE,diag=FALSE)
  return(igraph)
}

pa_igraph <- net.cal(pa_phylo)
tp_igraph <- net.cal(tp_phylo)

# subnetwork
network.cal <- function(phylo, igraph){
  library(igraph)
  library(funrar)
  df_x <- data.frame(t(otu_table(phylo)))
  dat <- matrix_to_stack(dat.net<-as.matrix(df_x), "value", "site", "species")
  # Removal of empty rows
  dat <- dat[which(dat$value > 0), ]
  sites <- rownames(df_x)
  df <- NULL
  for(i in sites){
    vids <- as.character(subset(dat,site == i)$species)
    tmp <- list(vids)
    if(is.null(df)) {df <- tmp} else { df<-append(df, tmp)} 
  }
  subgraphs <- lapply(df, function(x) induced.subgraph(igraph, as.character(unlist(x))))
  
  # 3.connectance (complexity)网络中物种之间实际发生的相互作用数之和（连接数之和）占总的潜在相互作用数（连接数）的比例，
  #等同于网络的密度（edge_density）,可以反映网络的复杂程度
  Complexity <-lapply(subgraphs, function(x) edge_density(x, loops=FALSE))# 同 graph.density;loops如果为TRUE,允许自身环（self loops即A--A或B--B）的存在
  Complexity <- as.vector(unlist(Complexity))
  
  # 8.Clustering coefficient：分局域聚类系数和全局聚集系数，是反映网络中节点的紧密关系的参数，也称为传递性。
  #整个网络的全局聚集系数C表征了整个网络的平均的“成簇性质”。
  Clustering.coefficient <- lapply(subgraphs, function(x) transitivity(x, type = c('global')))
  Clustering.coefficient <- as.vector(unlist(Clustering.coefficient))
  
  # 9介数中心性(Betweenness centralization)
  Betweenness.centralization <- lapply(subgraphs, function(x) centralization.betweenness(x)$centralization)
  Betweenness.centralization <- as.vector(unlist(Betweenness.centralization))
  
  # 13.closeness centralization
  Closeness.centralization <- lapply(subgraphs, function(x) mean(closeness(x, vids = V(x),mode = 'all',normalized = T,weights=NA)))
  Closeness.centralization <- as.vector(unlist(Closeness.centralization))
  
  subgraphs.indexes <- data.frame(Complexity, Clustering.coefficient,
                                  Betweenness.centralization, Closeness.centralization)
  rownames(subgraphs.indexes) = rownames(subgraphs.indexes)
  return(subgraphs.indexes)
}
pa_net_indexes <- network.cal(pa_phylo, pa_igraph)
tp_net_indexes <- network.cal(tp_phylo, tp_igraph)

net_indexes <- rbind(pa_net_indexes, tp_net_indexes)

net_indexes <- cbind(sample_data(meta_physeq)[,c('Region', 'Site')], net_indexes)


































library(microeco)
tp_phylo_0.1 <- microbiome::aggregate_rare(tp_phylo_rand, level = "OTU", 
                                           detection = 0, prevalence = 5/132)
pa_phylo_0.1 <- microbiome::aggregate_rare(pa_phylo, level = "OTU", 
                                           detection = 0, prevalence = 5/132)

tp_phylo_0.1
pa_phylo_0.1

meco_tp <- phyloseq2meco(tp_phylo_rand_net)
meco_pa <- phyloseq2meco(pa_phylo_0.1)

#calculate the abundance table
tp_1 <- meco_tp$cal_abund()
pa_1 <- meco_pa$cal_abund()

# trans_network class
## Use R base cor.test, slow
tp_1 <- trans_network$new(dataset = meco_tp, cal_cor = "base", taxa_level = "OTU", filter_thres = 0, cor_method = "spearman")
pa_1 <- trans_network$new(dataset = meco_pa, cal_cor = "base", taxa_level = "OTU", filter_thres = 0, cor_method = "spearman")
## return t1$res_cor_p list; one table: correlation; another: p value

## When the OTU number is large, use R WGCNA package to replace R base to calculate correlations
## require WGCNA package
tp_1 <- trans_network$new(dataset = meco_tp, cal_cor = "WGCNA", taxa_level = "OTU", filter_thres = 0, cor_method = "spearman")
pa_1 <- trans_network$new(dataset = meco_pa, cal_cor = "WGCNA", taxa_level = "OTU", filter_thres = 0, cor_method = "spearman")


## SparCC method, require SpiecEasi package
library(SpiecEasi)
## SparCC is very slow, so consider filtering more species with low abundance
tp_1 <- trans_network$new(dataset = meco_tp, cal_cor = "SparCC", taxa_level = "OTU", filter_thres = 0, SparCC_simu_num = 100)
pa_1 <- trans_network$new(dataset = meco_pa, cal_cor = "SparCC", taxa_level = "OTU", filter_thres = 0, SparCC_simu_num = 100)


# construct network; require igraph package
tp_1$cal_network(p_thres = 0.01, COR_optimization = TRUE)
pa_1$cal_network(p_thres = 0.01, COR_optimization = TRUE)
## return t1$res_network

#  Optionally
# use arbitrary coefficient threshold to construct network
tp_1$cal_network(p_thres = 0.01, COR_cut = 0.6)
pa_1$cal_network(p_thres = 0.01, COR_cut = 0.6)

# add modules in the network
tp_1$cal_module()
pa_1$cal_module()

# save network
# open the gexf file using Gephi(https://gephi.org/)
# require rgexf package
tp_1$save_network(filepath = "./meta_analysis/results/figs/tp_network_0.6.gexf")
pa_1$save_network(filepath = "./meta_analysis/results/figs/pa_network_0.6.gexf")


pa_igraph <- pa_1$res_network
tp_igraph <- tp_1$res_network

# extract the sub-network of sample 'S1'
samp_names <- colnames(meco_tp$otu_table)
subnets <- NULL
for (i in 1:length(samp_names)){
  tmp_sub <- tp_1$subset_network(node = meco_tp$otu_table[, 1, drop = FALSE] %>% .[.[, 1] !=0, , drop = FALSE] %>% rownames, rm_single = TRUE)
  if(is.null(subnets)) {
    subnets = list(tmp_sub)
  } else (
    subnets = append(subnets, tmp_sub)
  )
} 
names(subnets) <- samp_names

tmp_sub <- induced.subgraph(tp_igraph, node = meco_tp$otu_table[, 'Z1.H1', drop = FALSE] %>% .[.[, 1] !=0, , drop = FALSE] %>% rownames)



#determine the within-module degree z-score(Zi) and the participation coefficient (Pi)
z_cal <- function(net.obj, df_z = NULL){#net.obj indicates networks object
  require(igraph)
  edge.list <- as.data.frame(as_edgelist(net.obj, names = TRUE))
  colnames(edge.list) <- c('source', 'target')
  wtc <- cluster_fast_greedy(net.obj) #use fast greedy methods to extract networks structure
  m = max(wtc$membership)
  for (i in 1:m) {
    modul_nodes_id = wtc[[i]]
    num_ids <- length(modul_nodes_id)
    modul_edge_list <- edge.list[edge.list$source %in% modul_nodes_id & edge.list$target %in% modul_nodes_id, ]
    k <- NULL
    for (j in 1:num_ids) {
      tmp <- nrow(modul_edge_list[modul_edge_list$source %in% modul_nodes_id[j], ])
      if(is.null(k)){
        k = tmp
      }
      else{
        k = c(k, tmp)
      }
    }
    df_tmp = data.frame(ASV = modul_nodes_id, 
                        Module = rep(paste('Module', i, sep = '' ),
                                     num_ids), z = (k-mean(k))/sd(k))
    if(is.null(df_z)){
      df_z = df_tmp
    } else {
      df_z = rbind(df_z, df_tmp)
    }
  }
  return(df_z)
}


p_cal <- function(net.obj, df_p = NULL){#net.obj indicates networks object
  require(dplyr)
  require(igraph)
  edge.list <- as.data.frame(as_edgelist(net.obj, names = TRUE))
  colnames(edge.list) <- c('source', 'target')
  wtc <- cluster_fast_greedy(net.obj)
  m = max(wtc$membership)
  for (i in 1:m) {
    df <- NULL
    modul_nodes_id = wtc[[i]]
    num_ids <- length(modul_nodes_id)
    for (j in 1:m) {
      modul_edge_list <- edge.list[edge.list$source %in% wtc[[i]] & edge.list$target %in% wtc[[j]], ]
      k <- NULL
      for (ij in 1:num_ids) {
        tmp <- nrow(modul_edge_list[modul_edge_list$source %in% modul_nodes_id[ij], ])
        if(is.null(k)){
          k = tmp
        } else {
          k = c(k, tmp)
        }
      }
      names(k) <- paste('Module', j, sep = '')
      if(is.null(df)){
        df = k
      } else {
        df <- cbind(df, k)
      }
    }
    rownames(df) <- modul_nodes_id
    colnames(df) <- paste('Module', 1:m, sep='')
    if(is.null(df_p)){
      df_p <- df
    } else {
      df_p <- rbind(df_p, df)
    }
  }
  df_p <- data.frame(ASV = rownames(df_p), df_p)
  total.degree<-degree(net.obj)
  total.degree <- data.frame(ASV = names(total.degree), total.degree = total.degree)
  p_rev <- merge(df_p, total.degree, by = 'ASV', all.x = TRUE) %>% 
    mutate_each(funs((./total.degree)^2), starts_with("Module")) %>% 
    select(., starts_with('Module')) %>% 
    rowSums(.)
  p <- data.frame(ASV = rownames(df_p), p = 1- p_rev)
  p <- p[p$p != 1, ]
  return(p)
}



z <- z_cal(g)
p <- p_cal(g)

zp.table <- merge(z, p, all.y = TRUE)

#Z_P plot
assign_module_roles <- function(zp){
  zp <- na.omit(zp)
  zp$roles <- rep(0, dim(zp)[1])
  outdf <- NULL
  for(i in 1:dim(zp)[1]){
    df <- zp[i, ]
    if(df$z < 2.5){ #non hubs
      if(df$p < 0.05){
        df$roles <- "ultra peripheral"
      }
      else if(df$p < 0.620){
        df$roles <- "peripheral"
      }
      else if(df$p < 0.80){
        df$roles <- "non hub connector"
      }
      else{
        df$roles <- "non hub kinless"
      }
    }
    else { # module hubs
      if(df$p < 0.3){
        df$roles <- "provincial hub"
      }
      else if(df$p < 0.75){
        df$roles <- "connector hub"
      }
      else {
        df$roles <- "kinless hub"
      }
    }
    if(is.null(outdf)){outdf <- df}else{outdf <- rbind(outdf, df)}
  }
  return(outdf)
}

plot_roles <- function(node.roles, roles.colors=NULL){
  x1<- c(0, 0.05, 0.62, 0.8, 0, 0.30, 0.75)
  x2<- c(0.05, 0.62, 0.80, 1,  0.30, 0.75, 1)
  y1<- c(-Inf,-Inf, -Inf, -Inf,  2.5, 2.5, 2.5)
  y2 <- c(2.5,2.5, 2.5, 2.5, Inf, Inf, Inf)
  lab <- c("ultra peripheral", "peripheral", "non-hub connector", "non-hub kinless", "provincial", "hub connector","hub kinless")
  if(is.null(roles.colors)){roles.colors <- c("#E6E6FA", "#DCDCDC", "#F5FFFA", "#FAEBD7", "#EEE8AA", "#E0FFFF", "#F5F5DC")}
  p <- ggplot() + geom_rect(data=NULL, mapping=aes(xmin=x1, xmax=x2, ymin=y1,ymax=y2, fill=lab))
  p <- p + guides(fill=guide_legend(title="Topological roles"))
  p  <- p + scale_fill_manual(values = roles.colors)
  p <- p + geom_point(data=node.roles, aes(x=p, y=z,color=roles)) + theme_bw()
  p <- p + scale_x_continuous(expand = c(0, 0))
  p<-p+theme(strip.background = element_rect(fill = "white"),
             axis.title = element_text(size = 14, colour = 'black'),
             axis.text = element_text(size = 12, colour = 'black')) + 
    xlab("Participation Coefficient")+ 
    ylab(" Within-module connectivity z-score")
  return(p)
}

taxa.roles <- assign_module_roles(zp.table)
plot_roles(taxa.roles)

#extract PC1 for DOM, Climate and physicochemical factors
library(vegan)
dom_table<-decostand(env.table[ ,dom_vars], 'standardize')
climate_table<-decostand(env.table[ ,climate_vars], 'standardize')
physicochemical_table<-decostand(env.table[ ,physicochemical_vars], 'standardize')

sol<-rda(dom_table,trace = F)
sol$CA$eig[1]/sum(sol$CA$eig)
dom_PC1<-scores(sol)$sites[,1]
scores(sol)$species

sol<-rda(climate_table,trace = F)
sol$CA$eig[1]/sum(sol$CA$eig)
climate_PC1<-scores(sol)$sites[,1]
scores(sol)$species

sol<-rda(physicochemical_table,trace = F)
sol$CA$eig[1]/sum(sol$CA$eig)
physicochemical_PC1<-scores(sol)$sites[,1]
scores(sol)$species

#correlation between core otu taxa with PC1 for DOM, Climate and physicochemical 
pc1 <-data.frame(dom_PC1, climate_PC1, physicochemical_PC1)
map.s275.k <- data.frame(scale(env.table[, c('S275_295', 'MAP', 'K')]))
g.otu.table <- data.frame(t(otu_table(prune_taxa(names(V(g)), core.phy.rel))))
g.otu.table <- sqrt(g.otu.table)
cor.coretaxa.pc1 <- cor.mat.cal(g.otu.table, map.s275.k, 'spearman')
cor.coretaxa.pc1.pvalue0.05 <- subset(cor.coretaxa.pc1, AdjPvalue < 0.05)
nrow(cor.coretaxa.pc1.pvalue0.05)
nrow(cor.coretaxa.pc1)
cor.coretaxa.pc1.pvalue0.05 <- cor.coretaxa.pc1[cor.coretaxa.pc1$Taxa %in% unique(cor.coretaxa.pc1.pvalue0.05$Taxa), ]
cor.coretaxa.pc1.pvalue0.05 <- cor.coretaxa.pc1.pvalue0.05[ ,1:3]
cor.coretaxa.pc1.pvalue0.05$Correlation <- abs(cor.coretaxa.pc1.pvalue0.05$Correlation)

terrain.data <- reshape(cor.coretaxa.pc1.pvalue0.05, idvar = "Taxa", timevar = "Env", direction = "wide")


taxa.roles.module <- taxa.roles[ ,c('ASV', 'Module', 'roles')]
colnames(taxa.roles.module) <- c('Taxa', 'Module', 'roles')
terrain.data <- merge(terrain.data, taxa.roles.module, by = 'Taxa', all.x = T )
terrain.data[1:5,1:6]
colnames(terrain.data) <- c('OTU', 'S275_295','MAP', 'K', 'Module', 'roles')
terrain.data <- subset(terrain.data, OTU %in% tips.name)
terrain.data <- terrain.data[complete.cases(terrain.data), ]


library(ggtern)
ggtern(data = terrain.data, aes(S275_295, MAP, K, color = roles, shape = roles)) + 
  theme_rgbw()+ 
  geom_point()+
  scale_color_manual(values = mycol)+
  scale_shape_manual(values=c(1,rep(c(0:2,5:6,9:10,11:12,14), times=6)))+
  labs(x = "S275_295", y= "MAP", z= "K")

ggtern(data = terrain.data, mapping = aes(x = S275_295, y = MAP, z = K, 
                                          color = Module)) +
  scale_shape_discrete(name = "roles") +
  geom_point(alpha = 0.75, size = 1.2) +
  scale_color_manual(values = mycol)+
  theme_rgbw() + 
  theme_showsecondary() +
  theme_showarrows()
  
  
  labs(title = "USDA Textural Classification Chart",
         + fill = "Textural Class", color = "Textural Class")



taxa_roles<-c('ultra peripheral', 'peripheral', 'non hub connector', 'non hub kinless',
              'provincial hub', 'connector hub', 'kinless hub')
df <- NULL
for(i in taxa_roles){
  tmp.phy<- prune_taxa(as.vector((subset(taxa.roles, roles== i))$ASV), core.phy)
  tmp.table <- as(taxa_sums(otu_table(tmp.phy)), 'matrix')
  tmp.table <- data.frame(OTU = rownames(tmp.table), tmp.table[,1])
  if(is.null(df)){
    df <- tmp.table
  } else {
    df <- merge(df, tmp.table, by = 'OTU', all = T)
  }
  df[is.na(df)] <- 0
}
colnames(df) <- c('OTU', 'ultra peripheral', 'peripheral', 'non hub connector',
                  'non hub kinless', 'provincial hub', 'connector hub', 'kinless hub')

write.csv(df, './data/taxa_roles_data/taxa.roles.table.csv')
taxa.role.phy <- prune_taxa(as.vector(df$OTU), core.phy)
#write the tree, otu table and taxonomy table
taxa.roles.tree <- phy_tree(taxa.role.phy)
taxa.roles.taxonomy <- tax_table(taxa.role.phy)
taxa.roles.otu.table <- otu_table(taxa.role.phy)
taxa.roles.ref.seqs <- refseq(taxa.role.phy)

write.csv(taxa.roles.taxonomy, './data/taxa_roles_data/taxonomy.csv')
write.csv(taxa.roles.otu.table, './data/taxa_roles_data/otu.table.csv')
write.tree(taxa.roles.tree, './data/taxa_roles_data/tree.nwk')
writeXStringSet(taxa.roles.ref.seqs, './data/taxa_roles_data/ref.seqs.fasta', append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")


















tail(df)
#'Module roles
#' This function assigns roles to features identified in sub communities using
#' two metrics that is: within-module degree which measures how well a particular feature is connected to others in the same subcommunity (module.)
#' The second metric is among-module connectivity which measures how a feature is linked
#' to other modules in the network. Features are classified as ultra peripherals,
#' peripherals, provincial, connectors, kinless, module hubs, or non hubs.
#' @param comm_graph : Graph object returned by `co_occurence_network` function.
#' @examples
#' taxa.roles <- module.roles(g)
#' p <- plot_roles(taxa.roles)
#' print(p)
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @references Guimera, Roger, and Luis A Nunes Amaral. 2005. “Functional Cartography of Complex Metabolic Networks.”
#' Nature 433 (7028). NIH Public Access: 895.
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},  Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#' @export module.roles
#' @export plot_roles


module.roles <- function(comm_graph){
  td <- network_degree(comm_graph)
  wmd <- within_module_degree(comm_graph)
  z <- zscore(wmd)
  amd <- among_module_connectivity(comm_graph)
  pc <- participation_coeffiecient(amd, td)
  zp <- data.frame(z,pc)
  nod.roles <- assign_module_roles(zp)
  return(nod.roles)
}

# find total degree for each of the features in the graph

network_degree <- function(comm_graph){
  ki_total <-NULL
  net_degree <- igraph::degree(comm_graph)
  for(i in 1:length(V(comm_graph))){
    ki <- net_degree[i]
    tmp <- data.frame(taxa=names(ki), total_links=ki)
    if(is.null(ki_total)){ki_total<-tmp} else{ki_total <- rbind(ki_total, tmp)}
  }
  return(ki_total)
}


#compute within-module degree for each of the features

within_module_degree <- function(comm_graph){
  mods <- igraph::get.vertex.attribute(comm_graph, "module")
  vs <- as.list(V(comm_graph))
  modvs <- data.frame(taxon = names(vs), mod = mods)
  sg1 <- decompose.graph(comm_graph,mode="strong")
  df <- data.frame()
  for(mod in unique(modvs$mod)){
    mod_nodes <- subset(modvs$taxon,modvs$mod==mod)
    neighverts <- unique(unlist(sapply(sg1,FUN=function(s){if(any(V(s)$name %in% mod_nodes)) V(s)$name else NULL})))
    g3 <- induced.subgraph(graph=comm_graph,vids=neighverts)
    mod_degree <- igraph::degree(g3)
    for(i in mod_nodes){
      ki <- mod_degree[which(names(mod_degree)==i)]
      tmp <- data.frame(module=mod, taxa=names(ki), mod_links=ki)
      df <- rbind(df,tmp)
    }
  }
  return(df)
}
#calculate the degree (links) of each node to nodes in other modules.

among_module_connectivity <- function(comm_graph){
  mods <- igraph::get.vertex.attribute(comm_graph, "module")
  vs <- as.list(V(comm_graph))
  modvs <- data.frame("taxa"= names(vs), "mod"=mods)
  df <- data.frame()
  for(i in modvs$taxa){
    for(j in modvs$taxa){
      if(are_adjacent(graph=comm_graph, v1=i , v2=j)){
        mod <- subset(modvs$mod, modvs$taxa==j)
        tmp <- data.frame(taxa=i, taxa2=j, deg=1, mod_links=mod)
        df <- rbind(df, tmp)
      }
    }
  }
  out <- aggregate(list(mod_links=df$deg), by=list(taxa=df$taxa, module=df$mod), FUN=sum)
  return(out)
}

#compute within-module degree z-score which
#measures how well-connected a node is to other nodes in the module.

zscore <- function(mod.degree){
  ksi_bar <- aggregate(mod_links ~ module, data=mod.degree, FUN = mean)
  ksi_sigma <- aggregate(mod_links ~ module, data=mod.degree, FUN = sd)
  z <- NULL
  for(i in 1:dim(mod.degree)[1]){
    mod_mean <- ksi_bar$mod_links[which(ksi_bar$module == mod.degree$module[i])]
    mod_sig <- ksi_sigma$mod_links[which(ksi_bar$module == mod.degree$module[i])]
    z[i] <- (mod.degree$mod_links[i] - mod_mean)/mod_sig
  }
  z <- data.frame(row.names=rownames(mod.degree), z, module=mod.degree$module)
  return(z)
}

#The participation coefficient of a node measures how well a  node is distributed
# in the entire network. It is close to 1 if its links are uniformly
#distributed among all the modules and 0 if all its links are within its own module.

participation_coeffiecient <- function(mod.degree, total.degree){
  p <- NULL
  for(i in total.degree$taxa){
    ki <- subset(total.degree$total_links, total.degree$taxa==i)
    taxa.mod.degree <- subset(mod.degree$mod_links, mod.degree$taxa==i)
    p[i] <- 1 - (sum((taxa.mod.degree)**2)/ki**2)
  }
  p <- as.data.frame(p)
  return(p)
}


assign_module_roles <- function(zp){
  zp <- na.omit(zp)
  zp$roles <- rep(0, dim(zp)[1])
  outdf <- NULL
  for(i in 1:dim(zp)[1]){
    df <- zp[i, ]
    if(df$z < 2.5){ #non hubs
      if(df$p < 0.05){
        df$roles <- "ultra peripheral"
      }
      else if(df$p < 0.620){
        df$roles <- "peripheral"
      }
      else if(df$p < 0.80){
        df$roles <- "non hub connector"
      }
      else{
        df$roles <- "non hub kinless"
      }
    }
    else { # module hubs
      if(df$p < 0.3){
        df$roles <- "provincial hub"
      }
      else if(df$p < 0.75){
        df$roles <- "connector hub"
      }
      else {
        df$roles <- "kinless hub"
      }
    }
    if(is.null(outdf)){outdf <- df}else{outdf <- rbind(outdf, df)}
  }
  return(outdf)
}

plot_roles <- function(node.roles, roles.colors=NULL){
  x1<- c(0, 0.05, 0.62, 0.8, 0, 0.30, 0.75)
  x2<- c(0.05, 0.62, 0.80, 1,  0.30, 0.75, 1)
  y1<- c(-Inf,-Inf, -Inf, -Inf,  2.5, 2.5, 2.5)
  y2 <- c(2.5,2.5, 2.5, 2.5, Inf, Inf, Inf)
  lab <- c("ultra peripheral","peripheral" ,"non-hub connector","non-hub kinless","provincial"," hub connector","hub kinless")
  if(is.null(roles.colors)){roles.colors <- c("#E6E6FA", "#DCDCDC", "#F5FFFA", "#FAEBD7", "#EEE8AA", "#E0FFFF", "#F5F5DC")}
  p <- ggplot() + geom_rect(data=NULL, mapping=aes(xmin=x1, xmax=x2, ymin=y1,ymax=y2, fill=lab))
  p <- p + guides(fill=guide_legend(title="Topological roles"))
  p  <- p + scale_fill_manual(values = roles.colors)
  p <- p + geom_point(data=node.roles, aes(x=p, y=z,color=module)) + theme_bw()
  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("Participation Coefficient")+ylab(" Within-module connectivity z-score")
  return(p)
}

taxa.roles <- module.roles(comm_graph)
plot_roles(taxa.roles)


devtools::install_github('umerijaz/microbiomeSeq')
library(microbiomeSeq)

# explore the ecological processes underlying community assembly
## NST
library(NST)

### total community
#### read bacterial otu table
otu_table <- as.matrix(t(otu_table(tp_physeq)))

#### set the group
group <- data.frame(c(rep('all', nrow(otu_table))))
rownames(group)<-rownames(otu_table)

#### determine the NST of each taxa
tnst <- tNST(comm = otu_table, group = group, dist.method = "jaccard",
             abundance.weighted = TRUE, rand = 100, output.rand = T,
             nworker = 4, null.model= "PF", between.group = F,
             SES = T, RC = T)

nst.bt <- nst.boot(nst.result = tnst, group = NULL, rand = 99,
                   trace = TRUE, two.tail = FALSE, out.detail = T,
                   between.group = FALSE, nworker = 1)

NST.total <- data.frame(NST = (nst.bt$detail$NST.boot$all))

#### plot
library(ggplot2)
boxplot(NST.total*100,
        xlab = "NST (%)",
        #ylab = "Ozone",
        col = "#1b9e77",
        border = "black",
        horizontal = F,
        notch = TRUE
)

## carbon cycling community
#### read bacterial otu table
otu_table <- as.matrix(t(otu_table(tp_C_phylo)))

#### set the group
group <- data.frame(c(rep('all', nrow(otu_table))))
rownames(group)<-rownames(otu_table)

#### determine the NST of each taxa
tnst <- tNST(comm = otu_table, group = group, dist.method = "jaccard",
             abundance.weighted = TRUE, rand = 100, output.rand = T,
             nworker = 4 , null.model= "PF", between.group = F,
             SES = T, RC = T)

nst.bt <- nst.boot(nst.result = tnst, group = NULL, rand = 99,
                   trace = TRUE, two.tail = FALSE, out.detail = T,
                   between.group = FALSE, nworker = 1)

NST_carbon <- data.frame(NST = (nst.bt$detail$NST.boot$all))

#### plot
library(ggplot2)
boxplot(NST_carbon*100,
        xlab = "NST (%)",
        #ylab = "Ozone",
        col = '#d95f02',
        border = "black",
        horizontal = T,
        notch = TRUE 
)

nst.carbon

# NULL MODEL
null1 <- trans_nullmodel$new(meco_df, taxa_number = 1000, add_data = meco_df$sample_table)

# use pH as the test variable
null1$cal_mantel_corr(use_env = "MAP")
# return null1$res_mantel_corr
# plot the mantel correlogram
null1$plot_mantel_corr()

# null model run 500 times
null1$cal_ses_betampd(runs=99, abundance.weighted = TRUE)
# return t1$res_ses_betampd
null1$cal_ses_betampd





# NULL model analysis
### extract the most 1000 abundant ASVs
filter.abun <- function(phylo, N) {
  require(phyloseq)
  comun <- otu_table(phylo)
  rowmean <-sapply(1:nrow(comun), function(x) mean(comun[x, ]))
  comun <- comun[order(rowmean, decreasing = TRUE), ]
  taxa_list <- rownames(comun)[1:N]
  otu_tab <- comun[rownames(comun) %in% taxa_list, ]
  otu_tab <- data.frame(t(otu_tab))
  return(otu_tab)
}
### filter the tree with the most 1000 abundant ASVs
filter.tree <- function(phylo, otu_tab) {
  require(phyloseq)
  require(ape)
  tree <- phy_tree(phylo)
  tips<-colnames(otu_tab)
  tree<-keep.tip(tree, tips)
}

### extract the most 1000 abundant ASVs
commu <- filter.abun(tp_physeq, 1000)
ncol(commu)
sum(commu)/sum(otu_table(tp_physeq))
commu.tree <- filter.tree(tp_physeq, commu)

commu_carbon <- filter.abun(tp_C_phylo, 1000)
ncol(commu_carbon)
sum(commu_carbon)/sum(otu_table(tp_C_phylo))
commu.tree_carbon <- filter.tree(tp_C_phylo, commu_carbon)



### determinet the NRI and NTI
library(picante)
mpd.all<-mpd(commu, cophenetic(commu.tree), abundance.weighted = T)
ses.mpd.all <- ses.mpd(commu, cophenetic(commu.tree), 
                       null.model="taxa.labels", abundance.weighted = T, runs = 99)
ses.mntd.all <- ses.mntd(commu, cophenetic(commu.tree), 
                         null.model="taxa.labels", abundance.weighted = T, runs = 99)

env_phy_index <- data.frame(NRI = (ses.mpd.all$mpd.obs.z)*(-1),
                            NTI = (ses.mntd.all$mntd.obs.z)*(-1),
                            sample_data(tp_physeq))

### calculate the average LCBD for each site
phy_table_sel_agg <-  env_phy_index %>%
  dplyr::select(-c(3:7, 10:13)) %>%
  group_by(Site, Sitegroup1) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))

### determine the relationships between LCBD and envs using linear regression modes
phy_vars <- c("NRI", "NTI")
mode <- lapply(phy_vars, function(x) {
  lm(substitute(i ~ MAP, list(i = as.name(x))), data = phy_table_sel_agg)})
sum.mode <- lapply(mode, broom::glance)
### normality test using Shapiro-Wilk test 
res <- lapply(mode, residuals)
norm_test <- lapply(res, shapiro.test)
phy_norm_results <- data.frame(
  variables = phy_vars, 
  w = sapply(norm_test, "[[", "statistic"), 
  pvalue = sapply(norm_test, "[[", "p.value")
)
phy_norm_results

### extract the standardized regression coefficients
sd.coeff <- lapply(mode, QuantPsyc::lm.beta)
### arrange the table for plot
LCBD <- c(rep('LCBD', length(vars)))
sd.coeff <- sapply(sd.coeff, function(x){as.numeric(x[1])})
r.squared <- sapply(sum.mode, "[[", "r.squared")
adj.r.squared <- sapply(sum.mode, "[[", "adj.r.squared")
pvalue <- sapply(sum.mode, "[[", "p.value")
sig <- cut(pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
results_phy_index <- data.frame(phy_vars, sd.coeff, r.squared, adj.r.squared, pvalue, sig)
results_phy_index

### carbon cycling 
mpd.carbon<-mpd(commu_carbon, cophenetic(commu.tree_carbon), abundance.weighted = T)
ses.mpd.carbon <- ses.mpd(commu_carbon, cophenetic(commu.tree_carbon), 
                          null.model="taxa.labels", abundance.weighted=T, runs=99)
ses.mntd.carbon <- ses.mntd(commu_carbon, cophenetic(commu.tree_carbon), 
                            null.model="taxa.labels", abundance.weighted=T, runs=99)

env_phy_carbon <- data.frame(NRI = (ses.mpd.carbon$mpd.obs.z)*(-1),
                            NTI = (ses.mntd.carbon$mntd.obs.z)*(-1),
                            sample_data(tp_physeq))

### calculate the average NRI and NTI for each site
library(dplyr)
phy_table_sel_carbon_agg <-  env_phy_carbon %>% 
  dplyr::select(-c(3:7, 10:13)) %>%
  group_by(Site, Sitegroup1) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))
phy_table_sel_carbon_agg$SUVA254 <- log(phy_table_sel_carbon_agg$SUVA254)

### determine the relationships between LCBD and envs using linear regression modes
phy_vars <- c("NRI", "NTI")
mode <- lapply(phy_vars, function(x) {
  lm(substitute(i ~ MAP, list(i = as.name(x))), data = phy_table_sel_carbon_agg)})
sum.mode <- lapply(mode, broom::glance)
### normality test using Shapiro-Wilk test 
res <- lapply(mode, residuals)
norm_test <- lapply(res, shapiro.test)
phy_norm_carbon_results <- data.frame(
  variables = phy_vars, 
  w = sapply(norm_test, "[[", "statistic"), 
  pvalue = sapply(norm_test, "[[", "p.value")
)
phy_norm_carbon_results

### extract the standardized regression coefficients
sd.coeff <- lapply(mode, QuantPsyc::lm.beta)
### arrange the table for plot
LCBD <- c(rep('LCBD', length(vars)))
sd.coeff <- sapply(sd.coeff, function(x){as.numeric(x[1])})
r.squared <- sapply(sum.mode, "[[", "r.squared")
adj.r.squared <- sapply(sum.mode, "[[", "adj.r.squared")
pvalue <- sapply(sum.mode, "[[", "p.value")
sig <- cut(pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
results_phy_carbon_index <- data.frame(phy_vars, sd.coeff, r.squared, adj.r.squared, pvalue, sig)
results_phy_carbon_index

### PLOT
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

p_map_NRI <- ggplot(phy_table_sel_agg, aes(MAP, NRI)) + 
  geom_point(shape = 19, colour = "#1b9e77", size=3.5, alpha=0.8)+
  geom_smooth(method = "lm", size = 1.5, se = T, colour = 'black') +
  scale_y_continuous(limits = c(-6, 6)) +
  ylab('NRI')+xlab('MAP (mm)') +
  mytheme
p_map_NTI <- ggplot(phy_table_sel_agg, aes(MAP, NTI)) + 
  geom_point(shape = 19, colour = "#1b9e77", size=3.5, alpha=0.8)+
  geom_smooth(method = "lm", size = 1.5, se = T, colour = 'black') +
  scale_y_continuous(limits = c(-2, 5)) +
  ylab('NTI')+xlab('MAP (mm)') +
  mytheme
p_SUVA_NRI <- ggplot(phy_table_sel_carbon_agg, aes(SUVA254, NRI)) + 
  geom_point(shape = 19, colour = "#1b9e77", size=3.5, alpha=0.8)+
  geom_smooth(method = "lm", size = 1.5, se = T, colour = 'black') +
  scale_y_continuous(limits = c(-3, 5)) +
  ylab('NRI')+xlab('MAP (mm)') +
  mytheme
p_SUVA_NTI <- ggplot(phy_table_sel_carbon_agg, aes(SUVA254, NTI)) + 
  geom_point(shape = 19, colour = "#1b9e77", size=3.5, alpha=0.8)+
  geom_smooth(method = "lm", size = 1.5, se = T, colour = 'black') +
  scale_y_continuous(limits = c(-2, 4)) +
  ylab('NTI')+xlab('MAP (mm)') +
  mytheme
### Beta_NTI function
Beta_NTI<-function(tree, comun){
  require(picante)
  ## make sure the names on the phylogeny are ordered the same as the names in otu table
  match.phylo.comun = match.phylo.data(tree, t(comun));
  ## calculate empirical betaMNTD
  beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.comun$data), 
                                           cophenetic(match.phylo.comun$phy), abundance.weighted = T));
  identical(colnames(match.phylo.comun$data), colnames(beta.mntd.weighted)); # just a check, should be TRUE
  identical(colnames(match.phylo.comun$data), rownames(beta.mntd.weighted)); # just a check, should be TRUE
  # calculate randomized betaMNTD
  beta.reps = 999; # number of randomizations
  rand.weighted.bMNTD.comp = array(c(-999), dim = c(ncol(match.phylo.comun$data), ncol(match.phylo.comun$data), beta.reps));
  for (rep in 1:beta.reps) {
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.comun$data),
                                                          taxaShuffle(cophenetic(match.phylo.comun$phy)),
                                                          abundance.weighted = T, exclude.conspecifics = F));
    print(c(date(), rep));
  }
  weighted.bNTI = matrix(c(NA), nrow = ncol(match.phylo.comun$data), ncol = ncol(match.phylo.comun$data));
  dim(weighted.bNTI);
  for (columns in 1:(ncol(match.phylo.comun$data)-1)) {
    for (rows in (columns+1):ncol(match.phylo.comun$data)) {
      
      rand.vals = rand.weighted.bMNTD.comp[rows, columns,];
      weighted.bNTI[rows, columns] = (beta.mntd.weighted[rows, columns] - mean(rand.vals)) / sd(rand.vals);
      rm("rand.vals");
    };
  };
  rownames(weighted.bNTI) = colnames(match.phylo.comun$data);
  colnames(weighted.bNTI) = colnames(match.phylo.comun$data);
  results<-as.dist(weighted.bNTI);
  return(results)
}

### RC_bray function
raup_crick= function(comun, reps=999){
  require(ecodist) 
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(comun)
  gamma<-ncol(comun)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data = NA, nrow = n_sites, ncol = n_sites, dimnames = list(row.names(comun), row.names(comun)))
  ##make the comun matrix into a new, pres/abs. matrix:
  ceiling(comun/max(comun))->comun.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur <- apply(comun.inc, MARGIN = 2, FUN = sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance <- apply(comun, MARGIN = 2, FUN = sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in 1:(nrow(comun)-1)){
    for(null.two in (null.one+1):nrow(comun)){
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1 <- rep(0, gamma)
        com2 <- rep(0, gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(comun.inc[null.one,]), replace = FALSE, prob = occur)]<-1
        com1.samp.sp = sample(which(com1>0), (sum(comun[null.one,])-sum(com1)), replace = TRUE, prob = abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2], com1.samp.sp[,1], FUN = sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        ##same for com2:
        com2[sample(1:gamma, sum(comun.inc[null.two,]), replace = FALSE, prob = occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(comun[null.two,])-sum(com2)),replace = TRUE, prob = abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2], com2.samp.sp[,1], FUN = sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.comun = rbind(com1, com2); # null.comun;
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.comun, method = 'bray-curtis');
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(comun[c(null.one, null.two),], method = 'bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis == obs.bray);
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis < obs.bray);
      rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      ##modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
      rc = (rc-.5)*2
      results[null.two,null.one] = round(rc, digits=2); ##store the metric in the results matrix
      print(c(null.one, null.two,date()));
    }; ## end null.two loop
  }; ## end null.one loop
  
  results <- as.dist(results)
  return(results)
}

### calculate the BetaNTI
bNTI <- Beta_NTI(commu.tree, commu)
bNTI_1000 <- as.matrix(bNTI)
write.csv(bNTI_1000,"E:/thermokast_lakes/water_microbes/tibet_dada2_asv/results/tables/bNTI_1000.csv", quote = F)

bNTI_carbon <- Beta_NTI(commu.tree_carbon, commu_carbon)
bNTI_carbon_1000 <- as.matrix(bNTI_carbon)
write.csv(bNTI_carbon_1000,"E:/thermokast_lakes/water_microbes/tibet_dada2_asv/results/tables/bNTI_carbon_1000.csv", quote = F)

bNTI <- read.csv("E:/thermokast_lakes/water_microbes/tibet_dada2_asv/results/tables/bNTI_1000.csv", row.names = 1)
bNTI_carbon <- read.csv("E:/thermokast_lakes/water_microbes/tibet_dada2_asv/results/tables/bNTI_carbon_1000.csv", row.names = 1)
key_var <- env_phy_carbon %>% dplyr::select('Site', 'Sample_Name') %>% 
  group_by(Site) %>% dplyr::summarize(samples = list(Sample_Name))

mean_bNTI <- NULL
for (i in 1:nrow(key_var)){
  dist = usedist::dist_subset(bNTI, key_var$samples[[i]])
  tmp_bNTI = mean(dist)
  if (is.null(mean_bNTI)) {
    mean_bNTI = tmp_bNTI
  } else {
    mean_bNTI = c(mean_bNTI, tmp_bNTI)
  }
}
mean_bNTI_carbon <- NULL
for (i in 1:nrow(key_var)){
  dist = usedist::dist_subset(bNTI_carbon, key_var$samples[[i]])
  tmp_bNTI = mean(dist)
  if (is.null(mean_bNTI_carbon)) {
    mean_bNTI_carbon = tmp_bNTI
  } else {
    mean_bNTI_carbon = c(mean_bNTI_carbon, tmp_bNTI)
  }
}

null_data <- cbind(NRI_total = phy_table_sel_agg$NRI, NTI_totall = phy_table_sel_agg$NTI,
                   NRI_carbon = phy_table_sel_carbon_agg$NRI, NTI_carbon = phy_table_sel_carbon_agg$NTI,
                   mean_bNTI_total = mean_bNTI, mean_bNTI_carbon = mean_bNTI_carbon, phy_table_sel_agg[, -c(3,4)])


mode <- lm(mean_bNTI ~ phy_table_sel_agg$MAP)
summary(mode)
plot(mean_bNTI ~ phy_table_sel_agg$MAP)


mode <- lm(mean_bNTI_carbon ~ phy_table_sel_agg$SUVA254)
summary(mode)
plot(mean_bNTI_carbon ~ phy_table_sel_agg$SUVA254)

### calculate the RCbray
rc_bray <- raup_crick(commu)
rc_bray_1000<-as.matrix(rc_bray)
write.csv(rc_bray_1000,"E:/thermokast_lakes/water_microbes/tibet_dada2_asv/results/tables/rc_bray_1000.csv", quote = F);

### determine the relative contribution of each process in shaping community
null.index <- function(bNTI.df){
  bNTI.dist <- as.vector(bNTI.df)
  number <- c(1:length(bNTI.dist))
  null.value <- data.frame(number, bNTI = bNTI.dist)
  heterogeneous_selection <- nrow(subset(null.value, bNTI > 2))
  homogeneous_selection <- nrow(subset(null.value, bNTI < (-2)))
  stochastic_process <- nrow(subset(null.value,bNTI < 2 & bNTI > (-2)))
  null.pro <- data.frame(homogeneous_selection, heterogeneous_selection, stochastic_process)
  process <- colnames(null.pro)
  Proportion <- as.numeric(null.pro[1, ]/sum(null.pro[1, ]))
  df <- data.frame(process,Proportion)
  return(df)
}

#determine the relative contribution of each process in shaping community
null_process <- null.index(bNTI)
null_process_carbon <- null.index(bNTI_carbon)


#transform the wide dataframe to long format
null.3col <- function(bNTI.df){
  bnti_mat <- as.matrix(bNTI.df)
  bNTI.mat <- matrix(0,nrow = 147, ncol = 147)
  bNTI.mat[lower.tri(bNTI.mat)] <- bnti_mat[lower.tri(bnti_mat, diag=TRUE)]
  bNTI.dist <- as.dist(bNTI.df)
  library(NST)
  null.value<-dist.3col(bNTI.dist)
  colnames(null.value)<-c('name1','name2','bNTI')
  return(null.value)
}

bNTI <- null.3col(Beta_NTI)

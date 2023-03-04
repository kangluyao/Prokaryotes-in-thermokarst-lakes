otu <- data.frame(t(otu_table(tp_phylo)))
tree <- phy_tree(tp_phylo)

core.phy.rel <- core(water.rel, detection = 0.01/100, prevalence = 10/188, include.lowest = TRUE)


#MPD.WEIGHTED.ABUNDANCE
library(picante)
mpd.uw <- mpd(otu, cophenetic(tree), 
                   abundance.weighted = F)


#determine the MNTD, NTI, beta MNTD and beta NTI
NTI <- ses.mntd(otu, cophenetic(tree), null.model = "taxa.labels", 
                abundance.weighted= F, runs=999)
MNTD <- NTI$mntd.obs
NTI <- -(NTI$mntd.obs.z)

beta.mntd <- comdistnt(core, cophenetic(tree), abundance.weighted = F )

#
dat<-data.frame(SES.MNTD = ses.mntd.coretaxa$mntd.obs.z, env.table[,-(1:3)])
model <- lm(SES.MNTD ~ ., data=dat)

#Beta_NTI
Beta_NTI<-function(tree, otu_table){
  require(picante)
  ## make sure the names on the phylogeny are ordered the same as the names in otu table
  match.phylo.otu = match.phylo.data(tree, otu_table);
  ## calculate empirical betaMNTD
  beta.mntd.unweighted = as.matrix(comdistnt(t(match.phylo.otu$data),
                                           cophenetic(match.phylo.otu$phy),
                                           abundance.weighted=F));
  identical(colnames(match.phylo.otu$data),colnames(beta.mntd.unweighted)); # just a check, should be TRUE
  identical(colnames(match.phylo.otu$data),rownames(beta.mntd.unweighted)); # just a check, should be TRUE
  # calculate randomized betaMNTD
  beta.reps = 999; # number of randomizations
  rand.unweighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),
                                                 ncol(match.phylo.otu$data),
                                                 beta.reps));
  for (rep in 1:beta.reps) {
    rand.unweighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),
                                                          taxaShuffle(cophenetic(match.phylo.otu$phy)),
                                                          abundance.weighted=F,exclude.conspecifics = F));
    print(c(date(),rep));
  }
  unweighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
  dim(unweighted.bNTI);
  for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
    for (rows in (columns+1):ncol(match.phylo.otu$data)) {
      
      rand.vals = rand.unweighted.bMNTD.comp[rows,columns,];
      unweighted.bNTI[rows,columns] = (beta.mntd.unweighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
      rm("rand.vals");
    };
  };
  rownames(unweighted.bNTI) = colnames(match.phylo.otu$data);
  colnames(unweighted.bNTI) = colnames(match.phylo.otu$data);
  results<-as.dist(unweighted.bNTI);
  return(results)
}


beta_nti <- Beta_NTI(tree, otu)
beta_nti1 <-as.matrix(beta_nti)
#write.csv(beta_nti1, './tables/beta_nti.csv', quote=F)

#RC_bray
raup_crick= function(spXsite, reps=999){
  require(ecodist)
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite/max(spXsite))->spXsite.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(spXsite, MARGIN=2, FUN=sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in 1:(nrow(spXsite)-1)){
    for(null.two in (null.one+1):nrow(spXsite)){
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        ##same for com2:
        com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.spXsite = rbind(com1,com2); # null.spXsite;
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.spXsite,method='bray-curtis');
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(spXsite[c(null.one,null.two),],method='bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      ##modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
      rc = (rc-.5)*2
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      print(c(null.one,null.two,date()));
    }; ## end null.two loop
  }; ## end null.one loop
  
  results<-as.dist(results)
  return(results)
}
total.otu <- as(otu_table(water_physeq), 'matrix')
total.otu.df <- t(total.otu)
rcbray <- raup_crick(total.otu.df)

null.index<-function(bNTI.df,rcbray.df){
  bnti_mat <- as.matrix(bNTI.df)
  bNTI.mat<-matrix(0,nrow = 147,ncol = 147)
  bNTI.mat[lower.tri(bNTI.mat)] <- bnti_mat[lower.tri(bnti_mat, diag=TRUE)]
  bNTI.dist<-as.dist(bNTI.mat)
  #
  rcbray.mat <- as.matrix(rcbray.df)
  RC.mat<-matrix(0,nrow = 147,ncol = 147)
  RC.mat[lower.tri(RC.mat)] <- rcbray.mat [lower.tri(rcbray.mat, diag=TRUE)]
  RC.dist<-as.dist(RC.mat)
  library(NST)
  null.value<-cbind(dist.3col(bNTI.dist),dist.3col(RC.dist)[3])
  colnames(null.value)<-c('name1','name2','bNTI','RC_bray')
  homogeneous_selection<-nrow(subset(null.value,bNTI>2))
  heterogeneous_selection<-nrow(subset(null.value,bNTI<(-2)))
  dispersal_limitation<-nrow(subset(null.value,bNTI<2 & bNTI>(-2) & RC_bray>0.95))
  homogenizing_dispersal<-nrow(subset(null.value,bNTI<2 & bNTI>(-2) & RC_bray<(-0.95)))
  undominated<-nrow(subset(null.value,bNTI<2 & bNTI>(-2) & RC_bray<0.95 & RC_bray>(-0.95)))
  null.pro<-data.frame(homogeneous_selection,heterogeneous_selection,dispersal_limitation,
                       homogenizing_dispersal,undominated)
  process<-colnames(null.pro)
  Proportion<-as.numeric(null.pro[1,]/sum(null.pro[1,]))
  df<-data.frame(process,Proportion)
  return(df)
}


null.fraction <- null.index(beta_nti,rcbray_core)
null.fraction$process <- c('Homogeneous selection', 'Heterogeneous selection',
                           'Dispersal limitation', 'Homogenizing dispersal', 
                           'Undominated')
library(ggplot2)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = gg_color_hue(n)
p.null<-ggplot(null.fraction,aes(x = process, y = Proportion*100, fill = process))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values= rev(cols))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 80))+
  labs(y = 'Proportion (%)', x = 'Process', fill = 'Process')+
  theme_bw()+
  theme(axis.title =element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(), legend.position = c(0.2, 0.75))
p.null


##

# create environmental distance matrices
DOC.dist = vegdist(scale(env.table$DOC), 'euclidean')
TC.dist = vegdist(scale(env.table$TC), 'euclidean')
TN.dist = vegdist(scale(env.table$TN), 'euclidean')
NH4_N.dist = vegdist(scale(env.table$NH4_N), 'euclidean')
NO3_N.dist = vegdist(scale(env.table$NO3_N), 'euclidean')
DIN.dist = vegdist(scale(env.table$DIN), 'euclidean')
DON.dist = vegdist(scale(env.table$DON), 'euclidean')
S275_295.dist = vegdist(scale(env.table$S275_295), 'euclidean')
SUVA254.dist = vegdist(scale(env.table$SUVA254), 'euclidean')
bix.dist = vegdist(scale(env.table$bix), 'euclidean')
fi.dist = vegdist(scale(env.table$fi), 'euclidean')
hix.dist = vegdist(scale(env.table$hix), 'euclidean')
a300.dist = vegdist(scale(env.table$a300), 'euclidean')
pH.dist = vegdist(scale(env.table$pH), 'euclidean')
DO.dist = vegdist(scale(env.table$DO), 'euclidean')
Conductivity.dist = vegdist(scale(env.table$Conductivity), 'euclidean')
Salinity.dist = vegdist(scale(env.table$Salinity), 'euclidean')
Temp.dist = vegdist(scale(env.table$Temp), 'euclidean')
Depth.dist = vegdist(scale(env.table$Depth), 'euclidean')
Ca.dist = vegdist(scale(env.table$Ca), 'euclidean')
K.dist = vegdist(scale(env.table$K), 'euclidean')
Mg.dist = vegdist(scale(env.table$Mg), 'euclidean')
Na.dist = vegdist(scale(env.table$Na), 'euclidean')
MAT.dist = vegdist(scale(env.table$MAT), 'euclidean')
MAP.dist = vegdist(scale(env.table$MAP), 'euclidean')

env.table <- data.frame(env.table)
dom_vars <- c('S275_295', 'SUVA254', 'TN', 'TC')
climate_vars <- c('MAP', 'MAT')
physicochemical_vars <- c('pH', 'Temp', 'Depth', 'DO', 'K', 'Conductivity')

pH_vars <- c('pH', 'Temp', 'Depth', 'DO')
conductivity_vars <- c('K', 'Ca', 'Na', 'Mg', 'Conductivity')

dom_distmat <- vegdist(scale(env.table[ ,dom_vars]), 'euclidean')
climate_vars_distmat <- vegdist(scale(env.table[ ,climate_vars]), 'euclidean')
physicochemical_distmat <- vegdist(scale(env.table[ ,physicochemical_vars]), 'euclidean')

conductivity_vars_distmat <- vegdist(scale(env.table[ ,conductivity_vars]), 'euclidean')
pH_vars_distmat <- vegdist(scale(env.table[ ,pH_vars]), 'euclidean')

# plot
all.distmat <- as.data.frame(cbind(beta_nti, beta.mntd.core, bray.distmat, 
                                   unifac.dismat, fun.dismat, climate_vars_distmat, 
                                   dom_distmat, physicochemical_distmat, pH_vars_distmat,
                                   conductivity_vars_distmat, DOC.dist, TC.dist, TN.dist,
                                   NH4_N.dist, NO3_N.dist, DIN.dist, DON.dist, S275_295.dist,
                                   SUVA254.dist, bix.dist, fi.dist, a300.dist, pH.dist,
                                   DO.dist, Conductivity.dist, Salinity.dist, Temp.dist,
                                   Depth.dist, Ca.dist, K.dist, Mg.dist, Na.dist, MAT.dist,
                                   MAP.dist))


varpart(beta_nti, env.table[ ,climate_vars], env.table[ ,dom_vars], env.table[ ,physicochemical_vars])



library(ggplot2)
library(ggpubr)
βNTI.plot <- ggplot(all.distmat, aes(x = S275_295.dist, y = beta_nti))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="βNTI",x="ΔS275_295")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(colour = 'black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
βNTI.plot

β.mntd.plot <- ggplot(all.distmat, aes(x =  S275_295.dist, y = beta.mntd.core))+
  geom_point(size = 3, alpha=0.1,colour = "DarkGreen")+
  #scale_color_manual(values=c("Tan4","OliveDrab3","DarkGreen"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc = 0.80,label.y.npc = 0.02)+
  theme_bw()+
  labs(y="βMNTD",x="ΔS275_295")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(colour = 'black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
β.mntd.plot


#fig3
fig3 <- plot_grid(p.null, β.mntd.plot, labels = 'auto', ncol = 2, 
                  label_x = .01, label_y = 0.99, 
                  hjust = 0, label_size=14,align = "v")


#据AIC逐步回归
library(MASS)
step <- stepAIC(model, direction="both",trace=F)
step$call
best.model <- lm(formula = SES.MNTD ~ DOC + TOC + TN + SUVA254 + a300 + bix + 
                   fi + hix + Conductivity + Salinity + Depth, data = dat)
library(relaimpo)
impor<-calc.relimp(best.model)
impor.ses.matd<-(impor$lmg)*100
Hemicellulose<-data.frame(variable=names(impor.Hemicellulose),importance=impor.Hemicellulose,stringsAsFactors = FALSE)
Hemicellulose$variable<-reorder(Hemicellulose$variable,Hemicellulose$importance)


library(randomForest)
library(mlbench)
library(caret)
##Extend Caret
customRF <- list(type = "Regression", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

# train model
control <- trainControl(method="repeatedcv", number=10, repeats=10)
tunegrid <- expand.grid(.mtry=c(1:10), .ntree=c(200, 400, 600, 800, 1000))
custom_bacteria <- train(SES.MNTD ~ .,
                         data=dat, method=customRF, metric='RMSE', 
                         tuneGrid=tunegrid, trControl=control, na.action=na.exclude)
summary(custom_bacteria)
plot(custom_bacteria)
d1<-varImp(custom_bacteria, scale = T)$importance
d1<-cbind(variable=rownames(d1),d1)
d1$variable<-reorder(d1$variable,d1$Overall)
library(ggplot2)
theme_set(theme_classic())
pd1<-ggplot(d1,aes(x=variable,y= Overall))+
  geom_bar(stat="identity",fill="steelblue",colour=NA)+
  scale_y_continuous(expand = c(0,0))+
  xlab('importance')+
  theme(axis.title.y =element_blank(),
        axis.text.y = element_text())+
  coord_flip()
pd1

Predicted.abun<-predict(custom_bacteria,data=diversity)
RF_bacteria_data<-data.frame(Actual_abun=diversity$Abundance,Predicted_abun=Predicted.abun)

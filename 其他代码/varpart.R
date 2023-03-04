#check collinearities among all variables， remove the varibles with vif > 10
ord <- capscale(tp_otu_hel_dist ~., env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_C_otu_hel_dist ~ DOC + TN + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Salinity + Ca + Mg + K + Na + comp1 + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX + HIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_C_otu_hel_dist ~ DOC + TN + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + Na + comp1 + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX + HIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_C_otu_hel_dist ~ DOC + TN + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + Na + comp1 + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_C_otu_hel_dist ~ DOC + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + Na + comp1 + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_C_otu_hel_dist ~ DOC + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + comp1 + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_C_otu_hel_dist ~ DOC + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_C_otu_hel_dist ~ DOC + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + 
                  comp2 + comp3 + comp4 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_C_otu_hel_dist ~ DOC + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Ca + Mg + K + comp2 + comp3 + comp4 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
names(ord.vif)
sel.env <- env_df[,names(ord.vif)]
##forward selection of environmental variables
mod1 <- capscale(tp_C_otu_hel_dist ~. , sel.env, add=T)
mod0 <- capscale(tp_C_otu_hel_dist ~1 , sel.env, add=T)
mod <- ordiR2step(mod0, scope = formula(mod1), perm.max = 999, trace = F)
Env.P <- anova.cca(mod,permutations = 999)[[4]][1]
Env.P
Env_sel <- mod$CCA$biplot
Env_C_sel <- sel.env[ , rownames(Env_sel)]

climate_vars <- Env_sel[ , c('MAP', 'MAT')]
DOM_vars <- Env_sel[ , c('SUVA254', 'S275_295', 'BIX', 'comp4', 'a300')]
phyche_vars <- Env_sel[ , c('pH', 'Depth', 'K', 'Ca', 'DON', 'NO3_N')]

mod_C <- varpart(tp_C_otu_hel_dist, ~ MAP + MAT, ~ SUVA254 + S275_295 + BIX + comp4 + a300,
                 ~pH + Depth + K + Ca + DON + NO3_N + NH4_N, data = Env_sel)
mod_C


#extract the prokaryotic community involve in nitrogen cycling on the TP
tp_nitrogen_phy <- subset_taxa(tp_phylo, OTU %in% meta_n_cycle_otus)
tp_nitrogen_phy <- prune_samples(sample_sums(tp_nitrogen_phy) > 0, tp_nitrogen_phy)
tp_N_otu <- as.matrix(t(otu_table(tp_nitrogen_phy)))
tp_N_otu_hel <- decostand(tp_N_otu, 'hellinger')
tp_N_otu_hel_dist<-vegdist(tp_N_otu_hel, 'bray',upper=F)
#extract the env table
env_tp <- sample_data(tp_nitrogen_phy)
env_df <- env_tp[,-c(1:14)]
env_df <- data.frame(decostand(env_df,"standardize"))
#check collinearities among all variables， remove the varibles with vif > 10
ord <- capscale(tp_N_otu_hel_dist ~., env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_N_otu_hel_dist ~ DOC + TN + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Salinity + Ca + Mg + K + Na + comp1 + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX + HIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_N_otu_hel_dist ~ DOC + TN + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + Na + comp1 + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX + HIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_N_otu_hel_dist ~ DOC + TN + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + Na + comp1 + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_N_otu_hel_dist ~ DOC + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + Na + comp1 + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_N_otu_hel_dist ~ DOC + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + comp1 + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_N_otu_hel_dist ~ DOC + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + 
                  comp2 + comp3 + comp4 + comp5 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_N_otu_hel_dist ~ DOC + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Conductivity + Ca + Mg + K + 
                  comp2 + comp3 + comp4 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
ord <- capscale(tp_N_otu_hel_dist ~ DOC + DON + NH4_N + NO3_N + S275_295 + SUVA254 + a300 + MAP +
                  MAT + Depth + DO + pH + Ca + Mg + K + comp2 + comp3 + comp4 + FluI + BIX, env_df)
ord.vif <- vif.cca(ord)
ord.vif
names(ord.vif)
sel.env <- env_df[,names(ord.vif)]
##forward selection of environmental variables
mod1 <- dbrda(tp_N_otu_hel_dist ~. , sel.env, add=T)
mod0 <- dbrda(tp_N_otu_hel_dist ~1 , sel.env, add=T)
mod <- ordiR2step(mod0, scope = formula(mod1), perm.max = 999, trace = F)
Env.P <- anova.cca(mod,permutations = 999)[[4]][1]
Env.P
Env_sel <- mod$CCA$biplot
Env_N_sel <- sel.env[ , rownames(Env_sel)]
climate_vars <- Env_sel[ , c('MAP', 'MAT')]
DOM_vars <- Env_sel[ , c('SUVA254', 'S275_295', 'BIX', 'comp4', 'a300')]
phyche_vars <- Env_sel[ , c('pH', 'Depth', 'K', 'Ca', 'DON', 'NO3_N')]
mod_N <- varpart(tp_N_otu_hel_dist, ~ MAP + MAT, ~ DON + NH4_N,
                 ~ pH + K + Mg + S275_295 + a300, data = Env_sel)
mod_N

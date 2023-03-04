# beta diversity analysis
## local contribution to beta diversity (LCBD) analysis for total community
library(adespatial)
beta_tax_div <- beta.div(t(as.matrix(otu_table(tp_physeq))), 
                         method = "hellinger", sqrt.D = FALSE, samp = TRUE, 
                         nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)
env_div <- data.frame(LCBD = beta_tax_div$LCBD, sample_data(tp_physeq))
### calculate the average LCBD for each site
library(dplyr)
env_div_agg <-  env_div %>% 
  dplyr::select(-c(2:6, 9, 12)) %>%
  group_by(Site, Sitegroup1) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))
#write.csv(env_div_agg, file = './tibet_dada2_asv/results/tables/env_div_agg.csv')

### determine the relationships between LCBD and envs using linear regression modes
vars <- c("MAT", "MAP", "DOC", "S275_295", "SUVA254", "a300", "BIX", "HIX",
          "TN", "NH4_N", "DO", "pH", "Conductivity", "Salinity",
          "K", "Ca", "Na", "Mg")
mode <- lapply(vars, function(x) {
  lm(substitute(LCBD ~ i, list(i = as.name(x))), data = env_div_agg)})
sum.mode <- lapply(mode, broom::glance)
### normality test using Shapiro-Wilk test 
res <- lapply(mode, residuals)
norm_test <- lapply(res, shapiro.test)
norm_results <- data.frame(
  variables = vars, 
  w = sapply(norm_test, "[[", "statistic"), 
  pvalue = sapply(norm_test, "[[", "p.value")
)
norm_results

### extract the standardized regression coefficients
sd.coeff <- lapply(mode, QuantPsyc::lm.beta)
### arrange the table for plot
LCBD <- c(rep('LCBD', length(vars)))
sd.coeff <- sapply(sd.coeff, function(x){as.numeric(x[1])})
r.squared <- sapply(sum.mode, "[[", "r.squared")
adj.r.squared <- sapply(sum.mode, "[[", "adj.r.squared")
pvalue <- sapply(sum.mode, "[[", "p.value")
sig <- cut(pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
results <- data.frame(vars, LCBD, sd.coeff, r.squared, adj.r.squared, pvalue, sig)
results

#model selection
library(MASS)
library(glmulti)
A1 <- glmulti(LCBD ~ MAP + MAT + S275_295 + SUVA254 + a300 + BIX + HIX +
                TN + Conductivity + Salinity + Mg + K + Na, data=env_div_agg,
              level=1, fitfunction=lm, crit="aicc", confsetsize= 2^13, plotty = F, trace = 0)
top <- weightable(A1)
###  models with values more than 2 units away are considered substantially 
### less plausible than those with AICc values closer to that of the best model. 
### refrence:Anderson, D. R. (2007). Model based inference in the life sciences: A primer on evidence. New York: Springer. 
top_1 <- top[top$aicc <= min(top$aicc) + 2,] # 
top_1

modes_inf <- NULL
for(i in 1:nrow(top_1)){
  rse_sum <- summary(A1@objects[[i]])
  adj.r.squared <- rse_sum$adj.r.squared # obtain the adjust r squared
  multicollinearity <- any(car::vif(A1@objects[[i]]) > 2) # check the multicollinearity
  tmp <- data.frame(adj.r.squared, multicollinearity)
  if(is.null(modes_inf)){
    modes_inf<-tmp
  } else {
    modes_inf <- rbind(modes_inf,tmp)
  } 
}
modes_inf <- cbind(top_1, modes_inf)
modes_inf

vpa.mod <- varpart(env_div_agg$LCBD, ~ env_div_agg$HIX,
                   ~ env_div_agg$MAP)
plot(vpa.mod)


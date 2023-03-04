#ComBat-Seq for batch adjustment on DNASeq count data
env_mat <- metadata[, c('Region', 'Site', 'batch', 'Sitegroup')]
batch = as.factor(env_mat$batch)
modcombat = model.matrix(~1, data = env_mat)
combat_otu_table = ComBat(dat = meta.otu.table, batch=batch, mod=modcombat,
                          par.prior=F, prior.plots=T)

# Now turn into a DGEList
y = DGEList(counts=meta.otu.table, remove.zeros = TRUE)
# Calculate the normalization factors
z = calcNormFactors(y, method='RLE')
# Check for division by zero inside `calcNormFactors`
if( !all(is.finite(z$samples$norm.factors)) ){
  stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
}
# Estimate dispersions
return(estimateTagwiseDisp(estimateCommonDisp(z)))


dds <- estimateSizeFactors(meta.otu.table)
pheno = pData(bladderEset)
edata = exprs(bladderEset)
mod = model.matrix(~as.factor(Region), data=env_mat)
mod0 = model.matrix(~1,data=env_mat)
n.sv = num.sv(meta.otu.table,mod,method="leek")
n.sv
svobj = sva(meta.otu.table,mod,mod0,n.sv=n.sv)
pValues = f.pvalue(meta.otu.table,mod,mod0)
qValues = p.adjust(pValues,method="BH")

meta.otu.table[1:5,1:5]
meta.otu.table1 <- as.matrix(meta.otu.table + 1)
batch <- metadata$batch
group <- metadata$Region1
cov1 <- metadata$Sitegroup1
cov2 <- metadata$Site1
covar_mat <- as.matrix(cov1, cov2)
adjusted_otu_table <- ComBat_seq(meta.otu.table, batch = batch)
adjusted_otu_table[1:5,1:5]
library(staRdom)
system.file()
system.file("extdata/EEMs", package = "staRdom")
folder <- system.file("extdata/EEMs/", package = "staRdom")  # folder containing example EEMs
library("eemR")
eem_list <- eem_read(folder, recursive = TRUE, import_function = eem_csv)  # in case you use your own
# data, just replace folder by a path. e.g.
# 'C:/folder/another folder' and change
# import_function according to instrument. eem_list
# <- eem_read(folder, import_function = 'cary')
eem_overview_plot(eem_list, spp = 9, contour = TRUE)
absorbance_path = system.file("extdata/absorbance", 
                              package = "staRdom")  # load example data, set a
# path without using system.file to use your own
# data e.g. 'C:/folder/another folder'
absorbance <- absorbance_read(absorbance_path, cores = cores)  # load csv or txt tables in folder
metatable <- system.file("extdata/metatable_dreem.csv", 
                         package = "staRdom")
meta <- read.table(metatable, header = TRUE, sep = ",", 
                   dec = ".", row.names = 1)
problem <- eem_checkdata(eem_list, absorbance, meta, 
                         metacolumns = c("dilution"), error = FALSE)
eem_list
absorbance_path = system.file("extdata/absorbance", 
                              package = "staRdom")  # load example data, set a path without using system.file to use your own data e.g. 'C:/folder/another folder'
summary(eem_list_test)


# import EEMs data
library(staRdom)
library(eemR)
cores <- detectCores(logical = FALSE)

eem_list <- eem_read("E:/thermokast_lakes/pond_data/EEMs/EEMs_t0/csvfile", 
                          recursive = TRUE, import_function = eem_csv)
#eem_overview_plot(eem_list_test, spp = 9, contour = TRUE)
blank <- eem_read("E:/thermokast_lakes/pond_data/EEMs/EEMs_t0/blank", 
                  recursive = TRUE, import_function = eem_csv)
absorbance <- read.csv('E:/thermokast_lakes/pond_data/EEMs/absorbance/absorbance_t0.csv',
                       header = T, stringsAsFactors = F)
meta <- read.csv('E:/thermokast_lakes/pond_data/EEMs/EEMs_t0/metatable_eem.csv', header = TRUE, row.names = 1)
problem <- eem_checkdata(eem_list,absorbance, meta, metacolumns = c("dilution"),error=FALSE)
# data correction
library(dplyr)
eem_list <- eem_remove_blank(eem_list, blank)
eem_list <- eem_remove_scattering(eem_list, "rayleigh", 1, 10) %>% 
  eem_remove_scattering("raman", 1, 10)
eem_list <- eem_ife_correction(eem_list, absorbance, cuvl = 1 )
eem_list <- eem_raman_normalisation(eem_list, blank)

#Remove blanks from sample set
eem_list <- eem_extract(eem_list, c("nano", "miliq", "milliq", "mq", "blank"),ignore_case = TRUE)
#Remove and interpolate scattering
remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
remove_scatter_width <- c(15,15,15,15)
eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter, remove_scatter_width =
                           remove_scatter_width)
eem_overview_plot(eem_list, spp=9, contour = TRUE)


eem_list <- eem_interp(eem_list, cores = cores, type = 1, extend = FALSE)
eem_overview_plot(eem_list, spp=9, contour = TRUE)
#Correct for dilution
dil_data <- meta["dilution"]
eem_list <- eem_dilution(eem_list,dil_data)
#Smooth data
eem4peaks <- eem_smooth(eem_list, n = 4, cores = cores)
summary(eem_list)


# indexes
bix <- eem_biological_index(eem4peaks, verbose = TRUE)
fi <- eem_fluorescence_index(eem4peaks, verbose = TRUE)
hix <- eem_humification_index(eem4peaks, scale = TRUE)
coble_peaks <- eem_coble_peaks(eem4peaks)

indices_peaks <- bix %>% full_join(coble_peaks, by = "sample") %>% 
  full_join(fi, by = "sample") %>% full_join(hix, by = "sample")
indices_peaks

write.csv(indices_peaks,file = 'E:/thermokast_lakes/pond_data/EEMs/EEMs_index_t0.csv')

#=============Creating a PARAFAC model==========
# minimum and maximum of numbers of components
dim_min <- 3
dim_max <- 7
nstart <- 50 # number of similar models from which best is chosen
maxit = 5000 # maximum number of iterations in PARAFAC analysis
ctol <- 10^-6 # tolerance in PARAFAC analysis
# calculating PARAFAC models, one for each number of components
pf1 <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE, 
                   const = c("uncons","uncons", "uncons"), maxit = maxit, 
                   nstart = nstart, ctol = ctol, cores = cores)

# same model but using non-negative constraints
pf1n <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), 
                    normalise = FALSE, const = c("nonneg","nonneg", "nonneg"), 
                    maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
# rescale B and C modes to a maximum fluorescence of 1 for each component
pf1 <- lapply(pf1, eempf_rescaleBC, newscale = "Fmax")
pf1n <- lapply(pf1n, eempf_rescaleBC, newscale = "Fmax")

# This plot is not shown, because the components violate the assumptions for fluorescence peaks (negative fluorescence). Please try, if you are interested.
eempf_compare(pf1, contour = TRUE)
eempf_compare(pf1n, contour = TRUE)

# check for correlation between components table
eempf_cortable(pf1n[[4]], normalisation = FALSE)
eempf_corplot(pf1n[[4]], progress = FALSE, normalisation = FALSE)

pf2 <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), 
                   normalise = TRUE, const = c("nonneg","nonneg", "nonneg"),
                   maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
# rescale B and C modes
pf2 <- lapply(pf2, eempf_rescaleBC, newscale = "Fmax")

eempf_plot_comps(pf2, contour = TRUE, type = 1)
# calculate leverage
cpl <- eempf_leverage(pf2[[4]])
# plot leverage (nice plot)
eempf_leverage_plot(cpl,qlabel=0.1)
# plot leverage, not so nice plot but interactive to select what to exclude
# saved in exclude, can be used to start over again with eem_list_ex <- eem_list %>%
eem_exclude(exclude) above
exclude <- eempf_leverage_ident(cpl,qlabel=0.1)
# samples, excitation and emission wavelengths to exclude, makes sense after calculation of leverage
exclude <- list("ex" = c(),
                "em" = c(),
                "sample" = c("Z5_H1","Z37_H4")
)

#A new PARAFAC model is then generated without outliers:
# exclude outliers if neccessary. if so, restart analysis
eem_list_ex <- eem_exclude(eem_list, exclude)
pf3 <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, maxit = maxit, nstart =
                     nstart, ctol = ctol, cores = cores)
pf3 <- lapply(pf3, eempf_rescaleBC, newscale = "Fmax")
eempf_plot_comps(pf3, contour = TRUE, type = 1)
eempf_leverage_plot(eempf_leverage(pf3[[4]]),qlabel=0.1)
#Recalculating the model with increased accuracy
ctol <- 10^-8 # decrease tolerance in PARAFAC analysis
nstart = 50 # increase number of random starts
maxit = 10000 # increase number of maximum interations
pf4 <- eem_parafac(eem_list_ex, comps = 6, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"),
                   maxit = maxit, nstart = nstart, ctol = ctol, output = "all", cores = cores)
pf4 <- lapply(pf4, eempf_rescaleBC, newscale = "Fmax")
eempf_convergence(pf4[[1]])
eempf_compare(pf4, contour = TRUE)
eempf_leverage_plot(eempf_leverage(pf4[[1]]))
eempf_corplot(pf4[[1]], progress = FALSE)


file.edit(system.file("EEM_simple_analysis.Rmd", package = "staRdom"))

fix(absorbance)

devtools::install_github("tidyverse/dplyr")
install.packages("devtools")
devtools::install_github("tidyverse/dplyr")
.libPaths("C:/Users/kangluyao/Documents/R/win-library/3.6")
install.packages("data.table")
devtools::install_github("Rdatatable/data.table")

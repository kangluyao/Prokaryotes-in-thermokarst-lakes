remotes::install_github(
  repo = "blasbenito/spatialRF", 
  ref = "main",
  force = TRUE,
  quiet = TRUE
)
#loading packages
library(spatialRF)
library(kableExtra)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(randomForestExplainer)
library(pdp)

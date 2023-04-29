# Prokaryotes in thermokarst lakes
This repository contains R codes and data for microbial community analysis of the manuscript 'Climatic factors dominate prokaryotic communities in the water of thermokarst lakes across the northern Hemispherical permafrost regions'

## Overview

------------------------------------------------------------------------

> Project repository associated with the following manuscripts

------------------------------------------------------------------------

* Kang et al., (2023).[Climatic factors dominate prokaryotic communities in the water of thermokarst lakes across the northern Hemispherical permafrost regions](https://psyarxiv.com/w8vhk/). (*in revision*). 

### Structure

```
.
├── script                                  <-- R codes for the statistical analysis
    ├── qiime2
    ├── 2_composition_analysis.R            <-- codes for taxonomic compositional analysis
    ├── 3_Comparison between PA & TP.R      <-- codes for comparison analysis between PA and TP
    ├── 5_cor_DOM_taxa_fun.R                <-- test the relationship between the dominant taxa and DOM properties
    ├── 6_LCBD_analysis.R                   <-- codes for LCBD analysis
    └── FARPROTAX_Carbon.R                  <-- codes for functional groups analysis
└── data1                     
    ├── qiime2                              <-- files and code to reproduce CogSci paper
        ├── thermokarst-lakes-rep-seqs.qza  <-- represent sequences output file from qiime2
        ├── thermokarst-lakes-table.qza     <-- otutable output from qiime2
    └── meta_data              
        ├── meta_otu_table.txt              <-- otu table file
        ├── taxonomy.txt                    <-- taxonomimc table file
        └── sample_data.txt                 <-- environmental table file

```

# Prokaryotes in thermokarst lakes

> This project repository associated with the following manuscripts:

* Kang *et al*., (2023).Climatic factors dominate prokaryotic communities in the water of thermokarst lakes across the northern Hemispherical permafrost regions. (**in revision**). 

## Structure

```
.
├── data1                                   <-- relevant data necessary for reproducing results 
    ├── qiime2                              <-- output files from the qiime2
        ├── thermokarst-lakes-rep-seqs.qza  <-- represent sequences output file from qiime2
        ├── thermokarst-lakes-table.qza     <-- otutable output from qiime2
        └── meta-analysis-DADA2.sh          <-- shell scripts for the bioinformatic analysis
    └── meta_data              
        ├── meta_otu_table.txt              <-- otu table file
        ├── taxonomy.txt                    <-- taxonomimc table file
        └── sample_data.txt                 <-- environmental table file
└── script                                  <-- R codes for the statistical analysis
    ├── composition_analysis.R              <-- codes for taxonomic compositional analysis
    ├── comparison between PA & TP.R        <-- codes for comparison analysis between pan-Arctic and Tibetan Plateau
    ├── cor_DOM_taxa_fun.R                  <-- test the relationship between the dominant taxa and DOM properties
    ├── LCBD_analysis.R                     <-- codes for LCBD analysis
    └── FARPROTAX.R                         <-- codes for functional groups analysis

```

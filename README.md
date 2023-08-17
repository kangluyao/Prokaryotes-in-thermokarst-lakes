# Prokaryotes in thermokarst lakes

> This project repository associated with the following manuscripts:

* Kang Luyao, Chen Leiyi, Li Ziliang, Wang Jianjun, Xue Kai, Deng Ye, Delgado-Baquerizo Manuel, Yang Yuanhe (2023). Patterns and drivers of prokaryotic communities in thermokarst lake water across permafrost regions in the Northern Hemisphere. (**in revision**). 

## Structure

```
.
├── data1
    ├── qiime2                              <-- output files from the qiime2
        ├── thermokarst-lakes-rep-seqs.qza  <-- represent sequences output file from qiime2
        └── thermokarst-lakes-table.qza     <-- otutable output from qiime2
    └── meta_data                           <-- relevant data necessary for reproducing results
        ├── meta_otu_table.txt              <-- otu table file
        ├── taxonomy.txt                    <-- taxonomimc table file
        └── sample_data.txt                 <-- environmental table file
└── script
    ├── composition_analysis.R              <-- codes for taxonomic compositional analysis
    ├── comparison between PA & TP.R        <-- codes for comparison analysis between pan-Arctic and Tibetan Plateau
    ├── cor_DOM_taxa_fun.R                  <-- test the relationship between the dominant taxa and DOM properties
    ├── LCBD_analysis.R                     <-- codes for LCBD analysis
    └── FARPROTAX.R                         <-- codes for functional groups analysis

```

# Prokaryotes in thermokarst lakes

> This project repository associated with the following manuscript:

* Luyao Kang, Leiyi Chen, Ziliang Li, Jianjun Wang, Kai Xue, Ye Deng, Manuel Delgado-Baquerizo, Yutong Song, Dianye Zhang, Guibiao Yang, Wei Zhou, Xuning Liu, Futing Liu and Yuanhe Yang* (2023). Patterns and drivers of prokaryotic communities in thermokarst lake water across Northern Hemisphere. (http://doi.org/10.1111/geb.13764 ). 

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
    ├── 1_composition_analysis.R            <-- codes for taxonomic compositional analysis
    ├── 2_cor_DOM_taxa_fun.R                <-- test the relationship between the dominant taxa and DOM properties
    ├── 3_comparison between PA & TP.R      <-- codes for comparison analysis between pan-Arctic and Tibetan Plateau
    ├── 4_LCBD_analysis.R                   <-- codes for LCBD analysis
    └── 5_FARPROTAX.R                       <-- codes for functional groups analysis

```

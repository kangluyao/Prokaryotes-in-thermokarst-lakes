###########################################
=============Tibet data analysis================
###########################################
conda activate qiime2-2020.8

#tibet_sequences
cd tibet-data
#import data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path tibet-fastq-manifest.csv \
  --output-path tibet-sequences.qza \
  --input-format PairedEndFastqManifestPhred33

#trim primers
qiime cutadapt trim-paired \
--i-demultiplexed-sequences tibet-sequences.qza \
--p-front-f GTGCCAGCMGCCGCGGTAA \
--p-front-r GGACTACHVGGGTWTCTAAT \
--o-trimmed-sequences tibet-sequences-trimmed.qza \
--verbose \
&> primer_trimming.log

#quality check using qiime2
qiime demux summarize \
  --i-data tibet-sequences-trimmed.qza \
  --o-visualization tibet-sequences-trimmed.qzv


#Joining reads
qiime vsearch join-pairs \
  --i-demultiplexed-seqs tibet-sequences-trimmed.qza \
  --o-joined-sequences tibet-sequences-trimmed-joined.qza

#quality check using fastqc
unzip tibet-sequences-trimmed-joined.qza -d tibet-sequences-trimmed-joined
gzip -d ./tibet-sequences-trimmed-joined/7db12470-fcb8-473a-80bb-a9ab5da5a26b/data/*.fastq.gz 
cat ./tibet-sequences-trimmed-joined/7db12470-fcb8-473a-80bb-a9ab5da5a26b/data/*.fastq > all.fastq
fastqc all.fastq


#quality control using DADA2 methods
qiime dada2 denoise-paired \
--i-demultiplexed-seqs tibet-sequences-trimmed.qza \
--p-trunc-len-f 225 \
--p-trunc-len-r 210 \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-max-ee-f 2 \
--p-max-ee-r 2 \
--p-trunc-q 2 \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 1 \
--p-n-threads 80 \
--o-table tibet-sequences-dada2-table.qza \
--o-representative-sequences tibet-sequences-dada2-rep-seqs.qza \
--o-denoising-stats tibet-sequences-dada2-stats.qza \
--verbose


#for alpine data, additional individual de novo clustering was performed.
cd tibet-data
qiime vsearch cluster-features-de-novo \
  --i-table sequence-data/tibet-sequences-dada2-table.qza \
  --i-sequences sequence-data/tibet-sequences-dada2-rep-seqs.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table tibet-table-dada2-dn-97.qza \
  --o-clustered-sequences tibet-rep-seqs-dada2-dn-97.qza  

#construct phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences tibet-rep-seqs-dada2-dn-97.qza \
  --p-n-threads 80 \
  --o-alignment aligned-tibet-rep-seqs-dn-97.qza \
  --o-masked-alignment masked-aligned-tibet-rep-seqs-dn-97.qza \
  --o-tree tibet-unrooted-tree-dn-97.qza \
  --o-rooted-tree tibet-rooted-tree-dn-97.qza \
  --verbose

cd ..

#import sequence and taxon file
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path Silva/SILVA_138_QIIME_release/rep_set/rep_set_16S_only/99/silva_138_99_16S.fna \
--output-path 99-meta-data/ref_99_seqs_16S.qza


qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path Silva/SILVA_138_QIIME_release/taxonomy/16S_only/99/majority_taxonomy_7_levels.txt \
--output-path 99-meta-data/ref_99_otus_16S_taxonomy.qza

#Extract reference reads
qiime feature-classifier extract-reads \
  --i-sequences ref_99_seqs_16S.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 250 \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads ref-seqs_16S_V4.qza


#Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs_16S_V4.qza \
  --i-reference-taxonomy ref_99_otus_16S_taxonomy.qza \
  --o-classifier classifier_v4.qza
cd ..

#==========Test the classifier===================
qiime feature-classifier classify-sklearn \
  --i-classifier classifier_v4.qza \
  --i-reads tibet-rep-seqs-dada2-dn-97.qza  \
  --o-classification tibet-rep-seqs-dada2-dn-97-taxonomy.qza


#===========otu table filter===============

#===Remove chloroplast and mitochondria=====   
qiime taxa filter-table \
--i-table tibet-table-dada2-dn-97.qza \
--i-taxonomy tibet-rep-seqs-dada2-dn-97-taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-table tibet-table-nomito-nochloro.qza


#==remove feature frequency less than 10 ====

qiime feature-table filter-features \
--i-table tibet-table-nomito-nochloro.qza \
--p-min-frequency 10 \
--o-filtered-table tibet-table-nomito-nochloro-frequency10.qza

qiime feature-table summarize \
--i-table tibet-table-nomito-nochloro-frequency10.qza \
--o-visualization tibet-table-nomito-nochloro-frequency10.qzv

#=========	Rarefaction	===============
qiime feature-table rarefy \
--i-table tibet-table-nomito-nochloro-frequency10.qza \
--p-sampling-depth 33859 \
--o-rarefied-table table-filtered-mindepth.qza

#=====	Export data =========
mkdir data_from_qiime2
qiime tools export \
--input-path table-filtered-mindepth.qza \
--output-path data_from_qiime2

# Convert biom format to tab-separated text format:
biom convert \
-i data_from_qiime2/feature-table.biom \
-o data_from_qiime2/otu_table.tsv \
--to-tsv

# Export representative sequences:
qiime tools export \
--input-path tibet-rep-seqs-dada2-dn-97.qza \
--output-path data_from_qiime2

# Export taxonomy table:
qiime tools export \
--input-path tibet-rep-seqs-dada2-dn-97-taxonomy.qza \
--output-path data_from_qiime2

# Export tree file:
qiime tools export \
--input-path tibet-rooted-tree-dn-97.qza \
--output-path data_from_qiime2



###########################################
=============Meta analysis=================
###########################################

#ENA_PRJEB40506
cd ENA_PRJEB40506

#import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ENA_PRJEB40506_manifest.csv \
  --output-path ENA_PRJEB40506.qza \
  --input-format SingleEndFastqManifestPhred33

#trim primers
qiime cutadapt trim-single \
  --i-demultiplexed-sequences ENA_PRJEB40506.qza \
  --p-adapter AACMGGATTAGATACCCKG...AGGGTTGCGCTCGTTG \
  --p-match-read-wildcards True \
  --p-discard-untrimmed True \
  --o-trimmed-sequences ENA_PRJEB40506-trimmed-primers.qza \
  --verbose \
&> primer_trimming.log

#quality check using qiime2
qiime demux summarize \
  --i-data ENA_PRJEB40506-trimmed-primers.qza \
  --o-visualization ENA_PRJEB40506-trimmed-primers.qzv

#quality control using DADA2 methods
qiime dada2 denoise-single \
--i-demultiplexed-seqs ENA_PRJEB40506.qza \
--p-trunc-len 280 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-pooling-method independent \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 1 \
--p-n-threads 1 \
--o-table ENA_PRJEB40506-dada2-table.qza \
--o-representative-sequences ENA_PRJEB40506-dada2-rep-seqs.qza \
--o-denoising-stats ENA_PRJEB40506-dada2-stats.qza \
--verbose

qiime feature-table summarize \
  --i-table ENA_PRJEB40506-dada2-table.qza \
  --o-visualization ENA_PRJEB40506-dada2-table.qzv 

cd..

#NCBI_PRJEB36731
cd NCBI_PRJEB36731

#import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path NCBI_PRJEB36731_manifest.csv \
  --output-path NCBI_PRJEB36731.qza \
  --input-format SingleEndFastqManifestPhred33
  
#trim peimers
qiime cutadapt trim-single \
  --i-demultiplexed-sequences NCBI_PRJEB36731.qza \
  --p-adapter GTGYCAGCMGCCGCGGTA...CCCCGYCAATTCMTTTRAGT \
  --p-match-read-wildcards True \
  --p-discard-untrimmed True \
  --o-trimmed-sequences NCBI_PRJEB36731-trimmed-primers.qza \
  --verbose \
&> primer_trimming.log

#quality check using qiime2
qiime demux summarize \
  --i-data NCBI_PRJEB36731-trimmed-primers.qza \
  --o-visualization NCBI_PRJEB36731-trimmed-primers.qzv

#quality control using DADA2 methods
qiime dada2 denoise-single \
--i-demultiplexed-seqs NCBI_PRJEB36731-trimmed-primers.qza \
--p-trunc-len 230 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-pooling-method independent \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 1 \
--p-n-threads 1 \
--o-table NCBI_PRJEB36731-dada2-table.qza \
--o-representative-sequences NCBI_PRJEB36731-dada2-rep-seqs.qza \
--o-denoising-stats NCBI_PRJEB36731-dada2-stats.qza \
--verbose

qiime feature-table summarize \
  --i-table NCBI_PRJEB36731-dada2-table.qza \
  --o-visualization NCBI_PRJEB36731-dada2-table.qzv  
cd ..

#ENA_PRJEB25188
cd ENA_PRJEB25188

#import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ENA_PRJEB25188_manifest.csv \
  --output-path ENA_PRJEB25188.qza \
  --input-format SingleEndFastqManifestPhred33
  
#trim peimers
qiime cutadapt trim-single \
  --i-demultiplexed-sequences ENA_PRJEB25188.qza \
  --p-adapter GTGYCAGCMGCCGCGGTAA...CCGYCAATTYMTTTRAGTTT \
  --p-match-read-wildcards True \
  --p-discard-untrimmed True \
  --o-trimmed-sequences ENA_PRJEB25188-trimmed-primers.qza \
  --verbose \
&> primer_trimming.log

#quality check using qiime2
qiime demux summarize \
  --i-data ENA_PRJEB25188-trimmed-primers.qza \
  --o-visualization ENA_PRJEB25188-trimmed-primers.qzv

#quality control using DADA2 methods
qiime dada2 denoise-single \
--i-demultiplexed-seqs ENA_PRJEB25188-trimmed-primers.qza \
--p-trunc-len 270 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-pooling-method independent \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 1 \
--p-n-threads 1 \
--o-table ENA_PRJEB25188-dada2-table.qza \
--o-representative-sequences ENA_PRJEB25188-dada2-rep-seqs.qza \
--o-denoising-stats ENA_PRJEB25188-dada2-stats.qza \
--verbose

qiime feature-table summarize \
  --i-table ENA_PRJEB25188-dada2-table.qza \
  --o-visualization ENA_PRJEB25188-dada2-table.qzv 
cd ..


#ENA_PRJEB18117
cd ENA_PRJEB18117

#import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ENA_PRJEB18117_manifest.csv \
  --output-path ENA_PRJEB18117.qza \
  --input-format SingleEndFastqManifestPhred33
  
#trim peimers
qiime cutadapt trim-single \
  --i-demultiplexed-sequences ENA_PRJEB18117.qza \
  --p-adapter CCTACGGGNGGCWGCAG...GACTACHVGGGTATCTAATCC \
  --p-match-read-wildcards True \
  --p-discard-untrimmed True \
  --o-trimmed-sequences ENA_PRJEB18117-trimmed-primers.qza \
  --verbose \
&> primer_trimming.log

#quality check using qiime2
qiime demux summarize \
  --i-data ENA_PRJEB18117-trimmed-primers.qza \
  --o-visualization ENA_PRJEB18117-trimmed-primers.qzv

#quality control using DADA2 methods
qiime dada2 denoise-single \
--i-demultiplexed-seqs ENA_PRJEB18117-trimmed-primers.qza \
--p-trunc-len 230 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-pooling-method independent \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 1 \
--p-n-threads 1 \
--o-table ENA_PRJEB18117-dada2-table.qza \
--o-representative-sequences ENA_PRJEB18117-dada2-rep-seqs.qza \
--o-denoising-stats ENA_PRJEB18117-dada2-stats.qza \
--verbose

qiime feature-table summarize \
  --i-table ENA_PRJEB18117-dada2-table.qza \
  --o-visualization ENA_PRJEB18117-dada2-table.qzv  
cd ..



#NCBI_SRP044372
cd NCBI_SRP044372

#import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path NCBI_SRP044372_manifest.csv \
  --output-path NCBI_SRP044372.qza \
  --input-format SingleEndFastqManifestPhred33

qiime demux summarize \
  --i-data NCBI_SRP044372.qza \
  --o-visualization NCBI_SRP044372.qzv
  
#trim peimers
qiime cutadapt trim-single \
  --i-demultiplexed-sequences NCBI_SRP044372.qza \
  --p-adapter ACGCGHNRAACCTTACC...ACGGGCRGTGWGTRCAA \
  --p-match-read-wildcards True \
  --p-discard-untrimmed True \
  --o-trimmed-sequences NCBI_SRP044372-trimmed-primers.qza \
  --verbose \
&> primer_trimming.log

#quality check using qiime2
qiime demux summarize \
  --i-data NCBI_SRP044372-trimmed-primers.qza \
  --o-visualization NCBI_SRP044372-trimmed-primers.qzv

#quality control using DADA2 methods

	#qiime dada2 denoise-pyro is used for Roche454 and Ion Torrent sequencing platform

qiime dada2 denoise-pyro \
--i-demultiplexed-seqs NCBI_SRP044372-trimmed-primers.qza \
--p-trunc-len 300 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-max-len 0 \
--p-chimera-method consensus \
--p-n-threads 0 \
--o-table NCBI_SRP044372-dada2-table.qza \
--o-representative-sequences NCBI_SRP044372-dada2-rep-seqs.qza \
--o-denoising-stats NCBI_SRP044372-dada2-stas.qza \
--verbose

qiime feature-table summarize \
  --i-table NCBI_SRP044372-dada2-table.qza \
  --o-visualization NCBI_SRP044372-dada2-table.qzv 

cd ..


#NCBI_SRS476588
cd NCBI_SRS476588
#import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path NCBI_SRS476588_manifest.csv \
  --output-path NCBI_SRS476588.qza \
  --input-format SingleEndFastqManifestPhred33

#trim peimers
qiime cutadapt trim-single \
  --i-demultiplexed-sequences NCBI_SRS476588.qza \
  --p-adapter ACGCGHNRAACCTTACC...ACGGGCRGTGWGTRCAA \
  --p-match-read-wildcards True \
  --p-discard-untrimmed True \
  --o-trimmed-sequences NCBI_SRS476588-trimmed-primers.qza \
  --verbose \
&> primer_trimming.log

#quality check using qiime2
qiime demux summarize \
  --i-data NCBI_SRS476588-trimmed-primers.qza \
  --o-visualization NCBI_SRS476588-trimmed-primers.qzv
 
#quality control using DADA2 methods
qiime dada2 denoise-pyro \
--i-demultiplexed-seqs NCBI_SRS476588-trimmed-primers.qza \
--p-trunc-len 300 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-max-len 0 \
--p-chimera-method consensus \
--p-n-threads 0 \
--o-table NCBI_SRS476588-dada2-table.qza \
--o-representative-sequences NCBI_SRS476588-dada2-rep-seqs.qza \
--o-denoising-stats NCBI_SRS476588-dada2-stas.qza \
--verbose

qiime feature-table summarize \
  --i-table NCBI_SRS476588-dada2-table.qza \
  --o-visualization NCBI_SRS476588-dada2-table.qzv 
  
cd ..

#ENA_PRJEB30970
cd ENA_PRJEB30970

#import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ENA_PRJEB30970_manifest.csv \
  --output-path ENA_PRJEB30970.qza \
  --input-format SingleEndFastqManifestPhred33
  
#trim peimers
qiime cutadapt trim-single \
  --i-demultiplexed-sequences ENA_PRJEB30970.qza \
  --p-adapter CCTACGGGNGGCWGCAG...GACTACHVGGGTATCTAATCC \
  --p-match-read-wildcards True \
  --p-discard-untrimmed True \
  --o-trimmed-sequences ENA_PRJEB30970-trimmed-primers.qza \
  --verbose \
&> primer_trimming.log

#quality check using qiime2
qiime demux summarize \
  --i-data ENA_PRJEB30970-trimmed-primers.qza \
  --o-visualization ENA_PRJEB30970-trimmed-primers.qzv
 
#quality filter
qiime quality-filter q-score \
 --i-demux ENA_PRJEB30970-trimmed-primers.qza \
 --o-filtered-sequences ENA_PRJEB30970-trimmed-primers-filter.qza \
 --o-filter-stats ENA_PRJEB30970-trimmed-primers-filter-stats.qza \
 --verbose \
&> quality-filter-q-score.log

#quality check using qiime2
qiime demux summarize \
  --i-data ENA_PRJEB30970-trimmed-primers-filter.qza \
  --o-visualization ENA_PRJEB30970-trimmed-primers-filter.qzv

#quality control using DADA2 methods
qiime dada2 denoise-pyro \
--i-demultiplexed-seqs ENA_PRJEB30970-trimmed-primers.qza \
--p-trunc-len 300 \
--p-trim-left 15 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-max-len 0 \
--p-chimera-method consensus \
--p-n-threads 0 \
--o-table ENA_PRJEB30970-dada2-table.qza \
--o-representative-sequences ENA_PRJEB30970-dada2-rep-seqs.qza \
--o-denoising-stats ENA_PRJEB30970-dada2-stas.qza \
--verbose

qiime feature-table summarize \
  --i-table ENA_PRJEB30970-dada2-table.qza \
  --o-visualization ENA_PRJEB30970-dada2-table.qzv
  
cd ..


#merge represent sequences and features tables
qiime feature-table merge \
  --i-tables ENA_PRJEB40506/ENA_PRJEB40506-dada2-table.qza \
  --i-tables ENA_PRJEB25188/ENA_PRJEB25188-dada2-table.qza \
  --i-tables NCBI_PRJEB36731/NCBI_PRJEB36731-dada2-table.qza \
  --i-tables ENA_PRJEB18117/ENA_PRJEB18117-dada2-table.qza \
  --i-tables NCBI_SRP044372/NCBI_SRP044372-dada2-table.qza \
  --i-tables NCBI_SRS476588/NCBI_SRS476588-dada2-table.qza \
  --i-tables ENA_PRJEB30970/ENA_PRJEB30970-dada2-table.qza \
  --i-tables tibet-data/tibet-sequences-dada2-table.qza \
  --o-merged-table thermokarst-lakes-table.qza


qiime feature-table merge-seqs \
  --i-data ENA_PRJEB40506/ENA_PRJEB40506-dada2-rep-seqs.qza \
  --i-data ENA_PRJEB25188/ENA_PRJEB25188-dada2-rep-seqs.qza \
  --i-data NCBI_PRJEB36731/NCBI_PRJEB36731-dada2-rep-seqs.qza \
  --i-data ENA_PRJEB18117/ENA_PRJEB18117-dada2-rep-seqs.qza \
  --i-data NCBI_SRP044372/NCBI_SRP044372-dada2-rep-seqs.qza \
  --i-data NCBI_SRS476588/NCBI_SRS476588-dada2-rep-seqs.qza \
  --i-data ENA_PRJEB30970/ENA_PRJEB30970-dada2-rep-seqs.qza \
  --i-data tibet-data/tibet-sequences-dada2-rep-seqs.qza \
  --o-merged-data thermokarst-lakes-rep-seqs.qza

#species assigned
#import sequence and taxonomy file
cd ..
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path Silva/SILVA_138_QIIME_release/rep_set/rep_set_16S_only/99/silva_138_99_16S.fna \
--output-path 99-meta-data/ref_99_seqs_16S.qza


qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path Silva/SILVA_138_QIIME_release/taxonomy/16S_only/99/majority_taxonomy_7_levels.txt \
--output-path 99-meta-data/ref_99_otus_16S_taxonomy.qza

#for global data, closed-reference mapping protocol was performed.
#Closed-reference clustering
qiime vsearch cluster-features-closed-reference \
  --i-table meta_data/thermokarst-lakes-table.qza \
  --i-sequences meta_data/thermokarst-lakes-rep-seqs.qza \
  --i-reference-sequences 99-meta-data/ref_99_seqs_16S.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table 99-meta-data/thermokarst-lakes-table-97.qza \
  --o-clustered-sequences 99-meta-data/thermokarst-lakes-rep-seqs-97.qza \
  --o-unmatched-sequences 99-meta-data/thermokarst-lakes-rep-seqs-unmatched-97.qza \
  --verbose
  

#Constructing phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 99-meta-data/thermokarst-lakes-rep-seqs-97.qza \
  --o-alignment 99-meta-data/aligned-thermokarst-lakes-rep-seqs-97.qza \
  --o-masked-alignment 99-meta-data/masked-aligned-thermokarst-lakes-rep-seqs-97.qza \
  --o-tree 99-meta-data/thermokarst-lakes-unrooted-tree-dn-97.qza \
  --o-rooted-tree 99-meta-data/thermokarst-lakes-rooted-tree-dn-97.qza \
  --verbose

cd ..

#Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads 99-meta-data/ref_99_seqs_16S.qza \
  --i-reference-taxonomy 99-meta-data/ref_99_otus_16S_taxonomy.qza \
  --o-classifier 99-meta-data/SILVA-138-SSURef-full-length-classifier.qza

  
#=== Test the classifier ====
qiime feature-classifier classify-sklearn \
  --i-classifier SILVA-138-SSURef-full-length-classifier.qza \
  --i-reads results/thermokarst-lakes-rep-seqs-97.qza \
  --o-classification results/thermokarst-lakes-dn-97-taxonomy.qza



#===========	otu table filter	===============

mkdir -p results/filter_data

cd results
qiime feature-table summarize \
--i-table thermokarst-lakes-table-97.qza \
--o-visualization thermokarst-lakes-table-97.qzv



###Remove chloroplast and mitochondria:   
qiime taxa filter-table \
--i-table thermokarst-lakes-table-97.qza \
--i-taxonomy thermokarst-lakes-dn-97-taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-table filter_data/thermokarst-lakes-table-nomito-nochloro.qza


###remove feature frequency less than 10

qiime feature-table filter-features \
--i-table filter_data/thermokarst-lakes-table-nomito-nochloro.qza \
--p-min-frequency 10 \
--o-filtered-table filter_data/thermokarst-lakes-table-nomito-nochloro-frequency10.qza

qiime feature-table summarize \
--i-table filter_data/thermokarst-lakes-table-nomito-nochloro-frequency10.qza \
--o-visualization filter_data/thermokarst-lakes-table-nomito-nochloro-frequency10.qzv

# minum frequency is 1070




#=========	Rarefaction	===============
qiime feature-table rarefy \
--i-table filter_data/thermokarst-lakes-table-nomito-nochloro-frequency10.qza \
--p-sampling-depth 1070 \
--o-rarefied-table filter_data/table-filtered-rarefy.qza


#=====	Export data =========
mkdir data_from_qiime2
mkdir -p filter_data/data_from_qiime2_unrarefied
qiime tools export \
--input-path filter_data/table-filtered-rarefy.qza \
--output-path filter_data/data_from_qiime2

# Convert biom format to tab-separated text format:
biom convert \
-i filter_data/data_from_qiime2/feature-table.biom \
-o filter_data/data_from_qiime2/otu_table.tsv \
--to-tsv


# Export representative sequences:
qiime tools export \
--input-path thermokarst-lakes-rep-seqs-97.qza \
--output-path filter_data/data_from_qiime2

# Export taxonomy table:
qiime tools export \
--input-path thermokarst-lakes-dn-97-taxonomy.qza \
--output-path filter_data/data_from_qiime2

# Export tree file:
qiime tools export \
--input-path tree/thermokarst-lakes-rooted-tree-dn-97.qza \
--output-path filter_data/data_from_qiime2


#=====	Export unrarefied data =========
mkdir -p filter_data/data_from_qiime2_unrarefied
qiime tools export \
--input-path filter_data/thermokarst-lakes-table-nomito-nochloro-frequency10.qza \
--output-path filter_data/data_from_qiime2_unrarefied

# Convert biom format to tab-separated text format:
biom convert \
-i filter_data/data_from_qiime2_unrarefied/feature-table.biom \
-o filter_data/data_from_qiime2_unrarefied/otu_table.tsv \
--to-tsv


# Export representative sequences:
qiime tools export \
--input-path thermokarst-lakes-rep-seqs-97.qza \
--output-path filter_data/data_from_qiime2_unrarefied

# Export taxonomy table:
qiime tools export \
--input-path thermokarst-lakes-dn-97-taxonomy.qza \
--output-path filter_data/data_from_qiime2_unrarefied

# Export tree file:
qiime tools export \
--input-path tree/thermokarst-lakes-rooted-tree-dn-97.qza \
--output-path filter_data/data_from_qiime2_unrarefied

#close the qiime2
conda deactivate

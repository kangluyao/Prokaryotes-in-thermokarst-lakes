setwd('E:/Tax4Fun2')
# Delete previous versions of Tax4Fun2
if("Tax4Fun2" %in% installed.packages()) remove.packages("Tax4Fun2")

# Download and install Tax4Fun2
download.file(url = "https://github.com/bwemheu/Tax4Fun2/releases/download/v1.1.6/Tax4Fun2_1.1.6.tar.gz", destfile = "Tax4Fun2_1.1.6.tar.gz", mode = "wb")
install.packages("Tax4Fun2_1.1.6.tar.gz", repos = NULL, source = T)
unlink("Tax4Fun2_1.1.6.tar.gz")

# Build the reference data and dependencies
library(Tax4Fun2)
buildReferenceData(path_to_working_directory = ".", install_suggested_packages = T)
buildDependencies(path_to_reference_data = "Tax4Fun2_ReferenceData_v2")

# Many people had issues with blastm on Windows
# Test you installation
if(tolower(Sys.info()[1]) == "windows") system("Tax4Fun2_ReferenceData_v2/blast_bin/bin/blastn.exe -version")

# Make a prediction
# Run the reference blast first
runRefBlast(path_to_otus ="E:/thermokast_lakes/water_microbes/meta_analysis/data/meta_data/meta_ref_seqs.fasta", 
              path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
              path_to_temp_folder = "E:/thermokast_lakes/water_microbes/meta_analysis/data/meta_data/Pred",
              database_mode = "Ref99NR",
              num_threads = 12, use_force = T)

# Calculate functional predictions
makeFunctionalPrediction(path_to_otu_table = "E:/thermokast_lakes/water_microbes/meta_analysis/data/meta_data/meta_otu_table.txt",
                           path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
                           path_to_temp_folder = "E:/thermokast_lakes/water_microbes/meta_analysis/data/meta_data/Pred",
                           database_mode = "Ref99NR",
                           normalize_by_copy_number = TRUE,
                           min_identity_to_reference = 0.97)


# Make a prediction
# Run the reference blast first
runRefBlast(path_to_otus ="E:/thermokast_lakes/water_microbes/tibet_dada2_asv/data/total_taxa_data/ref.seqs.fasta", 
            path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
            path_to_temp_folder = "E:/thermokast_lakes/water_microbes/tibet_dada2_asv/data/total_taxa_data/Pred",
            database_mode = "Ref99NR",
            num_threads = 12, use_force = T)

# Calculate functional predictions
makeFunctionalPrediction(path_to_otu_table = "E:/thermokast_lakes/water_microbes/tibet_dada2_asv/data/total_taxa_data/otu.table.txt",
                         path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
                         path_to_temp_folder = "E:/thermokast_lakes/water_microbes/tibet_dada2_asv/data/total_taxa_data/Pred",
                         database_mode = "Ref99NR",
                         normalize_by_copy_number = TRUE,
                         min_identity_to_reference = 0.97)

# Make a prediction
# Run the reference blast first
runRefBlast(path_to_otus ="E:/thermokast_lakes/water_microbes/OTU_97/data/module_13_taxa_data/ref.seqs.fasta", 
            path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
            path_to_temp_folder = "Pred",
            database_mode = "Ref99NR",
            num_threads = 12, use_force = T)

# Calculate functional predictions
makeFunctionalPrediction(path_to_otu_table = "E:/thermokast_lakes/water_microbes/OTU_97/data/module_13_taxa_data/otu.table.txt",
                         path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
                         path_to_temp_folder = "Pred",
                         database_mode = "Ref99NR",
                         normalize_by_copy_number = TRUE,
                         min_identity_to_reference = 0.97)


# Calculate functional redundnacy (the doc file needs to be updates. bug)
calculateFunctionalRedundancy(path_to_otu_table = "./module_data/tax4fun2/otu_table.txt", 
                              path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
                              path_to_temp_folder = "Pred", database_mode = "Ref99NR", 
                              min_identity_to_reference = 0.97)

# Have fun

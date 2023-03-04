#diff class using microeco package
library(microeco)
ps_genus <- phyloseq::tax_glom(meta_physeq, taxrank = 'Genus')
meco_genus_df <- phyloseq2meco(ps_genus)
#calculate the abundance table
m1_genus <- meco_genus_df$cal_abund()

m1_genus <- trans_diff$new(dataset = meco_genus_df, method = "lefse", 
                     group = "Region", alpha = 0.01, 
                     lefse_subgroup = NULL)
# t1$res_lefse is the LEfSe result
# t1$res_abund is the abundance information
m1_genus$plot_lefse_bar(LDA_score = 4)

m1_genus$plot_diff_abund(use_number = 1:30)


# clade_label_level 5 represent phylum level in this analysis
# require ggtree package
use_labels <- c("Firmicutes", "Actinobacteria", "Gammaproteobacteria", 
                "Alphaproteobacteria","Bacteroidetes", 
                "Verrucomicrobiales")

m1_genus$plot_lefse_cladogram(use_taxa_num = 200, use_feature_num = 50, 
                              select_show_labels = use_labels)
res_lefse <- write.table(m1$res_lefse, file = './meta_analysis/results/tables/res_lefse.txt')




m1_genus$plot_lefse_cladogram(use_taxa_num = 200,
                              filter_taxa = NULL,
                              use_feature_num = 50,
                              group_order = NULL,
                              clade_label_level = 4,
                              select_show_labels = use_labels,
                              only_select_show = F,
                              sep = "|",
                              branch_size = 0.2,
                              alpha = 0.2,
                              clade_label_size = 0.7,
                              node_size_scale = 1,
                              node_size_offset = 1,
                              annotation_shape = 22,
                              annotation_shape_size = 5) 
                              






ps_genus_sel <- subset_taxa(ps_genus, Order == "Caulobacterales" | 
                                    Family == "Beijerinckiaceae" | 
                                    Family == "Sphingomonadaceae" | 
                                    Family == "Burkholderiaceae" | 
                                    Class == "Verrucomicrobiae" |
                                    Family == "Frankiales" |
                                    Genus == "CandidatusAquiluna" |
                                    Genus == "CandidatusLimnoluna" |
                                    Genus == "CandidatusPlanktoluna" |
                                    Family == "Chitinophagaceae" |
                                    Order == "Sphingobacteriales" |
                                    Family == "Flavobacteriaceae" |
                                    Family == "Spirosomaceae" |
                                    Family == "Cyclobacteriaceae" |
                                    Class == "Bacilli")

phyloseq2lefse <- function(
        ps,
        covars,
        file.name = "lefse_data.txt",
        taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
        transpose.otus = TRUE
) {
        
        # grab the taxa table from the phyloseq object and coerce into a matrix
        tax.tbl <- as(phyloseq::tax_table(ps), "matrix")
        tax.tbl <- tax.tbl[, taxa.levels]
        tax.tbl.join <- as.data.frame(do.call(paste, c(as.data.frame(tax.tbl), sep = "|")))
        
        row.names(tax.tbl.join) <- row.names(tax.tbl)
        colnames(tax.tbl.join) <- NULL
        
        # grab the otu table from the phyloseq object
        otu.tbl <- as(phyloseq::otu_table(ps), "matrix")
        colnames(otu.tbl) <- NULL
        tax.otu.tbl <- cbind(tax.tbl.join, otu.tbl)
        
        # grab the sample table from the phyloseq object
        smpl.data <- as(phyloseq::sample_data(ps), "data.frame")
        smpl.data$Sample <- row.names(smpl.data)
        t.smpl.data <- t(smpl.data)
        t.smpl.data <- as.matrix(t.smpl.data[c("Sample", covars), ])
        t.smpl.data <- cbind(rownames(t.smpl.data), t.smpl.data)
        colnames(t.smpl.data) <- colnames(tax.otu.tbl)
        
        final.data <- rbind(t.smpl.data, tax.otu.tbl)
        write.table(final.data, file = file.name, row.names = F, col.names = F, sep = "\t", quote = FALSE)
}


phyloseq2lefse(
        ps = ps_genus_sel,
        covars = 'Region',
        file.name = "./meta_analysis/results/tables/lefse_data_sel_genus.txt",
        taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
        transpose.otus = TRUE
)




tp_diversity <- subset(meta_diversity, Region == 'Tibetan Plateau')
write.csv(tp_diversity, file = './tp_diversity.csv')


genus_tree <- phy_tree(ps_genus)
write.tree(genus_tree, "./meta_analysis/results/tables/lefse_tree.nwk")
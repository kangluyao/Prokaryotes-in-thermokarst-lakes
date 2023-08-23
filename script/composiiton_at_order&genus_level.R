# plot the community composition for TP and PA at order level

subphylo.rel <- tax_glom(meta_physeq_rel, taxrank = 'Order')
ntaxa(subphylo.rel)
# meta_physeq_rel <- microbiome::transform(meta_physeq, "compositional")
# meta.com.cla <- microbiome::aggregate_rare(meta_physeq_rel, level = "Order", 
#                                            detection = 1/100, prevalence = 10/100)

ra.tab <- otu_table(subphylo.rel)
sum(ra.tab[, 1])
subtaxa_tab <- tax_table(subphylo.rel)[, 4]

order_tax_table <- data.frame(subtaxa_tab, ra.tab) %>% 
  mutate(MRA = rowMeans(select(., colnames(ra.tab)))) %>%
  arrange(desc(MRA)) %>% dplyr::top_n(15, MRA) %>%
  select(., -c('MRA')) %>% 
  bind_rows(summarise_all(., ~if(is.numeric(.)) 1-sum(.) else "Others")) %>%
  mutate(Order = factor(Order, levels = Order)) %>%
  pivot_longer(cols = -c(Order), names_to = "Sample_name", values_to = 'rel_abun') %>%
  right_join(data.frame(Sample_name = rownames(metadata), Region = metadata$Region), by = c("Sample_name")) %>%
  select(., -c('Sample_name')) %>% 
  group_by(Region, Order) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))

print(order_tax_table, n = 32)


# plot the community composition for TP and PA at order level

subphylo.rel <- tax_glom(meta_physeq_rel, taxrank = 'Genus')
ntaxa(subphylo.rel)

ra.tab <- otu_table(subphylo.rel)
sum(ra.tab[, 1])
subtaxa_tab <- tax_table(subphylo.rel)[, 6]

order_tax_table <- data.frame(subtaxa_tab, ra.tab) %>% 
  mutate(MRA = rowMeans(select(., colnames(ra.tab)))) %>%
  arrange(desc(MRA)) %>% dplyr::top_n(15, MRA) %>%
  select(., -c('MRA')) %>% 
  bind_rows(summarise_all(., ~if(is.numeric(.)) 1-sum(.) else "Others")) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  pivot_longer(cols = -c(Genus), names_to = "Sample_name", values_to = 'rel_abun') %>%
  right_join(data.frame(Sample_name = rownames(metadata), Region = metadata$Region), by = c("Sample_name")) %>%
  select(., -c('Sample_name')) %>% 
  group_by(Region, Genus) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))

print(order_tax_table, n = 32)

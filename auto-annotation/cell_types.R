cell_data <- readr::read_tsv("data/hpa/kidney/cell_data.tsv")

cell_data %>% ggplot(aes(x=umap_x, y=umap_y)) +
  geom_point(aes(color=cluster))


cell_type_tissue <- readr::read_tsv("data/hpa/rna_single_cell_type_tissue.tsv") %>%
  rename(gene_id = Gene, symbol = `Gene name`, cell_type = `Cell type`,
         read_count = `Read count`)

cell_type_kidney <- cell_type_tissue %>%
  filter(Tissue=="kidney")


kidney_cell_expn <- targets %>%
  left_join(
    cell_type_tissue,
#    group_by(gene_id, symbol, Tissue, cell_type) %>%
#    summarize(read_count = sum(read_count), nTPM = mean(nTPM)) %>%
#      ungroup(),
    by=c("gene_id", "symbol"))

kidney_cell_expn %>%
  ggplot(aes(Cluster, symbol)) + geom_tile(aes(fill = nTPM),
                                             colour = "white") + scale_fill_gradient(low = "white",
                                                                                     high = "steelblue")+
  labs(x='', y = '')+
  theme_bw()+
  guides(fill = guide_legend(title = "Log2(TPM)"))+
  #theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())

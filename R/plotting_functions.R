gtex_boxplot <- function(dataset, focal_gene_id, focal_symbol, tissue) {
    dataset %>%
        dplyr::filter(gene_id == focal_gene_id) %>%
        ggplot(aes(x=tissue, y=value, group=tissue)) +
        geom_boxplot() +
        ggtitle(focal_symbol) +
        ylab("GTEx RPKMs") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

kpmp_boxplot <- function(dataset, focal_gene_id, focal_symbol, tissue) {
    dataset %>%
        dplyr::filter(symbol==focal_symbol) %>%
        ggplot(aes(x = cell_type, y = cpm, fill = structure)) +
        geom_boxplot() +
        ggtitle(paste(focal_symbol,
                      "gene expression in",
                      gene_annotation$metadata$tissue,
                      "from KPMP")) +
        xlab("Cell type") +
        ylab("log(CPM)") +
        theme(plot.title = element_text(hjust = 0,size = 40),
              axis.title = element_text(size=15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 24),
              axis.text.y = element_text(size = 24),
              legend.key.size = unit(60, 'points'),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 18),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
}

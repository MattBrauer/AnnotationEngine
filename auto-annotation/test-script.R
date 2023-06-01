library(TissueEnrich)
library(biomaRt)
library(tidyverse)

gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
output<-teEnrichment(inputGenes = gs)

gtex <- readr::read_tsv("data/hpa/rna_tissue_gtex.tsv") %>%
  rename(gene_id = Gene, symbol = `Gene name`)
hpa <- readr::read_tsv("data/hpa/rna_tissue_hpa.tsv") %>%
  rename(gene_id = Gene, symbol = `Gene name`)
consensus <- readr::read_tsv("data/hpa/rna_tissue_consensus.tsv") %>%
  rename(gene_id = Gene, symbol = `Gene name`)

normal <- readr::read_tsv("data/hpa/normal_tissue.tsv") %>%
  rename(gene_id = Gene, symbol = `Gene name`, cell_type = `Cell type`) %>%
  mutate(level = if_else(Level %in% c("Medium","Low","High"), Level, "ND"),
         level = ordered(level, levels = c("ND", "Low", "Medium", "High")),
         value = as.integer(level) - 1) %>%
  unite(tissue_cell, Tissue:cell_type, remove = FALSE)

targets %>%
  left_join(normal) %>%
  select(symbol, Tissue, cell_type, tissue_cell, value) %>%
  filter(Tissue == "kidney") %>%
  ggplot(aes(symbol, cell_type)) + geom_tile(aes(fill = value),
                                              colour = "white") + scale_fill_gradient(low = "white",
                                                                                      high = "steelblue")+
  labs(x='', y = '')+
  theme_bw()+
  guides(fill = guide_legend(title = "Log2(TPM)"))+
  #theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())





gs <- GeneSet(geneIds = targets$gene_id, organism="Homo Sapiens", geneIdType=ENSEMBLIdentifier())
output <- teEnrichment(inputGenes = gs)

gs <- GeneSet(geneIds = targets$symbol, organism="Homo Sapiens", geneIdType=SymbolIdentifier())
output <- teEnrichment(inputGenes = gs)


seEnrichmentOutput<-output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)
head(enrichmentOutput)

seExp<-output[[2]][["Kidney"]]
exp<-setNames(data.frame(assay(seExp), row.names = rowData(seExp)[,1]), colData(seExp)[,1])
exp$Gene<-row.names(exp)
exp<-exp %>% gather(key = "Tissue", value = "expression",1:(ncol(exp)-1))

tissue_enrichment <- exp %>%
  tibble::as_tibble() %>%
  pivot_wider(names_from = Tissue, values_from = expression) %>%
  arrange(Gene)

ggplot(exp, aes(Tissue, Gene)) + geom_tile(aes(fill = expression),
                                           colour = "white") + scale_fill_gradient(low = "white",
                                                                                   high = "steelblue")+
  labs(x='', y = '')+
  theme_bw()+
  guides(fill = guide_legend(title = "Log2(TPM)"))+
  #theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())


# mouse
gs<-GeneSet(geneIds=targets$symbol,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
output<-teEnrichment(inputGenes = gs,rnaSeqDataset = 3)
seEnrichmentOutput<-output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)

seExp<-output[[2]][["Kidney"]]
exp<-setNames(data.frame(assay(seExp), row.names = rowData(seExp)[,1]), colData(seExp)[,1])
exp$Gene<-row.names(exp)
exp<-exp %>% gather(key = "Tissue", value = "expression",1:(ncol(exp)-1))
ggplot(exp, aes(Tissue, Gene)) + geom_tile(aes(fill = expression),
                                           colour = "white") + scale_fill_gradient(low = "white",
                                                                                   high = "steelblue")+
  labs(x='', y = '')+
  theme_bw()+
  guides(fill = guide_legend(title = "Log2(TPM)"))+
  #theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())




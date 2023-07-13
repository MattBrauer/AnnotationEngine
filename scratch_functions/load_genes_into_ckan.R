library(biomaRt)
library(tidyverse)
library(ckanr)

sightline_org_id <- "acd65bfd-df6d-4297-ae33-ab837864b022"
sightline_catalog_pkg_id <- "6035361d-871c-4e24-a027-e044fe5d42a7"

mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

query_string <- "select * from sightline.target_concepts"

system.time({res <- noctua::dbExecute(con, query_string)})

output <- noctua::dbFetch(res)
noctua::dbClearResult(res)

allTargets <- output %>%
  tibble::as_tibble() %>%
  dplyr::select(gene, disease) %>%
  rename(gene_id = gene, context = disease)

symbol_map <- getBM(
  mart = mart,
  attributes = c('ensembl_gene_id','hgnc_symbol'),
  filter = 'ensembl_gene_id',
  values = allTargets$gene_id,
  uniqueRows = TRUE) %>%
  tibble::as_tibble() %>%
  rename(gene_id = ensembl_gene_id, symbol = hgnc_symbol)

allTargets <- allTargets %>%
  left_join(symbol_map) %>%
  dplyr::select(gene_id, symbol, context)

renalNewTargets <- toupper(c('Aak1', 'Apoc3', 'Arg1', 'Blvrb', 'Cgnl1', 'Ckm',
                             'Cldn10', 'Col4a', 'Cubn', 'Dab2', 'Fnip1', 'G6pc',
                             'Gatm', 'Ift140', 'Lcn8', 'Lrp2', 'Miox', 'Mitf',
                             'Mtm1', 'Ninl', 'Nphp3', 'Nphs1', 'Pkd1', 'Pkd2',
                             'Pkhd1', 'Rnf128', 'Rpl3l', 'Slc22a1', 'Slc22a2',
                             'Slc25a45', 'Slc34a1', 'Slc34a3', 'Slc47a1',
                             'Slc6a19', 'Slc7a9', 'Svil', 'Umod', 'Smad6'))

newTargets <- getBM(
  mart = mart,
  attributes = c('ensembl_gene_id','hgnc_symbol'),
  filter = 'hgnc_symbol',
  values = renalNewTargets,
  uniqueRows = TRUE) %>%
  tibble::as_tibble() %>%
  rename(gene_id = ensembl_gene_id,
         symbol = hgnc_symbol) %>%
  dplyr::select(gene_id, symbol) %>%
  mutate(context = "Renal New Targets 2023")

targets <- allTargets %>%
  bind_rows(newTargets) %>%
  readr::write_csv("data/resources/sightline_targets.csv")

ckanr::resource_create(package_id = sightline_catalog_pkg_id,
                       name = "genes",
                       description = "genes in the Sightline catalog",
                       upload = "data/resources/sightline_targets.csv",
                       url = ckan_url,
                       key = ckan_key)

genes <- targets %>%
  dplyr::select(gene_id, symbol) %>%
  distinct()

targets <- getBM(
  mart = mart,
  attributes = c('ensembl_gene_id','hgnc_symbol'),
  filter = 'ensembl_gene_id',
  values = target_genes,
  uniqueRows = TRUE) %>%
  tibble::as_tibble() %>%
  rename(gene_id = ensembl_gene_id, symbol = hgnc_symbol)


gtex <- readr::read_tsv("data/hpa/rna_tissue_gtex.tsv") %>%
  rename(gene_id = Gene, symbol = `Gene name`)
hpa <- readr::read_tsv("data/hpa/rna_tissue_hpa.tsv") %>%
  rename(gene_id = Gene, symbol = `Gene name`)
consensus <- readr::read_tsv("data/hpa/rna_tissue_gtex.tsv") %>%
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




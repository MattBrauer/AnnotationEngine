library(Gviz)
domains_edb <- suppressWarnings(proteins(edb_Hs,
                                    filter = GeneIdFilter(genelist$gene_id,
                                                          "=="),
                                    columns = c("tx_id", listColumns(edb_Hs,
                                                                     "protein_domain"))))
txs <- getGeneRegionTrackForGviz(
    edb_Hs, filter = GeneIdFilter(genelist$gene_id,
                                  "=="))

pdoms <- proteins(edb_Hs, filter = ~ tx_id %in% txs$transcript &
                      protein_domain_source == "pfam",
                  columns = c("protein_domain_id", "prot_dom_start",
                              "prot_dom_end"))

pdoms

features.gr <- txs

gene_id <- "ENSG00000115977"
query_string <- paste("select distinct concat(tag_chrom, '_', cast(tag_pos as varchar), '_', tag_ref, '_', tag_alt) as variant_id,",
                      "tag_chrom as chr, tag_pos as pos, tag_ref as ref, tag_alt as alt,",
                      "study_id, overall_score,",
                      "element_at(trait_efos, 1) as trait_efo, trait_reported,",
                      "beta, beta_ci_lower, beta_ci_upper, pval",
                      "from open_targets_genetics.d2v2g_scored where",
                      "pval < 5e-8 and",
                      "concat(tag_chrom, '_', cast(tag_pos as varchar), '_', tag_ref, '_', tag_alt) in",
                      #variant_annotation$metadata$variant_strings[[3]])
                      variant_annotation$metadata$variant_strings[[gene_id]])

res <- noctua::dbExecute(variant_annotation$metadata$con, query_string)

phewas_output <- noctua::dbFetch(res) %>%
    tibble::as_tibble() %>%
#    dplyr::left_join(tibble::enframe(variant_annotation$metadata$variants, name="variant_id", value="gene_id"), by="variant_id") %>%
#    dplyr::left_join(tibble::enframe(variant_annotation$metadata$genes, name="gene_id", value="symbol"), by="gene_id") %>%
    dplyr::group_by(symbol, variant_id) %>%
    dplyr::arrange(pval)
noctua::dbClearResult(res)

phewas_output %>%
    ggplot(aes(x=as.numeric(pos), y=as.factor(trait_reported))) +
    geom_point(aes(size = overall_score), color="red", alpha=0.01)




gr <- GRanges("2", IRanges(start = phewas_output$pos, width = 1, names = phewas_output$))

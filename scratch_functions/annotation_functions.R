#' Annotate set of genes
#'
#' @description
#'
#' @param genes A vector of Ensembl gene_ids
#' @param data_dir The location where output will be placed
#'
#' @return A list containing the annotation artifacts created, as files, in
#'  `data_dir`

#' @export
annotate_genes <- function(geneset, data_dir, session = session) {
    check_geneset()
    check_data_dir()
    
        output <- setNames(lapply(as.list(c(geneset, "geneset")), function(gene) {
        list("tables" = list(), "figures" = list())
        }),
             c(geneset, "geneset"))
    
}





#' Generate HTML table of mouse phenotypes by gene
#'
#' @description
#'
#' @param gene The human gene of interest
#'
#' @return A formatted HTML table of phenotypes for the mouse ortholog
#' @export
mouse_phenotypes <- function(gene) {
  readRDS(fs::path(data_dir, "opentargets/mouse_models.rds")) %>% # need to regenerate after data update
    dplyr::filter(targetFromSourceId == gene) %>%
    unnest(modelPhenotypeClasses) %>%
    dplyr::rename(pheID = id, pheLabel = label) %>%
    unnest(biologicalModels) %>%
    flextable::regulartable() %>%
    flextable::theme_box() %>%
    flextable::merge_v(c("symbol", "targetFromSourceId", "targetInModel",
                         "targetInModelMgiId", "targetInModelEnsemblId",
                         "modelPhenotypeId", "modelPhenotypeLabel")) %>%
    flextable::valign(valign = "top", part = "all") %>%
    flextable::set_table_properties(width = 1, layout = "autofit") %>%
    flextable::set_header_labels(targetFromSourceId = "ENSG",
                                 targetInModel = "mouse gene",
                                 targetInModelMgiId = "mouse gene MGI",
                                 targetInModelEnsemblId = "mouse gene ENS",
                                 modelPhenotypeId = "phenotype",
                                 id = "model ID",
                                 allelicComposition = "allelic composition",
                                 geneticBackground = "genetic background"
    ) %>%
    htmltools_value()
}

#' Generate HTML table of mouse phenotypes by gene
#'
#' @description
#'
#' @param gene The human gene of interest
#'
#' @return A formatted HTML table of phenotypes for the mouse ortholog
#' @export
mouse_phenotypes <- function(gene) {
  readRDS(fs::path(data_dir, "opentargets/mouse_models.rds")) %>%
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

#' Annotate set of genes
#'
#' @description Set up initial gene annotation data structure
#'
#' @param genelist A vector of Ensembl gene_ids
#' @param setname A character string (default "all") describing the gene grouping
#' @param tissue A character string (default "all") specifying the focal tissue
#' @param data_dir The location where output will be placed
#' @param edb An `edb` object, specific to species
#' @param con A DBI connection to the relevant database
#'
#' @return A list containing the annotation artifacts created, as files, in
#'  `data_dir`
#' @export
annotate_genes <- function(genelist,
                           setname = "all",
                           tissue = "all",
                           data_dir,
                           out_dir,
                           edb = edb_Hs,
                           con) {

    genes_edb <- suppressWarnings(genes(edb,
                                        filter = GeneIdFilter(genelist$gene_id,
                                                              "==")))

    warnings <- list()

    data_dir_check <- geneset_check <- FALSE

    if(fs::dir_exists(data_dir)) {
        data_dir_check <- TRUE
    } else {
        warning("data_dir does not exist")
    }

    is_valid <- genelist$gene_id %in% genes_edb$gene_id

    if(all(is_valid)) {
        geneset_check <- TRUE
    } else {
        warning("Not all genes in geneset map to gene_ids")
        warnings[["ensembl_annotation_failures"]] <- genelist$gene_id[!is_valid]
    }
    annotation <- mcols(genes_edb) %>%
        tibble::as_tibble() %>%
        dplyr::select(gene_id, symbol, description,
                      gene_id_version, gene_biotype, canonical_transcript)

    symbols <- annotation %>%
        dplyr::select(gene_id, symbol) %>%
        tibble::deframe()

    gene_annotation <- list(
        metadata = list(ensembl = annotation,
                        tissue = tissue,
                        setname = setname,
                        data_dir = data_dir,
                        out_dir = out_dir,
                        edb = edb,
                        genes = symbols,
                        gene_string = paste0('(',
                                             names(symbols) %>%
                                                 stringr::str_c("'", ., "'") %>%
                                                 stringr::str_flatten_comma(),
                                             ')'),
                        con = con),
        warnings = warnings,
        datasets = list()
            # name, description, resources
        )

    return(gene_annotation)
}

## Functions of class `annotate_gene_` take and return a gene_annotation
## object.
##
## Each function appends a new dataset to `gene_annotation$datasets`
##
## Datasets are lists with the following structure
##  - name, character
##  - description, character
##  - resources, list
##  - published, logical indicating whether the dataset has been published
##  - cleared, logical indicating whether the dataset's resources have been
##    deleted
##
## `resources` is a list, with each element of the following structure:
##  - dataset_name, character, name of dataset to which resource belongs
##  - resource_name, character, name of resource
##  - filename, character, complete local filepath of resource to be uploaded
##  - type, character indicating the file type of the resource
##  - published, logical indicating whether the resource has been published and
##    may now be cleared
##  - cleared, logical indicating whether the uploaded file has been deleted
##
##  Function needs to define the dataset name and description, and all
##    parameters of the resources



#' Generate gene expression annotations
#'
#' @description Generate gene annotation artifacts for tissue-specific expression
#'
#' @param gene_annotation A list as output by 'annotate_genes'
#'
#' @return A list as output by 'annotate_genes', with gene expression annotation
#' artifact locations added
#' @export
annotate_gene_tissue_rna_specificity <- function(gene_annotation) {

    dataset <- list(name = tolower(paste("tissue_rna_expn",
                                         stringr::str_replace_all(gene_annotation$metadata$setname,
                                                                  "[^[a-z][A-Z][0-9][-_]]",
                                                                  "_"), sep ="_")),
                    description = paste("Tissue-specific RNA expression of",
                                                ifelse(gene_annotation$metadata$setname=="all",
                                                       "all genes",
                                                       paste("genes from",
                                                             gene_annotation$metadata$setname))),
                    published = FALSE,
                    cleared = FALSE)
    ## GTEx
    # GTEx samples
    query_string <- paste("select SAMPID as sample, SMTS as Tissue, SMTSD as Subtissue",
                          "from gtex.samples")

    #    res <- noctua::dbExecute(gene_annotation$metadata$con, query_string)
    res <- mazer::athena_query(query_string, gene_annotation$metadata$con)
    samples <- res %>%
        tibble::as_tibble() %>%
        dplyr::rename(tissue = Tissue, subtissue = Subtissue)

    # GTEx TPMs
    query_string <- paste("select split_part(gene_id, '.', 1) as gene_id, symbol, sample, value",
                          "from gtex.v8 where measurement_type='gene_tpm'",
                          "and split_part(gene_id, '.', 1) in",
                          gene_annotation$metadata$gene_string)
    res <- mazer::athena_query(query_string, gene_annotation$metadata$con)
    tpm <- res %>%
        tibble::as_tibble() %>%
        left_join(samples, by="sample")

    ## boxplots of GTEx TPMs by tissue
    dataset$resources <- setNames(lapply(names(gene_annotation$metadata$genes),
                                         function(focal_gene_id) {
        gene <- list()

        focal_symbol <- gene_annotation$metadata$genes[[focal_gene_id]]

        plotfile <- fs::path(gene_annotation$metadata$out_dir,
                             paste0(focal_symbol, "_GTEx_tissue_expn_boxplot.png"))

        p <- gtex_boxplot(tpm, focal_gene_id, focal_symbol, tissue)

        png(plotfile)
        print(p)
        dev.off()

        gene <- list(filename = plotfile,
                     dataset_name = dataset$name,
                     resource_name = paste(focal_symbol, "gtex_tissue_expn_boxplot", sep="_"),
                     type = "png",
                     plot = p,
                     published = FALSE,
                     cleared = FALSE)
        return(gene)
    }),
    gene_annotation$metadata$genes)

    ## heatmap of gtex gene expression (from HPA) across geneset
    query_string <- paste("select gene_id, symbol, tissue, value",
                          "from hpa.gene_expression where measurement='nTPM'",
                          "and source='gtex'", #consensus'",
                          "and gene_id in",
                          gene_annotation$metadata$gene_string)
    #res <- noctua::dbExecute(gene_annotation$metadata$con, query_string)
    res <- mazer::athena_query(query_string, gene_annotation$metadata$con)

    tpm <- res %>%
        tibble::as_tibble() %>%
        dplyr::mutate(tissue = stringr::str_to_title(tissue))

    tpm <- tpm %>%
        left_join(tpm %>%
                      group_by(gene_id) %>%
                      summarize(median_tpm = median(value),
                                sd_tpm = sd(value)),
                  by="gene_id") %>%
        dplyr::mutate(expn = (value - median_tpm)/sd_tpm)

    ## Example use of `tpm`
    # tpm %>% readr::write_csv(fs::path(out_dir, "gtex_expn_renal_new_targets", ext="csv"))
    # ckanr::package_search("tissue_rna_expn_renal_new_targets_2023")
    # ckanr::resource_create(package_id = "d850149c-3113-4b12-801b-eb6c6c76f88f",
    #                        name="gtex_data_renal_new_targets",
    #                        description="GTEx expression data for renal new targets",
    #                        format="CSV", upload=fs::path(out_dir, "gtex_expn_renal_new_targets", ext="csv"))


    ## heatmap of consensus from GTEx and HPA TPMs by tissue
    plotfile <- fs::path(gene_annotation$metadata$out_dir,
                         stringr::str_replace_all(paste0(gene_annotation$metadata$setname,
                                                         "_GTEx_HPA_tissue_expn_heatmap.png"),
                                                  " ", "_"))

    plot_title <- ifelse(gene_annotation$metadata$setname=="all",
                   "Median-normalized tissue expression of all genes",
                   paste("Median-normalized tissue expression of genes from",
                         gene_annotation$metadata$setname))

    p <- tpm %>%
        dplyr::select(symbol, tissue, expn) %>%
        ggplot(aes(symbol, tissue)) +
        geom_tile(aes(fill = expn), colour = "white") +
        scale_fill_gradient(low = "white", high = "steelblue") +
        labs(x='', y = '') +
        theme_bw() +
        ggtitle(plot_title) +
        guides(fill = guide_legend(title = "expression",
                                   title.hjust = 0.5)) +
        theme(plot.title = element_text(hjust = 0.5,size = 40),
              axis.title = element_text(size=15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 24),
              axis.text.y = element_text(size = 18),
              legend.key.size = unit(60, 'points'),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 18),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())

    png(plotfile, width = 1920, height = 960)
    print(p)
    dev.off()

    dataset$resources$aggregate <- list(filename = plotfile,
                                        dataset_name = dataset$name,
                                        resource_name = "GTEx_tissue_expn_heatmap",
                                        type = "png",
                                        plot = p,
                                        published = FALSE,
                                        cleared = FALSE)

    gene_annotation$datasets[[dataset$name]] <- dataset

    return(gene_annotation)
}

#' Generate gene expression annotations
#'
#' @description Generate gene annotation artifacts for tissue-specific expression
#'
#' @param gene_annotation A list as output by 'annotate_genes'
#'
#' @return A list as output by 'annotate_genes', with gene expression annotation
#' artifact locations added
#' @export
annotate_gene_tissue_hpa_specificity <- function(gene_annotation) {

    dataset <- list(name = tolower(paste("tissue_hpa_expn",
                                         stringr::str_replace_all(gene_annotation$metadata$setname,
                                                                  "[^[a-z][A-Z][0-9][-_]]",
                                                                  "_"), sep ="_")),
                    description = paste("Tissue-specific protein expression of",
                                        ifelse(gene_annotation$metadata$setname=="all",
                                               "all genes",
                                               paste("genes from",
                                                     gene_annotation$metadata$setname))),
                    published = FALSE,
                    cleared = FALSE)
    ## HPA
    query_string <- paste("select *",
                          "from hpa.protein_detection where",
                          "gene_id in",
                          gene_annotation$metadata$gene_string)
    #res <- noctua::dbExecute(gene_annotation$metadata$con, query_string)
    res <- mazer::athena_query(query_string, gene_annotation$metadata$con)

    hpa <- res %>%
        tibble::as_tibble() %>%
        dplyr::mutate(protein = as.numeric(factor(protein,
                                       levels=c("Not detected",
                                                "Low",
                                                "Medium",
                                                "High"))) - 1)

    ## heatmap of HPA protein detection by tissue
    plotfile <- fs::path(gene_annotation$metadata$out_dir,
                         stringr::str_replace_all(paste0(gene_annotation$metadata$setname,
                                                         "_HPA_tissue_protein_heatmap.png"),
                                                  " ", "_"))

    plot_title <- ifelse(gene_annotation$metadata$setname=="all",
                         "Tissue protein detection of all genes",
                         paste("Tissue protein detection for genes from",
                               gene_annotation$metadata$setname))

    p <- hpa %>%
        dplyr::group_by(gene_id, symbol, tissue) %>%
        dplyr::summarize(max_detection = max(protein)) %>%
        dplyr::select(symbol, tissue, max_detection) %>%
        ggplot(aes(symbol, tissue)) +
        geom_tile(aes(fill = max_detection), colour = "white") +
        scale_fill_gradient(low = "white", high = "steelblue")+
        labs(x='', y = '') +
        theme_bw() +
        ggtitle(plot_title) +
        guides(fill = guide_legend(title = "Protein level",
                                   title.hjust = 0.5)) +
        #theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5,size = 40),
              axis.title = element_text(size=15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 24),
              axis.text.y = element_text(size = 18),
              legend.key.size = unit(60, 'points'),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 18),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())

    png(plotfile, width = 1920, height = 960)
    print(p)
    dev.off()

    dataset$resources$aggregate <- list(filename = plotfile,
                                        dataset_name = dataset$name,
                                        resource_name = "HPA_tissue_protein_heatmap",
                                        type = "png",
                                        published = FALSE,
                                        cleared = FALSE)

    gene_annotation$datasets[[dataset$name]] <- dataset

    return(gene_annotation)
}

#' Generate gene expression annotations
#'
#' @description Generate gene annotation artifacts for tissue-specific expression
#'
#' @param gene_annotation A list as output by 'annotate_genes'
#'
#' @return A list as output by 'annotate_genes', with gene expression annotation
#' artifact locations added
#' @export
annotate_gene_tissue_celltype_hpa_specificity <- function(gene_annotation) {

    dataset <- list(name = tolower(paste("tissue_celltype_hpa_expn",
                                         stringr::str_replace_all(gene_annotation$metadata$setname,
                                                                  "[^[a-z][A-Z][0-9][-_]]",
                                                                  "_"), sep ="_")),
                    description = paste("Tissue-specific celltype protein expression of",
                                        ifelse(gene_annotation$metadata$setname=="all",
                                               "all genes",
                                               paste("genes from",
                                                     gene_annotation$metadata$setname))),
                    published = FALSE,
                    cleared = FALSE)
    ## HPA
    query_string <- paste("select *",
                          "from hpa.protein_detection where",
                          "gene_id in",
                          gene_annotation$metadata$gene_string)
    #res <- noctua::dbExecute(gene_annotation$metadata$con, query_string)
    res <- mazer::athena_query(query_string, gene_annotation$metadata$con)

    hpa <- res %>%
        tibble::as_tibble() %>%
        dplyr::mutate(protein = as.numeric(factor(protein,
                                                  levels=c("Not detected",
                                                           "Low",
                                                           "Medium",
                                                           "High"))) - 1)

    ## heatmap of HPA protein detection within tissue
    plotfile <- fs::path(gene_annotation$metadata$out_dir,
                         stringr::str_replace_all(paste(gene_annotation$metadata$setname,
                                                        filter_tissue,
                                                        "HPA_tissue_protein_heatmap.png"),
                                                  " ", "_"))

    plot_title <- ifelse(gene_annotation$metadata$setname=="all",
                         "Tissue protein detection of all genes",
                         paste("Protein detection in", filter_tissue,
                               "for genes from", gene_annotation$metadata$setname))

    p <- hpa %>%
        dplyr::filter(tissue == filter_tissue) %>%
        dplyr::select(symbol, cell_type, protein) %>%
        ggplot(aes(symbol, cell_type)) +
        geom_tile(aes(fill = protein), colour = "white") +
        scale_fill_gradient(low = "white", high = "steelblue")+
        labs(x = '', y = '') +
        theme_bw() +
        ggtitle(plot_title) +
        guides(fill = guide_legend(title = "Protein level",
                               title.hjust = 0.5)) +
        theme(plot.title = element_text(hjust = 0.5,size = 40),
              axis.title = element_text(size=15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 24),
              axis.text.y = element_text(size = 18),
              legend.key.size = unit(60, 'points'),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 18),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())


    png(plotfile, width = 1920, height = 960)
    print(p)
    dev.off()

    dataset$resources$aggregate <- list(filename = plotfile,
                                        dataset_name = dataset$name,
                                        resource_name = "HPA_tissue_specific_protein_heatmap",
                                        type = "png",
                                        published = FALSE,
                                        cleared = FALSE)

    gene_annotation$datasets[[dataset$name]] <- dataset

    return(gene_annotation)
}

#' Generate gene expression annotations
#'
#' @description Build KPMP gene expression annotation artifacts
#'
#' @param gene_annotation A list as output by 'annotate_genes'
#'
#' @return A list as output by 'annotate_genes', with gene expression annotation
#' artifact locations added
#' @export
annotate_gene_cell_type_specificity_kpmp <- function(gene_annotation) {


    dataset <- list(name = tolower(paste("kpmp_kidney_celltype_rna_expn",
                                         stringr::str_replace_all(gene_annotation$metadata$setname,
                                                                  " ",
                                                                  "_"), sep ="_")),
                    description = paste("KPMP cell-type-specific RNA expression of",
                                                ifelse(gene_annotation$metadata$setname=="all",
                                                       "all genes",
                                                       paste("genes from",
                                                             gene_annotation$metadata$setname))),
                    published = FALSE,
                    cleared = FALSE)

    # KPMP single cell pseudobulk reads
    query_string <- paste("select * from kpmp.pseudobulk",
                          "where source ='pseudobulk-DGEList-cells'",
                          "and gene_id in",
                          gene_annotation$metadata$gene_string)
    #res <- noctua::dbExecute(gene_annotation$metadata$con, query_string)
    res <- mazer::athena_query(query_string, gene_annotation$metadata$con)

    cpm <- res

    kpmp_pseudobulk <- cpm %>%
        tibble::as_tibble() %>%
        dplyr::mutate(cpm = log(cpm + 5), # prior_count
                      structure = factor(tolower(structure),
                                         levels = c("vessels",
                                                    "interstitium",
                                                    "proximal tubules",
                                                    "intermediate tubules",
                                                    "distal tubules",
                                                    "collecting tubules",
                                                    "renal corpuscle"))) %>%
        dplyr::left_join(gene_annotation$metadata$ensembl %>% distinct(),
                         by = "gene_id") %>%
        dplyr::select(gene_id, symbol, sample, cpm, cell_abbrev,
                      cell_label, condition.l2, structure)

    cell_type_levels <- kpmp_pseudobulk %>%
        dplyr::group_by(structure) %>%
        dplyr::mutate(structure_cell_type = paste(structure, cell_abbrev, sep="_")) %>%
        arrange(cell_abbrev, .by_group=TRUE) %>%
        ungroup() %>%
        dplyr::select(structure_cell_type) %>%
        distinct() %>%
        tidyr::separate_wider_delim(structure_cell_type, delim="_", names = c("structure", "cell_type")) %>%
        pull(cell_type)


    kpmp_pseudobulk <- kpmp_pseudobulk %>%
        dplyr::mutate(cell_type = factor(cell_abbrev, levels = cell_type_levels))

    ## boxplots of gene expression by cell type and structure

    dataset$resources <- setNames(lapply(names(gene_annotation$metadata$genes),
                                         function(focal_gene_id) {
        gene <- list()

        focal_symbol <- gene_annotation$metadata$genes[[focal_gene_id]]

        plotfile <- fs::path(gene_annotation$metadata$out_dir,
                             paste0(focal_symbol, "_KPMP_celltype_expn_boxplot.png"))

        p <- kpmp_boxplot(kpmp_pseudobulk, focal_gene_id, focal_symbol, tissue)

        png(plotfile, width = 1920, height = 960)
        print(p)
        dev.off()

        gene <- list(filename = plotfile,
                     dataset_name = dataset$name,
                     resource_name = paste(focal_symbol, "kpmp_celltype_expn_boxplot", sep="_"),
                     type = "png",
                     published = FALSE,
                     cleared = FALSE)
        return(gene)
    }),
    gene_annotation$metadata$genes)

    ## heatmap of expression in kidney cell-types from KPMP
    median_expn <- kpmp_pseudobulk %>%
        left_join(kpmp_pseudobulk %>%
                      group_by(gene_id) %>%
                      summarize(median_cpm = median(cpm),
                                sd_cpm = sd(cpm)),
                  by="gene_id") %>%
        dplyr::mutate(expn = (cpm - median_cpm)/sd_cpm)

    median_expn %>% readr::write_csv(fs::path(out_dir))

    # # Example use of `median_expn`
    # median_expn %>% readr::write_csv(fs::path(out_dir, "kpmp_expn_renal_new_targets", ext="csv"))
    # ckanr::package_search("kpmp_kidney_celltype_rna_expn_renal_new_targets_2023")
    # # "2912fa34-a3b8-41a5-aedf-5198ab627888"
    # ckanr::resource_create(package_id = "2912fa34-a3b8-41a5-aedf-5198ab627888",
    #                        name="kpmp_data_renal_new_targets",
    #                        description="KPMP expression data for renal new targets",
    #                        format="CSV", upload=fs::path(out_dir, "kpmp_expn_renal_new_targets", ext="csv"))

    plotfile <- fs::path(gene_annotation$metadata$out_dir,
                         stringr::str_replace_all(paste0(gene_annotation$metadata$setname,
                                                         "_KPMP_kidney_celltype_expn_heatmap.png"),
                                                  " ", "_"))

    title <- paste("Cell-type specific tissue expression of",
                   ifelse(gene_annotation$metadata$setname=="all",
                     "all genes",
                    paste("genes from",
                          gene_annotation$metadata$setname)))
    p <- kpmp_pseudobulk %>%
        dplyr::rename(expn = cpm) %>%
        dplyr::select(symbol, cell_type, structure, expn) %>%
        ggplot(aes(x = symbol, y = cell_type, fill = expn)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "steelblue") +
        labs(x='', y = '') +
        theme_bw() +
        ggtitle(title) +
        guides(fill = guide_legend(title = "log(CPM)",
                                   title.hjust = 0.5)) +
        theme(plot.title = element_text(hjust = 0.5,size = 40),
              axis.title = element_text(size=15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 24),
              axis.text.y = element_text(size = 18),
              legend.key.size = unit(60, 'points'),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 18),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())


    png(plotfile, width = 1920, height = 960)
    print(p)
    dev.off()

    dataset$resources$aggregate <- list(filename = plotfile,
                                        dataset_name = dataset$name,
                                        resource_name = "KPMP_kidney_celltype_expn_heatmap",
                                        type = "png",
                                        published = FALSE,
                                        cleared = FALSE)

    gene_annotation$datasets[[dataset$name]] <- dataset

    return(gene_annotation)
}

#' Generate MAF tables
#'
#' @description Generate gene annotation artifacts for tissue-specific expression
#'
#' @param gene_annotation A list as output by 'annotate_genes'
#'
#' @return A list as output by 'annotate_genes', with gene expression annotation
#' artifact locations added
#' @export
annotate_gene_maf <- function(gene_annotation) {

    dataset <- list(name = tolower(paste("gene_maf_",
                                         stringr::str_replace_all(gene_annotation$metadata$setname,
                                                                  "[^[a-z][A-Z][0-9][-_]]",
                                                                  "_"), sep ="_")),
                    description = paste("Minor Allele Frequency for variants in",
                                        ifelse(gene_annotation$metadata$setname=="all",
                                               "all genes",
                                               paste("genes from",
                                                     gene_annotation$metadata$setname))),
                    published = FALSE,
                    cleared = FALSE)

    query_string <- paste("select chr_id, position, ref_allele, alt_allele, gene_id_any,",
                          "af.gnomad_afr as African, af.gnomad_amr as American,",
                          "af.gnomad_asj as Ashkenazi, af.gnomad_eas as East_Asian,",
                          "af.gnomad_fin as Finnish, af.gnomad_nfe as Non_Finnish_European,",
                          "af.gnomad_nfe_est as Estonian, af.gnomad_nfe_nwe as NW_European,",
                          "af.gnomad_nfe_onf as Other_Non_Finnish, af.gnomad_nfe_seu as S_European,",
                          "af.gnomad_oth as Other",
                          "from open_targets_genetics.variant_index where gene_id_any in",
                          gene_annotation$metadata$gene_string)

    #res <- noctua::dbExecute(gene_annotation$metadata$con, query_string)
    res <- mazer::athena_query(query_string, gene_annotation$metadata$con)

    maf_output <- noctua::dbFetch(res) %>%
        tibble::as_tibble() %>%
        dplyr::rename(gene_id = gene_id_any,
                      chr = chr_id, pos = position, ref = ref_allele, alt = alt_allele) %>%
        tidyr::unite("label", chr:alt) %>%
        dplyr::select(-gene_id)
    noctua::dbClearResult(res)


    internal_variants %>%
        left_join(maf_output, by="label")

    dataset$resources$maf_table <- list(filename = plotfile,
                                        dataset_name = dataset$name,
                                        resource_name = "maf_table",
                                        type = "csv",
                                        published = FALSE,
                                        cleared = FALSE)

    gene_annotation$datasets[[dataset$name]] <- dataset

    return(gene_annotation)
}


#' Generate HTML table of mouse phenotypes by gene
#'
#' @description Build gene annotation artifacts for mouse genes
#'
#' @param gene The human gene of interest
#'
#' @return A formatted HTML table of phenotypes for the mouse ortholog
#' @export
annotate_mouse_gene_phenotypes <- function(gene_annotation) {

    dataset <- list(name = tolower(paste("mouse_phenotypes",
                                         stringr::str_replace_all(gene_annotation$metadata$setname,
                                                                  "[^[a-z][A-Z][0-9][-_]]",
                                                                  "_"), sep ="_")),
                    description = paste("Mouse phenotypes of",
                                        ifelse(gene_annotation$metadata$setname=="all",
                                               "all genes",
                                               paste("genes from",
                                                     gene_annotation$metadata$setname))),
                    published = FALSE,
                    cleared = FALSE)

    textfile <- fs::path(gene_annotation$metadata$out_dir,
                         "mouse_phenotypes.csv")
    htmlfile <- fs::path(gene_annotation$metadata$out_dir,
                         "mouse_phenotypes.html")

    # query_string <- paste("select split_part(gene_id, '.', 1) as gene_id, symbol, sample, value",
    #                       "from gtex.v8 where measurement_type='gene_tpm'",
    #                       "and split_part(gene_id, '.', 1) in",
    #                       gene_annotation$metadata$gene_string)


    # TODO: replace with athena query to OpenTargets
  phenotype_table <- readRDS(fs::path(data_dir, "opentargets/mouse_models.rds")) %>% # need to regenerate after data update
      dplyr::filter(targetFromSourceId %in% names(gene_annotation$metadata$genes)) %>%
      dplyr::rename(phenotype_id = modelPhenotypeId,
                    phenotype = modelPhenotypeLabel,
                    gene_id = targetFromSourceId,
                    mouse_gene_id = targetInModelEnsemblId,
                    mgi_id = targetInModelMgiId,
                    mouse_symbol = targetInModel) %>%
      dplyr::left_join(tibble::enframe(gene_annotation$metadata$genes,
                                       name="gene_id",
                                       value = "symbol"),
                       by="gene_id") %>%
      unnest(modelPhenotypeClasses) %>%
      dplyr::rename(phenotype_class_id = id,
                    phenotype_class = label) %>%
      unnest(biologicalModels) %>%
      dplyr::rename(model_id = id) %>%
      dplyr::select(symbol, gene_id, mouse_symbol, mouse_gene_id,
                    phenotype_id, phenotype, phenotype_class_id, phenotype_class,
                    model_id, allelicComposition, geneticBackground,
                    ) %>%
      readr::write_csv(textfile)

  phenotype_presentation_table <- phenotype_table %>%
      flextable::regulartable() %>%
      flextable::theme_box() %>%
      flextable::merge_v(c("symbol", "gene_id", "mouse_symbol", "mouse_gene_id",
                           "phenotype_id", "phenotype",
                           "phenotype_class_id", "phenotype_class")) %>%
      flextable::valign(valign = "top", part = "all") %>%
      flextable::set_table_properties(width = 1, layout = "autofit") %>%
      flextable::set_header_labels(gene_id = "Ensembl ID",
                                   mouse_symbol = "mouse gene",
                                   mgi_id = "mouse gene MGI",
                                   mouse_gene_id = "mouse Ensembl ID",
                                   phenotype_id = "phenotype",
                                   model_id = "model ID",
                                   allelicComposition = "allelic composition",
                                   geneticBackground = "genetic background") %>%
      save_as_html(values = NULL,
                   path = htmlfile,
                   lang = "en",
                   title = dataset$description)


  dataset$resources$text <- list(filename = textfile,
                                      dataset_name = dataset$name,
                                      resource_name = "mouse_phenotypes_text",
                                      type = "csv",
                                      published = FALSE,
                                      cleared = FALSE)

  dataset$resources$html <- list(filename = htmlfile,
                                 dataset_name = dataset$name,
                                 resource_name = "mouse_phenotypes_html",
                                 type = "html",
                                 published = FALSE,
                                 cleared = FALSE)

  gene_annotation$datasets[[dataset$name]] <- dataset

  return(gene_annotation)

  # phenotype_html_table <- phenotype_presentation_table %>%
  #     htmltools_value()
}

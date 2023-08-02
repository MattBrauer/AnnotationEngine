#' Annotate set of variants
#'
#' @description Set up initial variant annotation data structure
#'
#' @param variantlist A tibble of variants from CKAN
#' @param setname A character string (default "all") describing the gene grouping
#' @param tissue A character string (default "all") specifying the focal tissue
#' @param data_dir The location where output will be placed
#' @param edb An `edb` object, specific to species
#'
#' @return A list containing the annotation artifacts created, as files, in
#'  `data_dir`
#' @export
annotate_variants <- function(variantlist,
                           setname = "all",
                           tissue = "all",
                           data_dir,
                           edb = edb_Hs) {

    genes_edb <- suppressWarnings(genes(edb))

    warnings <- list()

    data_dir_check <- geneset_check <- FALSE

    if(fs::dir_exists(data_dir)) {
        data_dir_check <- TRUE
    } else {
        warning("data_dir does not exist")
    }

    is_valid <- variantlist$gene_id %in% genes_edb$gene_id

    if(all(is_valid)) {
        geneset_check <- TRUE
    } else {
        warning("Not all variants have been mapped to valid gene_ids")
        warnings[["ensembl_annotation_failures"]] <- variantlist$label[!is_valid]
    }

    # TODO: V2G for unmapped variants

    annotation_granges <- genes_edb[genes_edb$gene_id %in% variantlist$gene_id[is_valid]]
    annotation <- mcols(annotation_granges) %>%
        tibble::as_tibble() %>%
        dplyr::select(gene_id, symbol, description,
                      gene_id_version, gene_biotype, canonical_transcript)

    symbols <- annotation %>%
        dplyr::select(gene_id, symbol) %>%
        distinct() %>%
        tibble::deframe()

    variants <- variantlist %>%
        dplyr::select(variant_id, gene_id) %>%
        distinct() %>%
        tibble::deframe()

    #TODO: map variants to canonical proteins to add lolliplot or dandelion plot
    variant_annotation <- list(
        metadata = list(ensembl = variantlist %>% left_join(annotation, by="gene_id"),
                        tissue = tissue,
                        setname = setname,
                        data_dir = data_dir,
                        out_dir = out_dir,
                        edb = edb,
                        genes = symbols,
                        variants = variants,
                        variant_granges = setNames(GRanges(seqnames=paste0("chr", variantlist$chr),
                                                  ranges=IRanges(start=variantlist$pos, width=1),
                                                  mcols=variantlist),
                                                  variantlist$label),
                        gene_string = paste0('(',
                                             names(symbols) %>%
                                                 stringr::str_c("'", ., "'") %>%
                                                 stringr::str_flatten_comma(),
                                             ')'),
                        variant_string = paste0('(',
                                                names(variants) %>%
                                                    stringr::str_c("'", ., "'") %>%
                                                    stringr::str_flatten_comma(),
                                                ')'),
                        con = noctua::dbConnect(noctua::athena(),
                                                profile_name = user_profile,
                                                s3_staging_dir = s3_staging_dir)),
        warnings = warnings,
        datasets = list()
    )

    return(variant_annotation)
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
annotate_variant_maf <- function(variant_annotation) {

    dataset <- list(name = tolower(paste("variant_maf_",
                                         stringr::str_replace_all(variant_annotation$metadata$setname,
                                                                  "[^[a-z][A-Z][0-9][-_]]",
                                                                  "_"), sep ="_")),
                    description = paste("Minor Allele Frequency for internal variants in",
                                        ifelse(variant_annotation$metadata$setname=="all",
                                               "all genes",
                                               paste("genes from",
                                                     variant_annotation$metadata$setname))),
                    published = FALSE,
                    cleared = FALSE)

    query_string <- paste("select concat(chr_id, '_', cast(position as varchar), '_', ref_allele, '_', alt_allele) as variant_id,",
                          "chr_id, position, ref_allele, alt_allele, gene_id_any,",
                          "af.gnomad_afr as African, af.gnomad_amr as American,",
                          "af.gnomad_asj as Ashkenazi, af.gnomad_eas as East_Asian,",
                          "af.gnomad_fin as Finnish, af.gnomad_nfe as Non_Finnish_European,",
                          "af.gnomad_nfe_est as Estonian, af.gnomad_nfe_nwe as NW_European,",
                          "af.gnomad_nfe_onf as Other_Non_Finnish, af.gnomad_nfe_seu as S_European,",
                          "af.gnomad_oth as Other",
                          "from open_targets_genetics.variant_index where",
                          "concat(chr_id, '_', cast(position as varchar), '_', ref_allele, '_', alt_allele) in",
                          variant_annotation$metadata$variant_string)

    res <- noctua::dbExecute(variant_annotation$metadata$con, query_string)
    maf_output <- noctua::dbFetch(res) %>%
        tibble::as_tibble() %>%
        dplyr::rename(chr = chr_id, pos = position, ref = ref_allele, alt = alt_allele) %>%
        left_join(tibble::enframe(variant_annotation$metadata$variants, name="variant_id", value="gene_id"), by="variant_id") %>%
        left_join(tibble::enframe(variant_annotation$metadata$genes, name="gene_id", value="symbol"), by="gene_id") %>%
        group_by(symbol)
    noctua::dbClearResult(res)

    maf_output_list <- setNames(maf_output %>% group_split(),
                                maf_output %>% group_keys() %>% pull(symbol))


    dataset$resources$maf_table <- list(filename = plotfile,
                                        dataset_name = dataset$name,
                                        resource_name = "maf_table",
                                        type = "csv",
                                        published = FALSE,
                                        cleared = FALSE)

    gene_annotation$datasets[[dataset$name]] <- dataset

    return(gene_annotation)
}


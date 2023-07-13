#' Prepare organization structure for publication
#'
#' @description Get CKAN entities for publication of annotation objects
#'
#' @param ckan_org_id The name or id of the organization to publish to
#' @param dataset_title The title or display name of the dataset to publish to
#' @param url The CKAN URL
#' @param key The CKAN key
#'
#' @return A "ckan_org" object, comprising a list of the names and ids of
#' datasets/packages in the organization
#' @export
prepare_for_publication <- function(ckan_org_id,
                                    url = ckan_url,
                                    key = ckan_key,
                                    dry_run = TRUE) {

    orgs <- try(ckanr::organization_show(id = ckan_org_id,
                                         include_datasets = TRUE,
                                         as = "table",
                                         url = ckan_url,
                                         key = ckan_key),
                silent = TRUE)

    # if org exists, get or create dataset id
    if(orgs$is_organization) {
        if(orgs$package_count > 0) {
            packages <- setNames(orgs$packages$id,
                                 orgs$packages$name)
        } else {
            packages <- NULL
        }
    } else {
        warning(paste("Need exactly one organization with the given identifier (found",
                      nrow(orgs), ")."))
        packages <- NULL
    }
    return(list(id = orgs$id,
                name = orgs$name,
                packages = packages,
                url = url,
                key = key))
}

#' If a package with the given name exists, return its ID, otherwise create it
#' and return its ID
#'
#' @description Get a package's ID or if not existing create and get ID
#'
#'
#' @export
get_or_create_package_id <- function(ckan_org_id,
                                     package_name,
                                     package_title,
                                     dry_run = dry_run) {

    packages <- ckanr::package_search(paste0("owner_org:",ckan_org_id))
    package_ids <- setNames(
        unlist(lapply(packages$results, function(pkg) {
            pkg$id
        })),
        unlist(lapply(packages$results, function(pkg) {
            pkg$name
        })))

    if(package_name %in% names(package_ids)) {
        pkg_id <- package_ids[package_name]
    } else {
        if(dry_run) {
            pkg_id <- "0000"
        } else {
            pkg <- ckanr::package_create(name = package_name,
                                         title = package_title,
                                         owner_org = ckan_org_id)
            pkg_id <- pkg$id
        }
    }
    return(pkg_id)
}

#' Look for groups with given names. If a group with the name exists, return its
#' ID, otherwise create it and return its ID
#'
#' @description Get a vector of group ID's, creating where missing
#'
#' @param group_names Vector of names of potential groups
#' @param ckan_org CKAN organization object
#'
#' @return A named vector of group ids
#'
#' @export
get_or_create_group_ids <- function(group_names,
                                    ckan_org,
                                    dry_run = TRUE) {

    group_names <- setNames(paste0(tolower(group_names), "_annotation"),
                            group_names)

    groups <- ckanr::group_list(limit = 1000,
                                url = ckan_org$url,
                                key = ckan_org$key,
                                as = "list")

    existing_groups <- unlist(lapply(groups, function(group) {
        if(group$name %in% group_names) {
            return(setNames(group$id, names(group_names)[group_names==group$name]))
        }
    }))

    needed_groups <- unlist(lapply(as.list(names(group_names)[!group_names %in% names(existing_groups)]),
                            function(group_name) {
                                if(dry_run) {
                                    res <- list(id = "new", name = group_names[group_name])
                                } else {
                                    res <- ckanr::group_create(name = unname(group_names[group_name]),
                                                               title = group_name,
                                                               url = ckan_org$url,
                                                               key = ckan_org$key,
                                                               as = "list")
                                }
                                return(setNames(res$id, group_name))
                            }))

    return(c(existing_groups, needed_groups))

}

#' Publish annotations from a set of genes to CKAN
#'
#' @description Publish artifacts from annotation object to CKAN
#'
#' @param annotation An annotation structure as returned by gene
#'  annotation or variant annotation functions
#' @param ckan_org A list with the organization id and list of packages
#' @param update Whether to overwrite (default) or fail if resource exists
#' @param dry_run If TRUE, output publication steps without executing
#'
#' @return The status of the publication
#' @export
publish_annotation <- function(annotation,
                                    ckan_org,
                                    update = TRUE,
                                    dry_run = FALSE,
                                    prefix = NULL) {


    annotation$datasets <- lapply(annotation$datasets,
                                       function(dataset) {
                                           if(!is.null(prefix)) dataset$description <- paste(prefix, dataset$description, sep="_")
                                           dataset$id <- get_or_create_package_id(ckan_org$id,
                                                                                  package_name = dataset$name,
                                                                                  package_title = dataset$description,
                                                                                  dry_run = dry_run)
                                           return(dataset)
    })
    annotation$datasets <- lapply(annotation$datasets,
                  function(dataset) {
                      print(dataset$name)
                      if(!dataset$published & !dataset$cleared) {
                          existing_resources <- ckanr::package_show(dataset$id,
                                                      as="table")$resources
                          current_resources <- dataset$resources
                          dataset$resources <- lapply(dataset$resources,
                                                      function(resource) {
                                                          if(!resource$published & !resource$cleared) {
                                                              # publish resource
                                                              if(resource$resource_name %in% existing_resources$name & update) {
                                                                  res_id <- existing_resources$id[existing_resources$name==resource$resource_name]
                                                                  res <- try(ckanr::resource_update(id = res_id,
                                                                                                path = resource$filename),
                                                                             silent = TRUE)
                                                              } else {
                                                                  res <- ckanr::resource_create(package_id = dataset$id,
                                                                                                name = resource$resource_name,
                                                                                                upload = resource$filename)
                                                              }
                                                              if(class(res) == "ckan_resource") {
                                                                resource$published <- TRUE
                                                              }
                                                          }
                                                          if(!resource$cleared) {
                                                              # clear resource
                                                              unlink(resource$filename)
                                                              resource$cleared <- !fs::file_exists(resource$filename)
                                                              resource$filename <- NULL
                                                          }
                                                          return(resource)
                                                      })
                          dataset$published <- all(unlist(lapply(dataset$resources, function(resource) resource$published)))
                          dataset$cleared <- all(unlist(lapply(dataset$resources, function(resource) resource$cleared)))
                      }
                      return(dataset)
    })
    return(annotation)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/artifact_publication_functions.R
\name{prepare_for_publication}
\alias{prepare_for_publication}
\title{Prepare organization structure for publication}
\usage{
prepare_for_publication(
  ckan_org_id,
  url = ckan_url,
  key = ckan_key,
  dry_run = TRUE
)
}
\arguments{
\item{ckan_org_id}{The name or id of the organization to publish to}

\item{url}{The CKAN URL}

\item{key}{The CKAN key}

\item{dataset_title}{The title or display name of the dataset to publish to}
}
\value{
A "ckan_org" object, comprising a list of the names and ids of
datasets/packages in the organization
}
\description{
Get CKAN entities for publication of annotation objects
}

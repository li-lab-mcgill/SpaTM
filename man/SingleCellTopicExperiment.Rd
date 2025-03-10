% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellTopicExperiment.R
\name{SingleCellTopicExperiment}
\alias{SingleCellTopicExperiment}
\title{Constructor for SingleCellTopicExperiment}
\usage{
SingleCellTopicExperiment(
  sce,
  guided = FALSE,
  labels = NULL,
  K = NULL,
  hvg = NULL,
  verbal = FALSE
)
}
\arguments{
\item{sce}{A `SingleCellExperiment` object containing single-cell RNA-seq data.}

\item{guided}{Logical. If `TRUE`, topic modeling is guided by metadata labels.}

\item{labels}{A character string specifying the column name in `colData(sce)` that contains metadata labels.}

\item{K}{Integer. The number of topics to use. Defaults to 10 if not provided.}

\item{hvg}{Integer or NULL. If numeric, filters for the top `hvg` highly variable genes.}

\item{verbal}{Logical. If `TRUE`, prints progress messages.}
}
\value{
A `SingleCellTopicExperiment` object with initialized topic modeling components.
}
\description{
Initializes a `SingleCellTopicExperiment` object from a `SingleCellExperiment` object,
optionally filtering for highly variable genes (HVG) and setting topic modeling priors.
}

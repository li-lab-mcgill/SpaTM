# Spatial Topic Experiment


#### Extending spe object
#' SpatialTopicExperiment Class and Constructor
#'
#' The `SpatialTopicExperiment` class extends `SpatialExperiment`, incorporating topic modeling
#' for spatial transcriptomics data. It follows the same logic as `SingleCellTopicExperiment` (@seealso [SingleCellTopicExperiment()]),
#' adapting it for spatial data.
#'
#' @slot theta A matrix storing cell-topic distributions.
#' @slot phi A matrix storing gene-topic distributions.
#' @slot alphaPrior A matrix storing topic priors per cell.
#' @slot betaPrior A matrix storing gene priors per topic.
#' @slot nwk A matrix storing topic-gene counts.
#' @slot ndk A matrix storing topic-cell counts.
#'
#' @export
setClass(
  "SpatialTopicExperiment",
  contains="SpatialExperiment",
  slots=c(theta="matrix",
          phi = "matrix",
          alphaPrior = "matrix",
          betaPrior = "matrix",
          nwk = "matrix",
          ndk = "matrix")
)

##### Constructor/Prep Function ####
#' Constructor for SpatialTopicExperiment
#'
#' Initializes a `SpatialTopicExperiment` object, setting up topic modeling structures
#' for spatial transcriptomics data. It follows a similar process as `SingleCellTopicExperiment`.
#'
#' @param spe A `SpatialExperiment` object.
#' @param guided Logical; if `TRUE`, the initialization is guided by cell labels.
#' @param labels A character string specifying the column in `colData(spe)` containing cell labels.
#' @param K Number of topics. If `NULL` and `guided` is `FALSE`, defaults to `K = 10`.
#' @param hvg Integer specifying the number of highly variable genes to retain.
#' @param verbal Logical; if `TRUE`, prints progress messages.
#'
#' @return A `SpatialTopicExperiment` object initialized with topic modeling parameters.
#'
#' @export
#' @import scuttle SpatialExperiment
#' @importFrom scran modelGeneVar getTopHVGs
#' @importClassesFrom Matrix dgCMatrix
SpatialTopicExperiment <- function(spe,guided = FALSE,
                                      labels = NULL,
                                      K = NULL,
                                      hvg = NULL,
                                      verbal = FALSE,
                                      balanced = FALSE){
  if (!is(counts(spe),'dgCMatrix')){
    message('Converting counts to sparse matrix (dgCmatrix). \nSpaTM is currently not compatible with other matrix formats.')
    counts(spe) <- as(counts(spe),'dgCMatrix')
  }
  if (!guided & is.null(K)){
    message('No guided or topic number provided. Setting topics to K = 10. \nIf this is an error please re-run the function and assign a number of topics or a guide variable.\n')
    K <- 10
  }
  if (is.numeric(hvg)){
   if(verbal){ cat(paste("Filtering Top ",hvg," Highly variable genes\n",sep = ''))}
    ## Calculate HVG and filter
    if (!(any(names(assays(spe)) == 'logcounts'))){
      if(verbal){cat('Calculating log-normalized counts for HVG identification\n')}
      spe <- logNormCounts(spe)
    }
    dec <- modelGeneVar(spe)
    top_hvgs <- getTopHVGs(dec,n = hvg)
    spe <- spe[top_hvgs,]
  }
  spte <- new("SpatialTopicExperiment",spe)
  if (guided){
    if(!is.null(labels)){
      if (!is.factor(colData(spte)[,labels])){
        warning("Label needs to be a factor. Converting to factor.\n")
        colData(spte)[,labels] <- as.factor(colData(spte)[,labels])
      }
      spte$int_celltype <- as.numeric(colData(spte)[,labels])
      if (balanced){
        if(verbal){message("Adjusting alpha priors to accomodate for class imbalance")}
      }
      alphaPrior(spte) <- build_alpha(spte,labels,balanced = balanced)
      K <- ncol(alphaPrior(spte))
    } else {
      stop("labels argument must be defined if guided parameter is set to TRUE")
    }
  } else {
    alphaPrior(spte) <- build_alpha(spte,labels,K = K,balanced = balanced)
  }
  betaPrior(spte) <- build_beta(spte,K)

 if(verbal){ cat("Initialized Prior Matrices\n")}

  ndk(spte) <- matrix(0,ncol(spte),K)
  nwk(spte) <- matrix(0,nrow(spte),K)
  theta(spte) <- matrix(0,ncol(spte),K)
  phi(spte) <- matrix(0,nrow(spte),K)
  if(verbal){cat("Initialized Sufficient Statistics\n")}

  ## Assuming the labels are "Celltype"
  rownames(theta(spte)) <- rownames(ndk(spte)) <- rownames(alphaPrior(spte)) <- colnames(spte)
  rownames(phi(spte)) <- rownames(nwk(spte)) <- rownames(betaPrior(spte)) <- rownames(spte)

  colnames(ndk(spte)) <- colnames(nwk(spte)) <- colnames(betaPrior(spte)) <-
    colnames(alphaPrior(spte))
  colnames(theta(spte)) <- colnames(phi(spte)) <- colnames(alphaPrior(spte))
  ## Assuming the labels are "Celltype"

  spte$int_cell <- 1:ncol(spte)
  rowData(spte)$gene_ints <- 1:nrow(spte)
  return(spte)
}

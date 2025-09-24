## SingleCellTopicExperiment

# TODO finish documenting arguments


#### Extending SCE object
#' SingleCellTopicExperiment Class
#'
#' A class that extends `SingleCellExperiment` to incorporate topic modeling
#' components for single-cell RNA-seq analysis.
#'
#' @slot theta A matrix representing cell-topic distributions.
#' @slot phi A matrix representing gene-topic distributions.
#' @slot alphaPrior A matrix representing prior information for cell-topic distributions.
#' @slot betaPrior A matrix representing prior information for gene-topic distributions.
#' @slot nwk A matrix storing gene-topic counts.
#' @slot ndk A matrix storing cell-topic counts.
#'
#' @importFrom methods setClass
#' @export
setClass(
  "SingleCellTopicExperiment",
  contains="SingleCellExperiment",
  slots=c(theta="matrix",
          phi = "matrix",
          alphaPrior = "matrix",
          betaPrior = "matrix",
          nwk = "matrix",
          ndk = "matrix")
)

##### Constructor/Prep Function ####
#' Constructor for SingleCellTopicExperiment
#'
#' Initializes a `SingleCellTopicExperiment` object from a `SingleCellExperiment` object,
#' optionally filtering for highly variable genes (HVG) and setting topic modeling priors.
#'
#' @param sce A `SingleCellExperiment` object containing single-cell RNA-seq data.
#' @param guided Logical. If `TRUE`, topic modeling is guided by metadata labels.
#' @param labels A character string specifying the column name in `colData(sce)` that contains metadata labels.
#' @param K Integer. The number of topics to use. Defaults to 10 if not provided.
#' @param hvg Integer or NULL. If numeric, filters for the top `hvg` highly variable genes.
#' @param verbal Logical. If `TRUE`, prints progress messages.
#' @param balanced Logical. If `TRUE`, adapts guided priors to be 1/label proportion to address class imbalances. Only applied if guided = TRUE. default is FALSE.
#'
#' @return A `SingleCellTopicExperiment` object with initialized topic modeling components.
#'
#' @import scuttle SingleCellExperiment
#' @importFrom scran modelGeneVar getTopHVGs
#' @importClassesFrom Matrix dgCMatrix
#' @export
SingleCellTopicExperiment <- function(sce,guided = FALSE,
                                      labels = NULL,
                                      K = NULL,
                                      hvg = NULL,
                                      verbal = FALSE,
                                      balanced = FALSE,
                                      disease = NULL,
                                      disease_type = NULL){
  if (!is(counts(sce),'dgCMatrix')){
    message('Converting counts to sparse matrix (dgCmatrix). \nSpaTM is currently not compatible with other matrix formats.')
    counts(sce) <- as(counts(sce),'dgCMatrix')
  }
  if (!guided & is.null(K)){
    message('No guided or topic number provided. Setting topics to K = 10. \nIf this is an error please re-run the function and assign a number of topics or a guide variable.\n')
    K <- 10
  }
  if (is.numeric(hvg)){
    if(verbal){cat(paste("Filtering Top ",hvg," Highly variable genes\n",sep = ''))}
    ## Calculate HVG and filter
    if (!(any(names(assays(sce)) == 'logcounts'))){
      if(verbal){cat('Calculating log-normalized counts for HVG identification\n')}
      sce <- logNormCounts(sce)
    }
    dec <- modelGeneVar(sce)
    top_hvgs <- getTopHVGs(dec,n = hvg)
    sce <- sce[top_hvgs,]
  }
  scte <- new("SingleCellTopicExperiment",sce)
  if (guided){
    if(!is.null(labels)){
      if (!is.factor(colData(scte)[,labels])){
        warning("Label needs to be a factor. Converting to factor.\n")
        colData(scte)[,labels] <- as.factor(colData(scte)[,labels])
      }
      scte$int_celltype <- as.numeric(colData(scte)[,labels])
      if (balanced){
        if(verbal){message("Adjusting alpha priors to accomodate for class imbalance")}
      }
      if(verbal){
        if(!is.null(disease) & !is.null(disease_type)){
          cat("Running disease-informed initialization\n")
        }
      }
      alphaPrior(scte) <- build_alpha(te = scte,
                                      labels = labels,
                                      K = NULL,
                                      balanced = balanced,
                                      disease = disease,
                                      disease_type = disease_type)
      K <- ncol(alphaPrior(scte))


    } else {
      stop("labels argument must be defined if guided parameter is set to TRUE")
    }
  } else {
    alphaPrior(scte) <- build_alpha(scte,labels,K)
  }
  betaPrior(scte) <- build_beta(scte,K)

  if(verbal){cat("Initialized Prior Matrices\n")}

  ndk(scte) <- matrix(0,ncol(scte),K)
  nwk(scte) <- matrix(0,nrow(scte),K)
  theta(scte) <- matrix(0,ncol(scte),K)
  phi(scte) <- matrix(0,nrow(scte),K)
  if(verbal){cat("Initialized Sufficient Statistics\n")}

  ## Assuming the labels are "Celltype"
  rownames(theta(scte)) <- rownames(ndk(scte)) <- rownames(alphaPrior(scte)) <- colnames(scte)
  rownames(phi(scte)) <- rownames(nwk(scte)) <- rownames(betaPrior(scte)) <- rownames(scte)

  colnames(ndk(scte)) <- colnames(nwk(scte)) <- colnames(betaPrior(scte)) <-
    colnames(alphaPrior(scte))
  colnames(theta(scte)) <- colnames(phi(scte)) <- colnames(alphaPrior(scte))

  scte$int_cell <- 1:ncol(scte)
  rowData(scte)$gene_ints <- 1:nrow(scte)
  return(scte)
}

## SingleCellTopicExperiment




#### Extending SCE object
#' @rdname SingleCellTopicExperiment
#' @exportClass SingleCellTopicExperiment scGTM
#' @importFrom SingleCellExperiment SingleCellExperiment
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

SingleCellTopicExperiment <- function(sce,guided = FALSE,
                                      labels = NULL,
                                      K = NULL,
                                      hvg = NULL,
                                      verbal = FALSE){
  if (is.numeric(hvg)){
    if(verbal){cat(paste("Filtering Top ",hvg," Highly variable genes\n",sep = ''))}
    ## Calculate HVG and filter
    if (!(any(names(assays(sce)) == 'logcounts'))){
      if(verbal){cat('Calculating log-normalized counts for HVG identification\n')}
      spe <- logNormCounts(sce)
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
      alphaPrior(scte) <- build_alpha(scte,labels)
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
  
  if(verbal){cat("Initialized Sufficient Statistics\n")}
  
  ## Assuming the labels are "Celltype"
  rownames(ndk(scte)) <- rownames(alphaPrior(scte)) <- colnames(scte)
  rownames(nwk(scte)) <- rownames(betaPrior(scte)) <- rownames(scte)
  colnames(ndk(scte)) <- colnames(nwk(scte)) <- colnames(betaPrior(scte)) <-
    colnames(alphaPrior(scte))
  scte$int_cell <- 1:ncol(scte)
  rowData(scte)$gene_ints <- 1:nrow(scte)
  return(scte)
}


#### Print function ####



#Print Contents
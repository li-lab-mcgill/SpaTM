##Topic Experiment Utils

## Class Union to simplify utility functions
setClassUnion('TopicExperiment',c('SingleCellTopicExperiment','SpatialTopicExperiment'))

##Update object - Deprecated
update_scte <- function(scte,labels = NULL, guided = FALSE, K = 10){
  scte <- buildPriors(scte,labels,guided,K)
  scte <- build_SufStats(scte,ncol(alphaPrior(scte)))
  return(scte)
}

#### Getters & Setters ####

#### Alpha Prior ####

setGeneric("alphaPrior", function(te) standardGeneric("alphaPrior"))

setGeneric("alphaPrior<-", function(te, value) standardGeneric("alphaPrior<-"))

#' Get Alpha Prior Matrix
#'
#' Retrieves the alpha prior matrix (cell-by-topic) from a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @return The alpha prior matrix.
#' @export
setMethod(
  f = "alphaPrior",
  signature = c("TopicExperiment"),
  function(te) {
    return(te@alphaPrior)
  }
)

#' Set Alpha Prior Matrix
#'
#' Sets the alpha prior matrix (cell-by-topic) for a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @param value The new alpha prior matrix.
#' @return The updated TopicExperiment object.
#' @export
setReplaceMethod(
  f = "alphaPrior",
  signature = c("TopicExperiment"),
  function(te,value){
    te@alphaPrior <- value
    return(te)
  }
)

#### Beta Prior ####

setGeneric("betaPrior", function(te) standardGeneric("betaPrior"))

setGeneric("betaPrior<-", function(te, value) standardGeneric("betaPrior<-"))

#' Get Beta Prior Matrix
#'
#' Retrieves the beta prior matrix (gene-by-topic) from a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @return The beta prior matrix.
#' @export
setMethod(
  f = "betaPrior",
  signature = "TopicExperiment",
  function(te) {
    return(te@betaPrior)
  }
)

#' Set Beta Prior Matrix
#'
#' Sets the beta prior matrix (gene-by-topic) for a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @param value The new beta prior matrix.
#' @return The updated TopicExperiment object.
#' @export
setReplaceMethod(
  f = "betaPrior",
  signature = "TopicExperiment",
  function(te,value){
    te@betaPrior <- value
    return(te)
  }
)

## Cell by Topics (ndk) ####

setGeneric("ndk", function(te) standardGeneric("ndk"))

setGeneric("ndk<-", function(te, value) standardGeneric("ndk<-"))

#' Get Cell-Topic Matrix (ndk)
#'
#' Retrieves the cell-topic matrix (ndk) from a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @return The cell-topic matrix (ndk).
#' @export
setMethod(
  f = "ndk",
  signature = "TopicExperiment",
  function(te) {
    return(te@ndk)
  }
)

#' Set Cell-Topic Matrix (ndk)
#'
#' Sets the cell-topic matrix (ndk) for a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @param value The new cell-topic matrix (ndk).
#' @return The updated TopicExperiment object.
#' @export
setReplaceMethod(
  f = "ndk",
  signature = "TopicExperiment",
  function(te,value){
    te@ndk <- value
    return(te)
  }
)

## Gene by Topics (nwk) ####

setGeneric("nwk", function(te) standardGeneric("nwk"))

setGeneric("nwk<-", function(te, value) standardGeneric("nwk<-"))

#' Get Gene-Topic Matrix (nwk)
#'
#' Retrieves the gene-topic matrix (nwk) from a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @return The gene-topic matrix (nwk).
#' @export
setMethod(
  f = "nwk",
  signature = "TopicExperiment",
  function(te) {
    return(te@nwk)
  }
)

#' Set Gene-Topic Matrix (nwk)
#'
#' Sets the gene-topic matrix (nwk) for a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @param value The new gene-topic matrix (nwk).
#' @return The updated TopicExperiment object.
#' @export
setReplaceMethod(
  f = "nwk",
  signature = "TopicExperiment",
  function(te,value){
    te@nwk <- value
    return(te)
  }
)

## Theta ####

setGeneric("theta", function(te) standardGeneric("theta"))

setGeneric("theta<-", function(te, value) standardGeneric("theta<-"))

#' Get Theta Matrix (Cell-Topic Distributions)
#'
#' Retrieves the theta matrix (cell-topic distributions)
#' from a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @return The theta matrix.
#' @export
setMethod(
  f = "theta",
  signature = "TopicExperiment",
  function(te) {
    return(te@theta)
  }
)

#' Set Theta Matrix (Cell-Topic Distributions)
#'
#' Sets the theta matrix (cell-topic distributions) for a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @param value The new theta matrix.
#' @return The updated TopicExperiment object.
#' @export
setReplaceMethod(
  f = "theta",
  signature = "TopicExperiment",
  function(te,value){
    te@theta <- value
    return(te)
  }
)

## Phi ####

setGeneric("phi", function(te) standardGeneric("phi"))

setGeneric("phi<-", function(te, value) standardGeneric("phi<-"))

#' Get Phi Matrix (Gene-Topic Distributions)
#'
#' Retrieves the phi matrix (gene-topic distributions) from a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @return The phi matrix.
#' @export
setMethod(
  f = "phi",
  signature = "TopicExperiment",
  function(te) {
    return(te@phi)
  }
)

#' Set Phi Matrix (Gene-Topic Distributions)
#'
#' Sets the phi matrix (gene-topic distributions) for a TopicExperiment object.
#'
#' @param te A TopicExperiment object.
#' @param value The new phi matrix.
#' @return The updated TopicExperiment object.
#' @export
setReplaceMethod(
  f = "phi",
  signature = "TopicExperiment",
  function(te,value){
    te@phi <- value
    return(te)
  }
)

#### Build Param Functions ####
#' Build Alpha Prior Matrix
#'
#' Constructs the alpha prior matrix for topic modeling. If labels are provided,
#' it creates a sparse matrix based on the provided labels; otherwise, it assigns
#' a default Dirichlet prior.
#'
#' @param te A SingleCellTopicExperiment or SpatialTopicExperiment object.
#' @param labels Character string specifying the column name in `colData(te)`
#'        that contains labels for guiding topic assignment.
#' @param K Integer specifying the number of topics (ignored if labels are provided).
#'
#' @return A matrix representing the alpha prior.
#' @import Matrix
#' @export
build_alpha <- function(te,labels = NULL,K = 10,disease = NULL,disease_type = NULL){
  if (is.null(labels)){
    a <- matrix(1,ncol(te),K)
    colnames(a) <- paste("Topic-",1:K,sep='')
    return(a)
  } else{
    D <- ncol(te)
    K <- length(unique(colData(te)[,labels]))


    ## Dirichlet parameters
    covariates  <- unique(colData(te)[,labels])
    cell_id <- match(colData(te)[,labels],covariates)
    alpha <- 1
    a <- sparseMatrix(i = seq_along(cell_id),
                      j = cell_id,
                      x = alpha,
                      dims = c(D, K))
    a <- as.matrix(a)

    colnames(a) <- covariates

    ## If disease info provided - update matrix accordingly


    if (!is.null(disease) & !is.null(disease_type)){
      a <- a[,rep(1:ncol(a),each =length(unique(disease)))]

      a <- BuildDiseaseAlpha(a,disease)
      colnames(a) <- paste(rep(covariates,
                               each = length(disease_type)),
                               unique(disease_type),
                               sep='_')
    }
    a[a == 0] <- 0.01
    rownames(a) <- colnames(te)
    a <- as.matrix(a)
    return(a)

  }


}

#' Build Beta Prior Matrix
#'
#' Constructs the beta prior matrix for topic modeling, assigning a default prior.
#'
#' @param te A SingleCellTopicExperiment or SpatialTopicExperiment object.
#' @param K Integer specifying the number of topics.
#'
#' @return A matrix representing the beta prior.
#' @export
build_beta <- function(te,K){
  return(matrix(0.001, nrow = nrow(te),ncol = K))
}


#' Build Priors for Topic Modeling
#'
#' Initializes the alpha and beta prior matrices for a topic modeling experiment.
#'
#' @param te A SingleCellTopicExperiment or SpatialTopicExperiment object.
#' @param labels Character string specifying the column name in `colData(te)` for topic guidance.
#' @param guided Logical indicating whether topic modeling should be guided by labels.
#' @param K Integer specifying the number of topics (ignored if guided is TRUE and labels are provided).
#'
#' @return The input `te` object with alpha and beta prior matrices initialized.
#' @export
buildPriors <- function(te,labels = NULL,guided = FALSE,K = 10){
  if (guided){
    if(!is.null(labels)){
      alphaPrior(te) <- build_alpha(te,labels)
      K <- ncol(alphaPrior(te))
    } else {
      stop("labels argument must be defined if guided parameter is set to TRUE")
    }
  } else {
    alphaPrior(te) <- build_alpha(te,labels,K)
  }
  betaPrior(te) <- build_beta(te,K)
  return(te)
}

#' Initialize Sufficient Statistics for Topic Modeling
#'
#' Initializes the sufficient statistics matrices (`ndk` and `nwk`) with zeros.
#'
#' @param te A SingleCellTopicExperiment or SpatialTopicExperiment object.
#' @param K Integer specifying the number of topics.
#'
#' @return The input `te` object with initialized sufficient statistics.
#' @export
build_SufStats <- function(te,K){
  ndk(te) <- matrix(0,ncol(te),K)
  nwk(te) <- matrix(0,nrow(te),K)
  return(te)
}

#' Compute Theta Matrix (Cell-Topic Distributions)
#'
#' Computes the theta matrix, representing the per-cell topic distributions.
#'
#' @param te A SingleCellTopicExperiment or SpatialTopicExperiment object.
#'
#' @return The input `te` object with the `theta` matrix computed.
#' @export
buildTheta <- function(te){
  theta(te) <- get_theta(ndk(te),alphaPrior(te))
  colnames(theta(te)) <- colnames(alphaPrior(te))
  rownames(theta(te)) <- colnames(te)
  return(te)
}


#' Compute Phi Matrix (Gene-Topic Distributions)
#'
#' Computes the phi matrix, representing the per-gene topic distributions.
#'
#' @param te A SingleCellTopicExperiment or SpatialTopicExperiment object.
#'
#' @return The input `te` object with the `phi` matrix computed.
#' @export
buildPhi <- function(te){
  phi(te) <- get_phi(nwk(te),betaPrior(te))
  colnames(phi(te)) <- colnames(alphaPrior(te))
  rownames(phi(te)) <- rownames(te)
  return(te)
}


##Update Show Methods
#### Print Class ####
#' @export
setMethod("show",
          signature = "SingleCellTopicExperiment", function(object) {
            callNextMethod()  # Calls ParentClass's show() first
            cat("TopicPriors: ",
                paste('Alpha ',
                      ifelse(all(object@alphaPrior == 0),'(Empty)',''),
                      ' Beta' ,
                      ifelse(all(object@betaPrior == 0),'(Empty)',''),
                      sep = ''),'\n')
            cat("SuffStats:",
                paste('NDK ',
                      ifelse(all(object@ndk == 0),'(Empty)',''),
                      ' NWK ',
                      ifelse(all(object@ndk == 0),'(Empty)',''),
                      sep = ''),'\n')
            cat("ExpTopics:",
                paste('Theta ',
                      ifelse(all(object@theta == 0),'(Empty)',''),
                      ' Phi ',
                      ifelse(all(object@phi == 0),'(Empty)',''),
                      sep = ''),'\n')  # Adds additional details
          })


#### Print Class ####
#' @export
setMethod("show",
          signature = "SpatialTopicExperiment", function(object) {
            callNextMethod()  # Calls ParentClass's show() first
            cat("TopicPriors: ",
                paste('Alpha ',
                      ifelse(all(object@alphaPrior == 0),'(Empty)',''),
                      ' Beta' ,
                      ifelse(all(object@betaPrior == 0),'(Empty)',''),
                      sep = ''),'\n')
            cat("SuffStats:",
                paste('NDK ',
                      ifelse(all(object@ndk == 0),'(Empty)',''),
                      ' NWK ',
                      ifelse(all(object@ndk == 0),'(Empty)',''),
                      sep = ''),'\n')
            cat("ExpTopics:",
                paste('Theta ',
                      ifelse(all(object@theta == 0),'(Empty)',''),
                      ' Phi ',
                      ifelse(all(object@phi == 0),'(Empty)',''),
                      sep = ''),'\n')  # Adds additional details
          })


###Subset Topic Matrices ####
#' @export
setMethod("[", c("SingleCellTopicExperiment", "ANY", "ANY"),
          function(x, i, j, ..., drop=FALSE) {

            if (missing(i)) i <- TRUE
            if (missing(j)) j <- TRUE
            x <- callNextMethod()

            ## Update Gene-Topic Components
            ii <- SingleCellExperiment:::.convert_subset_index(i, rownames(x))
            betaPrior(x) <- betaPrior(x)[ii,,drop=FALSE]
            nwk(x) <- nwk(x)[ii,,drop=FALSE]
            phi(x) <- phi(x)[ii,,drop=FALSE]

            ## Update Sample-Topic Components
            jj <- SingleCellExperiment:::.convert_subset_index(j, colnames(x))
            alphaPrior(x) <- alphaPrior(x)[jj,,drop=FALSE]
            ndk(x) <- ndk(x)[jj,,drop=FALSE]
            theta(x) <- theta(x)[jj,,drop=FALSE]

            return(x)
          })

#' @export
setMethod("[", c("SpatialTopicExperiment", "ANY", "ANY"),
          function(x, i, j, ..., drop=FALSE) {
                if (missing(i)) i <- TRUE
                if (missing(j)) j <- TRUE
                x <- callNextMethod()

                ## Update Gene-Topic Components
                ii <- SingleCellExperiment:::.convert_subset_index(i, rownames(x))
                betaPrior(x) <- betaPrior(x)[ii,,drop=FALSE]
                nwk(x) <- nwk(x)[ii,,drop=FALSE]
                phi(x) <- phi(x)[ii,,drop=FALSE]

                ## Update Sample-Topic Components
                jj <- SingleCellExperiment:::.convert_subset_index(j, colnames(x))
                alphaPrior(x) <- alphaPrior(x)[jj,,drop=FALSE]
                ndk(x) <- ndk(x)[jj,,drop=FALSE]
                theta(x) <- theta(x)[jj,,drop=FALSE]

                return(x)
})

## Update dimnames
#' @export
setReplaceMethod("dimnames", c("SingleCellTopicExperiment", "list"),
                 function(x, value)
                 {
                   x <- callNextMethod()
                   #Update gene names
                   rownames(x@betaPrior) <- value[[1L]]
                   rownames(x@nwk) <- value[[1L]]
                   rownames(x@phi) <- value[[1L]]
                   #update sample names
                   rownames(x@alphaPrior) <- value[[2L]]
                   rownames(x@ndk) <- value[[2L]]
                   rownames(x@theta) <- value[[2L]]
                   return(x)})

#' @export
setReplaceMethod("dimnames", c("SpatialExperiment", "list"),
                 function(x, value)
                 {
                   # Call the default dimnames replacement from SummarizedExperiment
                   x <- callNextMethod()

                   #Update gene names
                   rownames(x@betaPrior) <- value[[1L]]
                   rownames(x@nwk) <- value[[1L]]
                   rownames(x@phi) <- value[[1L]]
                   #update sample names
                   rownames(x@alphaPrior) <- value[[2L]]
                   rownames(x@ndk) <- value[[2L]]
                   rownames(x@theta) <- value[[2L]]
                   if (!is.null(spatialCoords(x))) {
                     rownames(spatialCoords(x)) <- value[[2L]]  # Sync colnames with spatialCoords
                   }
                   return(x)})


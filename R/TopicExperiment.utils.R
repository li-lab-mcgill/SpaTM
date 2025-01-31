##Topic Experiment Utils


# Getters ####

#### Alpha Prior ####
setGeneric("alphaPrior", function(te) standardGeneric("alphaPrior"))
setGeneric("alphaPrior<-", function(te, values) standardGeneric("alphaPrior<-"))
setMethod(
  f = "alphaPrior",
  signature = c("SingleCellTopicExperiment"),
  function(te) {
    return(te@alphaPrior)
  }
)
setReplaceMethod(
  f = "alphaPrior",
  signature = c("SingleCellTopicExperiment"),
  function(te,values){
    te@alphaPrior <- values
    return(te)
  }
)

## Spatial
setMethod(
  f = "alphaPrior",
  signature = c("SpatialTopicExperiment"),
  function(te) {
    return(te@alphaPrior)
  }
)
setReplaceMethod(
  f = "alphaPrior",
  signature = c("SpatialTopicExperiment"),
  function(te,values){
    te@alphaPrior <- values
    return(te)
  }
)
#### Beta Prior ####
setGeneric("betaPrior", function(te) standardGeneric("betaPrior"))
setGeneric("betaPrior<-", function(te, values) standardGeneric("betaPrior<-"))
setMethod(
  f = "betaPrior",
  signature = "SingleCellTopicExperiment",
  function(te) {
    return(te@betaPrior)
  }
)
setReplaceMethod(
  f = "betaPrior",
  signature = "SingleCellTopicExperiment",
  function(te,values){
    te@betaPrior <- values
    return(te)
  }
)

## Spatial
setMethod(
  f = "betaPrior",
  signature = "SpatialTopicExperiment",
  function(te) {
    return(te@betaPrior)
  }
)
setReplaceMethod(
  f = "betaPrior",
  signature = "SpatialTopicExperiment",
  function(te,values){
    te@betaPrior <- values
    return(te)
  }
)


## Cell by Topics (ndk) ####
setGeneric("ndk", function(te) standardGeneric("ndk"))
setGeneric("ndk<-", function(te, values) standardGeneric("ndk<-"))
setMethod(
  f = "ndk",
  signature = "SingleCellTopicExperiment",
  function(te) {
    return(te@ndk)
  }
)

setReplaceMethod(
  f = "ndk",
  signature = "SingleCellTopicExperiment",
  function(te,values){
    te@ndk <- values
    return(te)
  }
)

##Spatial
setMethod(
  f = "ndk",
  signature = "SpatialTopicExperiment",
  function(te) {
    return(te@ndk)
  }
)

setReplaceMethod(
  f = "ndk",
  signature = "SpatialTopicExperiment",
  function(te,values){
    te@ndk <- values
    return(te)
  }
)
## Gene by Topics (nwk) ####
setGeneric("nwk", function(te) standardGeneric("nwk"))
setGeneric("nwk<-", function(te, values) standardGeneric("nwk<-"))
setMethod(
  f = "nwk",
  signature = "SingleCellTopicExperiment",
  function(te) {
    return(te@nwk)
  }
)
setReplaceMethod(
  f = "nwk",
  signature = "SingleCellTopicExperiment",
  function(te,values){
    te@nwk <- values
    return(te)
  }
)
## Spatial
setMethod(
  f = "nwk",
  signature = "SpatialTopicExperiment",
  function(te) {
    return(te@nwk)
  }
)
setReplaceMethod(
  f = "nwk",
  signature = "SpatialTopicExperiment",
  function(te,values){
    te@nwk <- values
    return(te)
  }
)

## Theta ####
setGeneric("theta", function(te) standardGeneric("theta"))
setGeneric("theta<-", function(te, values) standardGeneric("theta<-"))
setMethod(
  f = "theta",
  signature = "SingleCellTopicExperiment",
  function(te) {
    return(te@theta)
  }
)
setReplaceMethod(
  f = "theta",
  signature = "SingleCellTopicExperiment",
  function(te,values){
    te@theta <- values
    return(te)
  }
)
##Spatial
setMethod(
  f = "theta",
  signature = "SpatialTopicExperiment",
  function(te) {
    return(te@theta)
  }
)
setReplaceMethod(
  f = "theta",
  signature = "SpatialTopicExperiment",
  function(te,values){
    te@theta <- values
    return(te)
  }
)
##Phi ####
setGeneric("phi", function(te) standardGeneric("phi"))
setGeneric("phi<-", function(te, values) standardGeneric("phi<-"))
setMethod(
  f = "phi",
  signature = "SingleCellTopicExperiment",
  function(te) {
    return(te@phi)
  }
)
setReplaceMethod(
  f = "phi",
  signature = "SingleCellTopicExperiment",
  function(te,values){
    te@phi <- values
    return(te)
  }
)
##Spatial
setMethod(
  f = "phi",
  signature = "SpatialTopicExperiment",
  function(te) {
    return(te@phi)
  }
)
setReplaceMethod(
  f = "phi",
  signature = "SpatialTopicExperiment",
  function(te,values){
    te@phi <- values
    return(te)
  }
)

#### Build Param Functions ####

build_alpha <- function(te,labels = NULL,K = 10){
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
                      dims = c(D, K)) %>% as.matrix()
    a[a == 0] <- 0.01
    rownames(a) <- colnames(te)
    colnames(a) <- covariates
    a <- as.matrix(a)
    
    return(a)
    
  }
  
  
}
build_beta <- function(te,K){
  return(matrix(0.001, nrow = nrow(te),ncol = K))
}

buildPriors <- function(te,labels = NULL,guided = FALSE,K = 10){
  if (guided){
    if(!is.null(labels)){
      alphaPrior(te) <- build_alpha(te,labels)
      K <- ncol(alpha(te))
    } else {
      stop("labels argument must be defined if guided parameter is set to TRUE")
    }
  } else {
    alphaPrior(te) <- build_alpha(te,labels,K)
  }
  betaPrior(te) <- build_beta(te,K)
  return(te)
}

build_SufStats <- function(te,K){
  ndk(te) <- matrix(0,ncol(te),K)
  nwk(te) <- matrix(0,nrow(te),K)
  return(te)
}

buildTheta <- function(te){
  theta(te) <- get_theta(ndk(te),alphaPrior(te))
  colnames(theta(te)) <- colnames(alphaPrior(te))
  rownames(theta(te)) <- colnames(te)
  return(te)
}

buildPhi <- function(te){
  phi(te) <- get_phi(nwk(te),betaPrior(te))
  colnames(phi(te)) <- colnames(alphaPrior(te))
  rownames(phi(te)) <- rownames(te)
  return(te)
}
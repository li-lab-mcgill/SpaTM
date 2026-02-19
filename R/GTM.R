
#### Training ####
#' Guided Topic Model (GTM)
#'
#' This function trains a Guided Topic Model (GTM) using the input TopicExperiment (TE) object.
#' GTM optimizes topic distributions using prior knowledge provided through metadata labels.
#'
#' @param scte A `SingleCellTopicExperiment` or `SpatialTopicExperiment` object containing count data and prior matrices.
#' @param K Integer, the number of topics to infer.
#' @param D Integer, the number of cells (documents) in the dataset.
#' @param num_threads Integer, the number of threads to use for parallel computation. Default is 1.
#' @param maxiter Integer, the maximum number of iterations for training. Default is 100.
#' @param verbal Logical, whether to print progress messages. Default is TRUE.
#' @param zero_gamma Logical, whether to initialize gamma values to zero. Default is FALSE.
#' @param rand_gamma Logical, whether to initialize gamma values randomly. Default is TRUE.
#' @param burnin Maximum number of iterations to run during the STM LDA step for each sample (default is 0)
#'
#' @return A `SingleCellTopicExperiment` or `SpatialTopicExperiment` object with updated topic distributions.
#'
#' @import SingleCellExperiment
#' @export
GTM <- function(scte,K,D,num_threads = 1,maxiter = 100,verbal = TRUE,
                  zero_gamma = FALSE,
                  rand_gamma = TRUE,
                  thresh = 1e-8,burnin = 1, nk_init = TRUE){
  train_gtm(counts(scte),
              celltypes = scte$int_cell,
              genes = rowData(scte)$gene_ints,
              alpha = alphaPrior(scte),
              beta = betaPrior(scte),
              K,
              D,
              ndk(scte),
              nwk(scte),
              num_threads,
              maxiter,
              verbal,
              zero_gamma,
              rand_gamma,
              thresh,
              burnin,
              nk_init)
  return(scte)
}

# full-batch

# mini-batch

#### Prediction ####
#' Topic Inference using Pretrained Model
#'
#' This function infers topic distributions for new data using a pretrained topic model.
#'
#' @param scte A `SingleCellTopicExperiment` or `SpatialTopicExperiment` object containing count data.
#' @param num_threads Integer, the number of threads to use for parallel computation. Default is 1.
#' @param maxiter Integer, the maximum number of iterations for inference. Default is 50.
#' @param verbal Logical, whether to print progress messages. Default is TRUE.
#' @param phi A matrix representing the topic-gene distribution from a trained model.
#' @param burnin Maximum number of iterations to run during the STM LDA step for each sample (default is 0)
#'
#' @return A `SingleCellTopicExperiment` or `SpatialTopicExperiment` object with updated inferred topic distributions.
#'
#' @export
inferTopics <- function(scte,num_threads = 1,maxiter = 50,verbal = TRUE,
                         phi,
                        burnin = 1){
  infer_topics_cpp(counts(scte),
               celltypes = scte$int_cell,
               genes = rowData(scte)$gene_ints,
               alpha = 1,ncol(phi),ncol(scte),nrow(scte),ndk(scte),phi,1,maxiter,verbal,burnin)

  return(scte)
}





#' Cell type-specific Gene Expression Deconvolution
#'
#' This function infers Cell type-specific gene expression profiles from bulk data using a pre-trained topic model
#'
#' @param scte A `SingleCellTopicExperiment` or `SpatialTopicExperiment` object containing count data.
#' @param num_threads Integer, the number of threads to use for parallel computation. Default is 1.
#' @param maxiter Integer, the maximum number of iterations for inference. Default is 50.
#' @param verbal Logical, whether to print progress messages. Default is TRUE.
#' @param phi A matrix representing the topic-gene distribution from a trained model.
#'
#' @return A list containing the topic distributions and gene counts for each observation
#'
#' @export
inferGEx <- function(scte,num_threads = 1,maxiter = 50,verbal = TRUE,
                        phi){
  gex_res <- infer_gex_cpp(counts(scte),
                   celltypes = scte$int_cell,
                   genes = rowData(scte)$gene_ints,
                   alpha = 1,ncol(phi),ncol(scte),nrow(scte),ndk(scte),phi,1,maxiter,verbal)

  return(gex_res)
}


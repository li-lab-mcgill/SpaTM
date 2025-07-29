# TODO re-assess which helper functions should be contained in package as opposed to analysis repo

#### Training ####
#' Train a Relational Topic Model (RTM)
#'
#' This function trains a Relational Topic Model (RTM) on a `SpatialTopicExperiment` object.
#' The RTM incorporates relational data (cell neighborhoods) into the topic modeling process.
#'
#' @param scte A `SpatialTopicExperiment` object.
#' @param K An integer specifying the number of topics.
#' @param nbr_list A list specifying the neighbors of individual cells.
#' @param loss_fun An integer specifying the loss function to be used in RTM training (default: 1).
#' @param num_threads An integer indicating the number of threads to use for parallel computation (default: 1).
#' @param maxiter An integer specifying the maximum number of training iterations (default: 100).
#' @param lr A numeric specifying the learning rate for RTM training (default: 1e-5).
#' @param verbal A logical indicating whether to print progress messages (default: TRUE).
#' @param zero_gamma A logical specifying whether to initialize gamma values to zero (default: FALSE).
#' @param rand_gamma A logical indicating whether to randomly initialize gamma values (default: TRUE).
#' @param m_update A logical specifying if the classifier/regression component should be updated (default: TRUE)
#' @param burnin Maximum number of iterations to run during the STM LDA step for each sample (default is 0)
#'
#' @return A `SpatialTopicExperiment` object with updated RTM weights stored in `metadata(scte)[['RTM_weights']]`.
#'
#' @details
#' The RTM extends topic modeling by incorporating relational information, such as neighborhood information.
#' It uses a loss function specified by `loss_fun` to optimize the topic distributions.
#'
#' @seealso [GTM()] for the guided topic model.
#'
#' @importFrom S4Vectors metadata
#' @export
RTM <- function(scte,
                K,
                nbr_list,
                loss_fun = 1,
                num_threads = 1,
                maxiter = 100,
                lr = 1e-5,
                verbal = TRUE,
                zero_gamma = FALSE,
                rand_gamma = TRUE,
                m_update = T,
                burnin = 1){

  #Check if nbr_list is populated
  if (do.call('sum',lapply(nbr_list,nrow)) == 0){
    stop("Your current TopicExperiment does not have any neighbors detected.\nRe-run the get_nbrs function with a different neighbor criteria threshold to ensure positive examples are included.")
  }

  if (!is(scte,'TopicExperiment')){
    stop('This function requires a TopicExperiment object to run correctly. Please create a Spatial or SingleCell variant before using the RTM function.')
  }
  metadata(scte)[['RTM_weights']] <- train_RTM(counts(scte),
                                               scte$int_cell,
                                               rowData(scte)$gene_ints,
                                               nbr_list,
                                               alphaPrior(scte),
                                               betaPrior(scte),
                                               K,
                                               ncol(scte),
                                               ndk(scte),
                                               nwk(scte),
                                               num_threads = num_threads,
                                               maxiter = maxiter,
                                               verbal = verbal,
                                               zero_gamma = zero_gamma,
                                               rand_gamma = rand_gamma,
                                               thresh = 0.00001,
                                               lr =lr,
                                               rho = 50000,
                                               loss_fun = loss_fun,
                                               m_update = m_update,
                                               burnin)
  return(scte)
}







# TODO this assumes only one sample is present in the object. Need to tweak it to return a list of adjacency matrices
###predict full adjacency matrix
#' Get All Predictions from Relational Topic Model (RTM)
#'
#' This function computes predictions using the trained RTM weights stored in the metadata of a `SpatialTopicExperiment` object.
#'
#' @param spe A `SpatialTopicExperiment` object
#' @param loss_fun An integer specifying the loss function to be used for prediction.
#'
#' @return A matrix of predicted values based on the RTM model.
#'
#' @details
#' This function retrieves the trained RTM weights from `metadata(spe)[['RTM_weights']]` and applies them to
#' the topic distribution (`theta(spe)`) to compute relational topic predictions.
#'
#' @seealso [RTM()] for training the Relational Topic Model.
#'
#' @importFrom S4Vectors metadata
#' @export
get_all_pred <- function(spe,loss_fun){

  return(nbr_pred(theta(spe),metadata(spe)[['RTM_weights']],1,loss_fun))
}





#' Smooth labels based on spatial proximity
#'
#' This function smooths the labels of spots in a spatial experiment object based on their spatial proximity to other spots. It assigns new labels to spots based on the labels of their neighboring spots.
#'
#' If `k` is provided, the function computes the spatial distance matrix between spots using their array coordinates (`array_row`, `array_col`) and considers the `k` nearest neighbors for smoothing. Otherwise, the function uses a user-provided list of neighbors for each spot.
#'
#' @param spe A `SpatialTopicExperiment` object. Note that `spe` must include a `spatialCoords` slot with array coordinates (`array_row`, `array_col`).
#' @param labels A vector of current labels corresponding to each spot in `spe`.
#' @param nbr_list A list where each element corresponds to a spot and contains indices of its neighbors.
#' @param k An optional integer specifying the number of nearest neighbors to consider when smoothing the labels. If `NULL`, `nbr_list` is used instead.
#'
#' @return A character vector containing the new labels for each spot after smoothing.
#'
#' @import dplyr
#' @export
rtm_smooth <- function(spe,labels = NULL,nbr_list = NULL,k = NULL){

  if (is.null(labels)){
    stop('A vector of labels must be provided to apply label smoothing.')
  }

  new_labels <- as.character(labels)

  if (!is.null(k)){
    coldata <- colnames(colData(spe))
    if (!all(c('array_row','array_col') %in% coldata)){
      stop('Error: array_row and array_col not found in spe metadata.This function assumes that spatial array coordinates are stored under these two names.')
    }
    else {dist_mat <- get_dist_cpp(as.matrix(colData(spe)[,c('array_row','array_col')]))}
    for(i in 1:ncol(spe)){
      #Consider retaining current spot as a 'neighbor' to incorporate information in decision (tie-breaker)
      cur_nbr <- dist_mat[i,-i]
      cur_idx <- 1:ncol(spe)
      cur_idx <- cur_idx[-i]
      cur_idx <- cur_idx[order(cur_nbr,decreasing = FALSE)[1:k]]
      cur_labels <- labels[cur_idx]
      label_count <- table(cur_labels) %>%
        as.matrix()
      #print(label_count)
      label_count <- label_count/sum(label_count)

      new_labels[i] <- ifelse(any(label_count >= 0.5),
                              rownames(label_count)[which(label_count >= 0.5)],
                              new_labels[i])
    }
  } else if (!is.null(nbr_list)){
    for(i in 1:ncol(spe)){
      cur_nbr <- nbr_list[[i]]
      if(is.matrix(cur_nbr)){
        cur_nbr <- cur_nbr[cur_nbr[,2] == 1,1] + 1
      }
      cur_labels <- labels[cur_nbr]
      label_count <- table(cur_labels) %>%
        as.matrix()
      #print(label_count)
      label_count <- label_count/sum(label_count)

      new_labels[i] <- ifelse(any(label_count >= 0.5),
                              rownames(label_count)[which(label_count >= 0.5)],
                              new_labels[i])
    }
  } else {
    stop('rtm_smooth requires either k or nbr_list to apply label smoothing.')
  }

  return(new_labels)
}




#' Optimize clustering resolution to achieve a target number of clusters
#'
#' This function performs clustering optimization by adjusting the resolution parameter in the Louvain algorithm until the target number of clusters is reached. It iteratively adjusts the resolution parameter, running clustering with each iteration, until the desired number of clusters is obtained or the maximum number of iterations is reached.
#'
#' The function uses an adjacency matrix to construct a graph and performs clustering on the reduced-dimensional representation of the spatial experiment object. The algorithm used for clustering is Louvain by default, but other algorithms can be specified.
#'
#' @param spe A `SpatialTopicExperiment` object.
#' @param target_clusters The desired number of clusters to achieve.
#' @param max_iterations The maximum number of iterations to run for adjusting the resolution (default is 100).
#' @param X The name of the reduced dimension to use for clustering. Default is `'adj'`, which assumes an adjacency matrix is used for clustering.
#' @param k The number of nearest neighbors to consider when constructing the graph (default is 20).
#' @param alg The algorithm to use for clustering. The default is `'louvain'`, but other algorithms can be used (e.g., `'walktrap'`).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{clustering}{A vector of cluster assignments for each spot.}
#'     \item{resolution}{The resolution parameter used to achieve the target number of clusters.}
#'   }
#'
#' @import bluster
#' @importFrom scran clusterCells
#' @export
clust_optim <- function(spe, target_clusters, max_iterations = 100,
                        X = 'adj',k = 20,alg = 'louvain') {
  # Convert adjacency matrix to graph if necessary


  # Initialize resolution bounds and iteration counter
  min_res <- 0.01
  max_res <- 1.0
  iteration <- 0

  while (iteration < max_iterations) {
    # Calculate midpoint resolution
    res <- (min_res + max_res) / 2
    iteration <- iteration + 1

    # Run Louvain clustering with the current resolution
    clustering <-  clusterCells(spe,
                                use.dimred = X,
                                BLUSPARAM = SNNGraphParam(k = k,
                                                          cluster.fun = alg,
                                                          cluster.args = list(resolution = res)))
    # Get the number of clusters
    num_clusters <- length(unique(clustering))

    # Check if the result matches the target number of clusters
    if (num_clusters == target_clusters) {
      return(list(clustering = clustering, resolution = res))
    } else if (num_clusters < target_clusters) {
      # Increase resolution for more clusters
      min_res <- res
    } else {
      # Decrease resolution for fewer clusters
      max_res <- res
    }
  }

  warning("Maximum iterations reached without finding the exact number of clusters")
  print(res)
  print(num_clusters)
  return(list(clustering = clustering, resolution = res))
}

#' Min-Max Normalization
#'
#' This function performs min-max normalization on an input matrix, scaling the values to be between 0 and 1.
#'
#' @param adj_mat A numeric matrix to be normalized.
#'
#' @return A numeric matrix with values scaled between 0 and 1.
#'
#' @examples
#' # Example: Normalize a matrix
#' mat <- matrix(runif(100, 0, 10), nrow = 10, ncol = 10)
#' norm_mat <- minmax_norm(mat)
#' print(norm_mat)
#'
#' @export
minmax_norm <- function(adj_mat){
  min_val <- min(adj_mat)
  max_val <- max(adj_mat)
  adj_norm <- (adj_mat - min_val) / (max_val - min_val)
  return(adj_norm)
}





#' Get neighboring spots based on spatial distance and grouping
#'
#' This function identifies neighboring spots for each spot in a spatial experiment object based on their spatial coordinates. It can consider spots from the same sample and group, and allows for the inclusion of loss functions to modify neighbor selection.
#'
#' If `group_by` is provided, the function selects neighbors from the same sample and group, and may sample neighbors based on spatial proximity. If no `group_by` is given, neighbors are selected within a given distance threshold (`dist`) from the current spot, and optionally modified by a loss function.
#'
#' @param spe A `SpatialTopicExperiment` object.
#' @param samples The column name in `colData(spe)` representing sample identifiers.
#' @param cell_ids The column name in `colData(spe)` representing cell identifiers.
#' @param group_by The column name in `colData(spe)` to group the spots by. If `NULL`, no grouping is performed, and neighbors are selected based on distance alone.
#' @param dist The maximum spatial distance within which to consider spots as neighbors (default is 1).
#' @param loss_fun A numeric value to control the use of loss functions:
#'   - 0: No loss function applied.
#'   - 1: Applies a loss function that includes randomly selected negative neighbors based on spatial distance.
#'   - 2: Applies a Euclidean distance-based loss function.
#'
#' @return A list where each element corresponds to a spot, and contains its neighboring spots, represented by a matrix with indices and weights.
#'   If a loss function is used, additional information will be appended to the matrix.
#'
#' @import SpatialExperiment
#' @export
get_nbrs <- function(spe,samples,cell_ids,group_by = NULL,dist = 1,loss_fun = 1){
  coldata <- colnames(colData(spe))
  if (!all(c('array_row','array_col') %in% coldata)){
    stop('Error: array_row and array_col not found in spe metadata.This function assumes that spatial array coordinates are stored under these two names.')
  }
  else {dist_mat <- get_dist_cpp(as.matrix(colData(spe)[,c('array_row','array_col')]))}
  nbr_list <- list()

  if (!is.null(group_by)){
    for (i in colData(spe)[,cell_ids]){

      cur_sample <- colData(spe)[i,samples]
      cur_group <- colData(spe)[i,group_by]
      nbrs <- which(colData(spe)[,samples] == cur_sample & colData(spe)[,group_by] == cur_group) - 1
      #Remove current example
      nbrs <- nbrs[nbrs != (i-1)]
      n_sam <- 50
      if (length(nbrs) > n_sam){
        nbrs <- sample(nbrs,n_sam,prob = exp(-dist_mat[i,nbrs+1])/sum(exp(-dist_mat[i,nbrs+1])))
      }
      nbr_list[[i]] <- cbind(nbrs,rep(1,length(nbrs)))
      n <- length(nbrs)
     if (loss_fun == 1 & n > 0){
       nbr_list[[i]] <- rbind(nbr_list[[i]],
                              cbind(sample((1:ncol(spe))[-i],
                                           n,
                                           prob = 1/(1 + exp(-dist_mat[i,-i])))-1,
                                    rep(0,n)))
     }
    }
  }else {
    for (i in colData(spe)[,cell_ids]){
      cur_nbrs <- matrix(0,nrow = 0,ncol = 2)
      cur_sample <- colData(spe)[i,samples]
      rw <- colData(spe)$array_row[i]
      cl <- colData(spe)$array_col[i]
      nbrs <- which(dist_mat[i,] <= dist & colData(spe)[,samples] == cur_sample) - 1
      #Remove current example
      nbrs <- nbrs[nbrs != (i-1)]
      n <- length(nbrs)
      nbr_list[[i]] <- rbind(cur_nbrs,cbind(nbrs,rep(1,n)))
      if (loss_fun != 0 & n > 0){
        if (loss_fun == 1){
          nbr_list[[i]] <- rbind(nbr_list[[i]],
                                 cbind(sample((1:ncol(spe))[-i],
                                              n,
                                              prob = 1/(1 + exp(-dist_mat[i,-i])))-1,
                                       rep(0,n)))
        }
        else if (loss_fun == 2 ){
          ##Euclidean Distance Prediction
          nbr_list[[i]] <- rbind(nbr_list[[i]],
                                 cbind(sample((1:ncol(spe))[-i],
                                              n)-1,rep(0,n)))
          nbr_list[[i]] <- cbind(nbr_list[[i]],
                                 log(dist_mat[i,nbr_list[[i]][,1] + 1]))
          #
        }
      } else if (n == 0){
        nbr_list[[i]] <- rbind(cbind(cur_nbrs,rep(0,n)),cbind(nbrs,rep(1,n),rep(0,n)))
      }
    }
  }
  if (do.call('sum',lapply(nbr_list,nrow)) == 0){
    message('Warning. Your current parametrization returned 0 neighbors for all samples. Consider a different distance threshold or grouping by a covariate.')
  }
  return(nbr_list)
}

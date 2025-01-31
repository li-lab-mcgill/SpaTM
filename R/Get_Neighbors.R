library(spatialLIBD)
library(tidyverse)
library(Matrix)

# spe - SPE object
# samples - colData column name for sample  IDs
# cell_ids - colData column name for cell IDs
#cpp function that takes a matrix with two columns of coordinates and returns a matrix of distances
get_spots <- function(spe,samples,cell_ids,group_by = NULL,dist = 1,loss_fun = 0){
  coldata <- colnames(colData(spe))
  if (!all(c('array_row','array_col') %in% coldata)){
    stop('Error: array_row and array_col not found in spe metadata.This function assumes that spatial array coordinates are stored under these two names.')
  }
  else {dist_mat <- get_dist_cpp(as.matrix(colData(spe)[,c('array_row','array_col')]))}
  nbr_list <- list()

  if (!is.null(group_by)){
    for (i in colData(spe)[,cell_ids]){
      cur_nbrs <- matrix(0,nrow = 0,ncol = 2)
      cur_sample <- colData(spe)[i,samples]
      cur_group <- colData(spe)[i,group_by]
      candid_nbrs <- which(colData(spe)[,samples] == cur_sample & colData(spe)[,group_by] == cur_group) - 1
      #Remove current example
      cur_nbrs[cur_nbrs != (i-1)]
      n_sam <- 50
      if (length(cur_nbrs) > n_sam){
        cur_nbrs <- sample(cur_nbrs,n_sam,prob = exp(-dist_mat[i,cur_nbrs])/sum(exp(-dist_mat[i,cur_nbrs])))
      }
      nbr_list[[i]] <- rbind(cur_nbrs,rep(1,length(cur_nbrs)))
      n  <- length(cur_nbrs)
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

  return(nbr_list)
}


##########
rtm_smooth <- function(spe,labels,nbr_list,k = NULL){
  if (!is.null(k)){
    coldata <- colnames(colData(spe))
    if (!all(c('array_row','array_col') %in% coldata)){
      stop('Error: array_row and array_col not found in spe metadata.This function assumes that spatial array coordinates are stored under these two names.')
    }
    else {dist_mat <- get_dist_cpp(as.matrix(colData(spe)[,c('array_row','array_col')]))}
  }
  new_labels <- as.character(labels)

  if (!is.null(k)){
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
  } else{
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
  }
  
  return(new_labels)
}

###################3
check_nbr <- function(spe,nbr_set,cell_id,sample_id){
  spe$nbr <- NA
  spe$nbr[cell_id] <- 'Cur Cell'
  spe$nbr[nbr_set[nbr_set[,2] == 1,1] + 1] <- 'NBR'
  spe$nbr[nbr_set[nbr_set[,2] == 0,1] + 1] <- 'Neg. Sample'
  spe$nbr <- as.factor(spe$nbr)
  vis_clus(spe,
           sample_id,
           clustervar = 'nbr')
}

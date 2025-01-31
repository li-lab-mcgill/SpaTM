library(scran)
library(scater) 
library(tidyverse)   
library(Matrix) 
library(spatialLIBD)
library(DescTools)
library(cluster)


###predict full adjacency matrix
get_all_pred <- function(spe,rtm_weights,loss_fun){
  
  return(nbr_pred(theta(spe),rtm_weights,1,loss_fun))
}


## Build ROC Curve

build_ROC <- function(all_pred,grd_adj,nbr_list,full_adj = FALSE,
                      titles = "Def. Title"){
  roc_list <- list()
  pred_adj <- all_pred
  model_name <- titles
  perf_df <- matrix(0,nrow = 0,ncol = 5)
  for (thresh in seq(0,1,length.out = 250)){
    tn <- 0
    tp <- 0
    fn <- 0
    fp <- 0
    if (full_adj){
      pred_adj[pred_adj > thresh] <- 1
      pred_adj[pred_adj <= thresh] <- 0
      
      tp <- length(which((pred_adj + grd_adj) == 2))
      tn <- length(which((pred_adj + grd_adj) == 0))
      fp <- length(which((pred_adj - grd_adj) == 1))
      fn <- length(which((pred_adj - grd_adj) == -1))
    }
    else {
      for (cell in 1:ncol(pred_adj)){
        if (nrow(nbr_list[[cell]]) == 0){
          next
        }
        grd_adj <- nbr_list[[cell]][,2]
        nbr_ids <- nbr_list[[cell]][,1]+1
        sub_adj <-  pred_adj[cell,nbr_ids]
        sub_adj[sub_adj > thresh] <- 1
        sub_adj[sub_adj <= thresh] <- 0
        
        tp <- tp + length(which((sub_adj + grd_adj) == 2))
        tn <- tn + length(which((sub_adj + grd_adj) == 0))
        fp <- fp + length(which((sub_adj - grd_adj) == 1))
        fn <- fn + length(which((sub_adj - grd_adj) == -1))
        
        if (sum(c(tp,tn,fp,fn)) == 0){stop('Error with TP/TN/FP/FN Calculation!')}
      }
      
      
    }
    
    
    
    tpr <- tp/(tp+fn)
    fpr <- fp/(fp+tn)
    sens <- tp/(tp+fn) #recall
    precis <- tp/(tp+fp)
    perf_df <-rbind(perf_df,c(tpr,fpr,sens,precis,thresh))
  }
  roc_list[[model_name]] <- perf_df[,1:4]
  #print(DescTools::AUC(perf_df[,2],perf_df[,1]))
  
  
  return(roc_list)
}





plot_clusters <- function(spe,anno_label = '',test_sample = '',pred_list,clusters = 7,
                          model_name = 'Default Name'){
  #Run kmeans
  temp <- kmeans(pred_list,clusters)
  #store in spe
  spe$clust <- temp$cluster
  cur_ari <- aricode::ARI(spe$clust,colData(spe)[,anno_label])
  print(paste("ARI: ",round(cur_ari,3)))
  ari_plot <- vis_clus(spe,
                       clustervar = "clust",
                       sampleid = test_sample) + 
    labs(subtitle = paste(model_name," ARI: ",round(cur_ari,3),sep = ''))
  return(ari_plot)
}


rtm_cluster <- function(spe,adj_mat,method = 'Kmeans',clusters = 10){
  labels <- rep('',ncol(spe))
  if (method == 'Kmeans'){
    
  } else if (method == 'Louvain'){
    
  } else {
    stop("Only Kmeans or Louvain are accepted clustering methods")
  }
  
  return(labels)
}

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


library(bluster)

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

minmax_norm <- function(adj_mat){
  min_val <- min(adj_mat)
  max_val <- max(adj_mat)
  adj_norm <- (adj_mat - min_val) / (max_val - min_val)
  return(adj_norm)
}


get_perf <- function(spe,df,ground_truth,nbr_list,smooth = TRUE, k= NULL){
  ari_perf <- apply(df,2,function(a){
    if (smooth){
      a <- rtm_smooth(spe,a,nbr_list,k)
    }
    aricode::ARI(a[!is.na(a)],ground_truth[!is.na(a)])
  })
  nmi_perf <- apply(df,2,function(a){
    if (smooth){
      a <- rtm_smooth(spe,a,nbr_list,k)
    }
    aricode::NMI(a[!is.na(a)],ground_truth[!is.na(a)])
  })
  

  asw_perf <- apply(df,2,function(a){
    if (smooth){
      a <- rtm_smooth(spe,a,nbr_list,k)
      a <- as.numeric(a)
    }
    keep_idx <- which(!is.na(a))
    si_df <- cluster::silhouette(a[keep_idx],dist = dist(theta(spe_test[,keep_idx])))
    mean(si_df[,3])
  })
  out_df <- data.frame(ARI = ari_perf,
                       NMI = nmi_perf,
                       ASW = asw_perf)
  rownames(out_df) <- colnames(df)
  return(out_df)
}

plot_perf <- function(spe, df, perf, nbr_list, smooth = TRUE, k = NULL) {
  all_plots <- list()
  for (i in 1:ncol(df)) {
    spe$cur <- df[, i]
    if (smooth) {
      spe$cur <- rtm_smooth(spe, spe$cur, nbr_list, k)
    }
    plot_title <- paste(
      spe$sample_id[1],
      " ",
      colnames(df)[i],
      "\n",
      "ARI:",
      round(perf[i, 1], 3),
      " NMI:",
      round(perf[i, 2], 3),
      " ASW:",
      round(perf[i,3],3),
      sep = ''
    )
    all_plots[[i]] <- vis_clus(spe, spe$sample_id[1], clustervar = 'cur') +
        labs(title = plot_title)
  }
  return(all_plots)
}

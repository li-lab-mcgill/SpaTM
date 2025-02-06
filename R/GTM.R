
#### Training ####
GTM <- function(scte,K,D,num_threads = 1,maxiter = 100,verbal = TRUE,
                  zero_gamma = FALSE,
                  rand_gamma = TRUE){
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
              rand_gamma)
  return(scte)
}

# full-batch

# mini-batch

#### Prediction ####
inferTopics <- function(scte,num_threads = 1,maxiter = 50,verbal = TRUE,
                         phi){
  infer_topics_cpp(counts(scte),
               celltypes = scte$int_cell,
               genes = rowData(scte)$gene_ints,
               alpha = 1,ncol(phi),ncol(scte),nrow(scte),ndk(scte),phi,1,maxiter,verbal)

  return(scte)
}

#### Visualization ####

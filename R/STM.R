
library(torch)
library(Rcpp)
library(RcppArmadillo)
library(utils)
library(progress)
#Build mlp
build_mlp <- function(layers = 1,d_in,d_hidden,d_out,device = torch_device('cpu'),
                      dummy_topic = FALSE){
  ## Dummy Topic Edit
  if (dummy_topic){
    d_in <- d_in - 1
  }
  ###
  if (!layers %in% 0:3){
    stop('Error: MLP is set to use 3 hidden layers at most!')
  }

  else if (layers == 0){
    mlp <- nn_sequential(
      nn_linear(d_in,d_out),
      nn_softmax(dim = 2)
    )
  } else if (layers == 1){
    mlp <- nn_sequential(
      nn_linear(d_in,d_hidden),
      nn_relu(),
      nn_linear(d_hidden,d_out),
      nn_softmax(dim = 2)
    )
  } else if (layers == 2){
    mlp <- nn_sequential(
      nn_linear(d_in,d_hidden[1]),
      nn_relu(),
      nn_linear(d_hidden[1],d_hidden[2]),
      nn_relu(),
      nn_linear(d_hidden[2],d_out),
      nn_softmax(dim = 2)
    )
  } else if (layers == 3){
    mlp <- nn_sequential(
      nn_linear(d_in,d_hidden[1]),
      nn_relu(),
      nn_linear(d_hidden[1],d_hidden[2]),
      nn_relu(),
      nn_linear(d_hidden[2],d_hidden[3]),
      nn_relu(),
      nn_linear(d_hidden[3],d_out),
      nn_softmax(dim = 2)
    )
  }
  return(mlp)
}

mlp_data <- dataset(
  name = 'stm_data',
  initialize = function(x,y){
    self$x <- x
    self$y <- y
  },

  .getbatch = function(i){
    x <- self$x[i,]
    y <- self$y[i]
    return(list(x = x,y = y))
  },
  .length = function(){
    dim(self$y)
  }
)




# Define fitting function
fit_mlp <- function(mlp, train_dl, lr, max_epoch,class_prop) {
  opt <- optim_adam(mlp$parameters, lr = lr,weight_decay = lr/2)
  loss_fn <- nn_cross_entropy_loss(weight = class_prop)

  for (i in 1:max_epoch) {
    coro::loop(
      for (batch in train_dl) {
        x_train <- batch$x
        y_train <- batch$y
        # Train_Pred
        y_pred <- mlp(x_train)
        # %>%
        #   nnf_log_softmax(dim = 2)
        loss <- loss_fn(y_pred, y_train)

        # Backpropagation
        opt$zero_grad()
        loss$backward()

        # Update weights
        opt$step()
      }
    )
    #if (i %% 10 == 0){print(loss$item())}
  }

  return(mlp)
}

# Define fitting function
alt_fit_mlp <- function(mlp, x,y, lr, max_epoch,batch_size) {
  opt <- optim_adam(mlp$parameters, lr = lr,weight_decay = lr/2)
  loss_fn <- nn_cross_entropy_loss()


  for (i in 1:max_epoch) {
    batch_ids <- sample(1:dim(y))
    batches <- split(batch_ids,ceiling(seq_along(batch_ids)/batch_size))
    for (batch in batches){
      x_train <- x[batch,]
      y_train <- y[batch,]
      # Train_Pred
      y_pred <- mlp(x_train)
      # %>%
      #   nnf_log_softmax(dim = 2)
      loss <- loss_fn(y_pred, y_train)

      # Backpropagation
      opt$zero_grad()
      loss$backward()

      # Update weights
      opt$step()
    }

  }

  return(mlp)
}




pred_mlp <- function(mlp,
                     x,device = torch_device('cpu'),
                     dummy_topic = FALSE){
  ## Dummy Topic Edit

  if (dummy_topic){
    x <- x[,1:(ncol(x)-1)]
  }
  ####
  temp <- mlp(torch::torch_tensor(x,device = device))
  temp <- temp$to(dtype = torch_float()) %>%
    torch_squeeze()
  temp <- as.matrix(temp$cpu()) %>%
    apply(1,which.max)
  return(temp)

}

#

STM_Torch <- function(spe, num_threads = 1L,
                maxiter = 100L, verbal = TRUE,
                zero_gamma = TRUE,
                rand_gamma = FALSE, thresh = 0.00001,
                lr = 0.0001,
                mlp,
                mlp_layers,
                mlp_epoch = 1000,
                spe_val,
                device = torch_device('cpu'),
                balanced_class = TRUE,
                dummy_topic = FALSE){
  elbo_mat <- matrix(0,max_iter,5)
  old_elbo <- -9999999
  cur_elbo <- 0

  tensor_labels <- torch_tensor(spe$layer + 1,requires_grad = TRUE,device = device) %>%
    torch_squeeze()
  tensor_labels <- tensor_labels$to(dtype = torch_long())

  acc_mat <- matrix(0,maxiter,2)
  for(i in 1:maxiter){
    #print(paste("Iteration: ",i))
    #saveRDS(list(mlp = mlp,spe = spe),'cur_model_output.rds')
    mlp_parameters <- lapply(mlp$parameters,function(a){
      t(as.matrix(a$cpu()))})
    cur_elbo <- stm_torch_estep(counts(spe),
                                spe$int_cell,
                                rowData(spe)$gene_ints,
                                alphaPrior(spe),
                                betaPrior(spe),
                                K,
                                ncol(spe),
                                zero_gamma,
                                rand_gamma,
                                spe$layer,
                                ndk(spe),
                                nwk(spe),
                                elbo_mat,
                                mlp_layers,
                                mlp_parameters,
                                cur_elbo,
                                num_threads,
                                i-1,
                                dummy_topic) ##Added dummy topic

    #Calc N_mean
    #M STEP

    ndk_prop <- torch::torch_tensor(ndk(spe)/rowSums(ndk(spe)),device = device)
    ##Added Dummy Topic Edit
    if (dummy_topic){
      ndk_prop <- ndk_prop[,1:(K-1)]
    }
    ####

    train_dl <- mlp_data(ndk_prop,tensor_labels) %>%
      dataloader(1000,TRUE,num_workers = 0)
    if (balanced_class){
      class_prop <- as.numeric(table(spe$layer + 1)/ncol(spe))^(-1) %>%
        torch_tensor() %>%
        torch_squeeze()

    } else{
      num_layer <- length(unique(spe$layer))
      class_prop <- rep(1,num_layer) %>%
        torch_tensor() %>%
        torch_squeeze()

    }

    mlp <- fit_mlp(mlp,train_dl,lr,mlp_epoch,class_prop)
    if (i == maxiter){mlp <- fit_mlp(mlp,train_dl,lr,3000, class_prop)}
    cur_pred <- mlp(ndk_prop)
    ce_loss <- nnf_cross_entropy(cur_pred,tensor_labels)
    ce_loss <- as.numeric(ce_loss$cpu())
    #Calculate ELBO
    cur_elbo <- cur_elbo -ce_loss/ncol(counts(spe))
    dif <- abs(old_elbo-cur_elbo)


    elbo_mat[i,5] <- -ce_loss/ncol(counts(spe))
    # if (abs(cur_elbo - old_elbo) <= thresh & i < maxiter){
    #   mlp <- fit_mlp(mlp,train_dl,lr,3000,class_prop)
    #   print("ELBO converged early")
    #   break
    # }
    old_elbo <- cur_elbo
    #print(cur_elbo)
    elbo_mat[i,4] <- cur_elbo
    cur_elbo <- 0
    ## Check Accuracy

    spe <- buildTheta(spe)
    spe <- buildPhi(spe)


    spe <- scGTMPredict(spe,1,100,verbal = FALSE,phi(spe))
    spe <- buildTheta(spe)


    train_pred <- pred_mlp(mlp,theta(spe),device,dummy_topic)

    train_acc <- length(which(train_pred == as.numeric(spe$spatialLIBD)))*
      100/ncol(spe)
    #print(paste("train Acc: ",train_acc, sep = ''))

    ###
    spe_val <- scGTMPredict(spe_val,1,100,verbal = FALSE,phi(spe))
    spe_val <- buildTheta(spe_val)


    val_pred <- pred_mlp(mlp,theta(spe_val),device,dummy_topic)

    val_acc <- length(which(val_pred == as.numeric(spe_val$spatialLIBD)))*
      100/ncol(spe_val)
    #print(paste("Val Acc: ",val_acc, sep = ''))

    acc_mat[i,] <- c(train_acc,val_acc)
    progress_bar(i,maxiter,elbo_mat[i,4],train_acc,val_acc)
    if (dif < 0.0000000001){
      #write.csv(acc_mat,'debug_acc.csv')
      #write.csv(elbo_mat,'debug_elbo.csv')
      #TODO Find a nicer way to do the below segment
      return(mlp)
    }
  }
  #write.csv(acc_mat,'debug_acc.csv')
  #write.csv(elbo_mat,'debug_elbo.csv')
  #TODO Find a nicer way to do the below segment
  return(mlp)

}



STM_Torch_burnin <- function(spe, num_threads = 1L,
                       maxiter = 100L, verbal = TRUE,
                       zero_gamma = TRUE,
                       rand_gamma = FALSE, thresh = 0.00001,
                       lr = 0.0001,
                       mlp,
                       mlp_layers,
                       mlp_epoch = 1000,
                       spe_val){

  old_elbo <- -9999999
  cur_elbo <- 0
  elbo_mat <- matrix(0,max_iter,5)
  tensor_labels <- torch_tensor(spe$layer + 1,requires_grad = TRUE) %>%
    torch_squeeze()
  tensor_labels <- tensor_labels$to(dtype = torch_long())

  acc_mat <- matrix(0,maxiter,2)
  for(i in 1:maxiter){
    mlp_parameters <- lapply(mlp$parameters,function(a){
      t(as.matrix(a$cpu()))})
    cur_elbo <- stm_torch_estep(counts(spe),
                                spe$int_cell,
                                rowData(spe)$gene_ints,
                                alphaPrior(spe),
                                betaPrior(spe),
                                K,
                                ncol(spe),
                                zero_gamma,
                                rand_gamma,
                                spe$layer,
                                ndk(spe),
                                nwk(spe),
                                elbo_mat,
                                mlp_layers,
                                mlp_parameters,
                                cur_elbo,
                                num_threads,
                                i-1)
  }
}



## Full STM Wrapper
## TODO add class imbalance to GLM
## TODO Use MLP wrapper within this function for ease of use
## TODO Resolve Torch openMP bug!
STM <- function(spe, num_threads = 1L,
                      maxiter = 100L, verbal = TRUE,
                      zero_gamma = TRUE,
                      rand_gamma = FALSE, thresh = 0.00001,
                      lr = 0.0001,
                      mlp = NULL,
                      mlp_layers = NULL,
                      mlp_epoch = NULL,
                      spe_val = NULL,
                      device = NULL,
                      balanced_class = TRUE,
                      dummy_topic = FALSE,
                      burnin = 5){
  if (is.null(mlp)){
    message('Running GLM-STM')
    metadata(spe)[['STM_Weights']] <- train_stm(counts(spe),
                                                spe$int_cell,
                                                rowData(spe)$gene_ints,
                                                spe$layer,
                                                alphaPrior(spe),
                                                betaPrior(spe),
                                                K,
                                                ncol(spe),
                                                ndk(spe),
                                                nwk(spe),
                                                num_threads,
                                                maxiter,
                                                verbal,
                                                zero_gamma,
                                                rand_gamma,
                                                thresh,
                                                lr)
    return(spe)
  }
  else{
    message('Running Torch-STM')
    if (is.null(spe_val)){
      stop('Need to have a validation set to ensure there is no overfitting')
    }
    STM_Torch_burnin(spe_train,
                     num_threads,
                     burnin,
                     verbal,
                     zero_gamma,
                     rand_gamma,
                     thresh,
                     lr,
                     mlp,
                     mlp_layers,
                     mlp_epoch = 0,
                     spe_val)

    metadata(spe)[['STM_MLP']] <- STM_Torch(spe_train,
                                            num_threads,
                                            max_iter,
                                            zero_gamma,
                                            rand_gamma,
                                            thresh,
                                            lr,
                                            mlp,
                                            mlp_layers,
                                            mlp_epoch,
                                            spe_val)
  }

}

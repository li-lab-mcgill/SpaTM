# TODO hide the mlp step from the user to simplify code
# TODO write a function to save and load the spe-mlp objects

#' Build a Multi-Layer Perceptron (MLP) Model
#'
#' This function constructs a multi-layer perceptron (MLP) model with up to 3 hidden layers. The model can be used for tasks such as classification, and the number of hidden layers and their sizes are configurable.
#' If a dummy topic is used, the input dimension is adjusted accordingly.
#'
#' @param layers Integer indicating the number of hidden layers. Must be between 0 and 3 (default is 1).
#' @param d_in Integer specifying the input dimension (number of input features).
#' @param d_hidden A vector of integers specifying the dimensions of the hidden layers. The length of the vector determines the number of hidden layers.
#'   For example, if `layers = 2`, `d_hidden` should have two elements (e.g., `c(64, 32)`).
#' @param d_out Integer specifying the output dimension (e.g., the number of classes or topics).
#' @param device The device on which to run the model. Defaults to `"cpu"`. Can be set to `"cuda"` for GPU acceleration if available.
#' @param dummy_topic Logical indicating whether to include a dummy topic in the input (default is `FALSE`).
#'   If `TRUE`, the input dimension (`d_in`) is reduced by 1 to account for the dummy topic.
#'
#' @return An MLP model (of class `nn_module`) for the specified architecture. The model consists of sequential layers, including:
#'   - `nn_linear`: Linear transformations between layers.
#'   - `nn_relu`: ReLU activation functions.
#'   - `nn_softmax`: Softmax activation function applied to the output layer.
#'
#' @examples
#' # Example 1: Build a 1-hidden-layer MLP with 100 input features
#' # 50 hidden units, and 10 output units.
#' mlp_model <- build_mlp(layers = 1, d_in = 100, d_hidden = 50, d_out = 10)
#'
#' # Example 2: Build a 2-hidden-layer MLP with 100 input features,
#' # 64 and 32 hidden units, and 10 output units.
#' mlp_model <- build_mlp(layers = 2, d_in = 100,
#'  d_hidden = c(64, 32), d_out = 10)
#'
#' # Example 3: Build a 3-hidden-layer MLP with a dummy topic.
#' mlp_model <- build_mlp(layers = 3, d_in = 100,
#' d_hidden = c(64, 32, 16), d_out = 10, dummy_topic = TRUE)
#'
#' @import torch
#' @export
build_mlp <- function(layers = 1,d_in,d_hidden,d_out,device = NULL,
                      dummy_topic = FALSE){
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop('You do not have torch installed. Please install it before trying to use the STM-Torch workflow')
  }
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

#' MLP Data Dataset
#'
#' This function defines a dataset for the MLP model, including methods for initialization, getting batches, and getting the length of the dataset.
#'
#' @param x A tensor representing the input data.
#' @param y A tensor representing the target data.
#'
#' @return A dataset object for the MLP model.
#'
#' @import torch
#' @export
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
#' Fit MLP Model
#'
#' This function trains a multi-layer perceptron (MLP) model using the provided training data loader, learning rate, maximum number of epochs, and class proportions.
#'
#' @param mlp An MLP model (of class `nn_module`).
#' @param train_dl A data loader for the training data.
#' @param lr Learning rate for the optimizer.
#' @param max_epoch Maximum number of epochs for training.
#' @param class_prop Class proportions for weighting the loss function.
#'
#' @return The trained MLP model.
#'
#' @import torch
#' @importFrom coro loop
#' @export
fit_mlp <- function(mlp, train_dl, lr, max_epoch,class_prop) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop('You do not have torch installed. Please install it before trying to use the STM-Torch workflow')
  }
  opt <- optim_adam(mlp$parameters, lr = lr,weight_decay = lr/2)
  loss_fn <- nn_cross_entropy_loss(weight = class_prop)

  for (i in 1:max_epoch) {
    coro::loop(
      for (batch in train_dl) {
        opt$zero_grad()
        x_train <- batch$x
        y_train <- batch$y
        # Train_Pred
        y_pred <- mlp(x_train)
        # %>%
        #   nnf_log_softmax(dim = 2)
        loss <- loss_fn(y_pred, y_train)

        # Backpropagation

        loss$backward()

        # Update weights
        opt$step()
      }
    )
    #if (i %% 10 == 0){print(loss$item())}
  }

  return(mlp)
}



## TODO Add GLM Prediction function
#' Predict with MLP Model
#'
#' This function generates predictions using a trained multi-layer perceptron (MLP) model and the provided input data.
#'
#' @param mlp A trained MLP model (of class `nn_module`).
#' @param x A matrix representing the input data.
#' @param device The device on which to run the model. Defaults to `"cpu"`.
#' @param dummy_topic Logical indicating whether to include a dummy topic in the input (default is `FALSE`).
#'
#' @return A vector of predicted class labels.
#'
#' @import torch
#' @export
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
#' Train a Spatial Topic Model (SpaTM-S) using Torch
#'
#' This function trains a SpaTM-S using a multi-layer perceptron (MLP) model. Note that this function should only be used by the \link{STM} wrapper function for ease-of-use.
#' @param spe A `SpatialTopicExperiment` object.
#' @param labels A vector of spatial (or phenotype) labels for the cells.
#' @param num_threads Number of threads to use for parallel computation (default is 1).
#' @param maxiter Maximum number of iterations for model training (default is 100).
#' @param verbal Logical indicating whether to display progress messages (default is `TRUE`).
#' @param zero_gamma Logical indicating whether to initialize the gamma parameters to zero (default is `TRUE`).
#' @param rand_gamma Logical indicating whether to randomize the gamma parameters (default is `FALSE`).
#' @param thresh Convergence threshold for stopping the model training (default is 0.00001).
#' @param lr Learning rate for optimization (default is 0.0001).
#' @param mlp A model definition for the multi-layer perceptron (MLP).
#' @param mlp_layers A vector specifying the number of units in each layer of the MLP.
#' @param mlp_epoch The number of epochs for training the MLP.
#' @param spe_val A validation set used for monitoring overfitting in SpaTM-S (default `NUll`).
#' @param device The device for training (e.g., "cpu" or "cuda") when using SpaTM-S. If `NULL`, the default device is used. Note that we are still working on reliably GPU enabling the code.
#' @param balanced_class Logical indicating whether to balance the class distribution during model training (default is `TRUE`).
#' @param burnin Maximum number of iterations to run during the STM LDA step for each sample (default is 0)
#' @param dummy_topic Logical indicating whether to add a dummy topic for regularization (default is `FALSE`).
#'
#' @return The trained MLP model (of class `nn_module`).
#'
#' @seealso \link{STM}
#'
#' @import torch dplyr SingleCellExperiment
#' @export
STM_Torch <- function(spe,
                      labels,
                      num_threads = 1L,
                      maxiter = 100L,
                      verbal = TRUE,
                      zero_gamma = TRUE,
                      rand_gamma = FALSE,
                      thresh = 0.00001,
                      lr = 0.0001,
                      mlp,
                      mlp_layers,
                      mlp_epoch = 1000,
                      spe_val = NULL,
                      val_labels = NULL,
                      device = torch_device('cpu'),
                      balanced_class = TRUE,
                      burnin = 1,
                      dummy_topic = FALSE,
                      final_run_iter = 100){
  elbo_mat <- matrix(0,maxiter,5)
  old_elbo <- -9999999
  cur_elbo <- 0
  K <- ncol(alphaPrior(spe))
  tensor_labels <- torch_tensor(labels + 1,requires_grad = TRUE,device = device) %>%
    torch_squeeze()
  tensor_labels <- tensor_labels$to(dtype = torch_long())

  acc_mat <- matrix(0,maxiter,2)
  for(i in 1:maxiter){
    mlp_parameters <- lapply(mlp$parameters,function(a){
      t(as.matrix(a))})

    cur_elbo <- stm_torch_estep(counts(spe),
                                spe$int_cell,
                                rowData(spe)$gene_ints,
                                alphaPrior(spe),
                                betaPrior(spe),
                                K,
                                ncol(spe),
                                zero_gamma,
                                rand_gamma,
                                labels,
                                ndk(spe),
                                nwk(spe),
                                elbo_mat,
                                mlp_layers,
                                mlp_parameters,
                                cur_elbo,
                                num_threads,
                                burnin,
                                dummy_topic,
                                i-1) ##Added dummy topic

    #Calc N_mean
    #M STEP
    spe <- buildTheta(spe)
    spe <- buildPhi(spe)


    spe <- inferTopics(spe,1,100,verbal = FALSE,phi(spe),burnin = burnin)
    spe <- buildTheta(spe)


    ndk_prop <- torch::torch_tensor(theta(spe),device = device)
    ##Added Dummy Topic Edit
    if (dummy_topic){
      ndk_prop <- ndk_prop[,1:(K-1)]
    }
    ####
    train_dl <- mlp_data(ndk_prop,tensor_labels) %>%
      dataloader(1000,TRUE,num_workers = 0)
    if (balanced_class){
      class_prop <- as.numeric(table(labels + 1)/ncol(spe))^(-1) %>%
        torch_tensor() %>%
        torch_squeeze()

    } else{
      num_layer <- length(unique(labels))
      class_prop <- rep(1,num_layer) %>%
        torch_tensor() %>%
        torch_squeeze()

    }

    mlp <- fit_mlp(mlp,train_dl,lr,mlp_epoch,class_prop)
    if (i == maxiter){mlp <- fit_mlp(mlp,train_dl,lr,final_run_iter, class_prop)}
    cur_pred <- mlp(ndk_prop)
    ce_loss <- nnf_cross_entropy(cur_pred,tensor_labels)
    ce_loss <- as.numeric(ce_loss$cpu())
    #Calculate ELBO
    cur_elbo <- cur_elbo -ce_loss/ncol(counts(spe))
    dif <- abs(old_elbo-cur_elbo)


    elbo_mat[i,5] <- -ce_loss/ncol(counts(spe))
    if (dif <= thresh & i < maxiter){
      mlp <- fit_mlp(mlp,train_dl,1e-5,final_run_iter,class_prop)
      if (verbal){print("ELBO converged early")}
      break
    }
    old_elbo <- cur_elbo
    #print(cur_elbo)
    elbo_mat[i,4] <- cur_elbo
    cur_elbo <- 0
    ## Check Accuracy


    train_pred <- pred_mlp(mlp,theta(spe),device,dummy_topic)

    train_acc <- length(which(train_pred == as.numeric(labels+1)))*
      100/ncol(spe)

    ###
    if (!is.null(spe_val)){
      spe_val <- inferTopics(spe_val,1,100,verbal = FALSE,phi(spe),burnin = burnin)
      spe_val <- buildTheta(spe_val)


      val_pred <- pred_mlp(mlp,theta(spe_val),device,dummy_topic)

      val_acc <- length(which(val_pred == as.numeric(val_labels)))*
        100/ncol(spe_val)
    } else {
      val_acc <- NA
    }

    acc_mat[i,] <- c(train_acc,val_acc)
    if (verbal){ progress_bar(i,maxiter,elbo_mat[i,4],train_acc,val_acc)}

  }
  #TODO Find a nicer way to do the below segment
  return(mlp)

}


#' Burn-in Phase for SpaTM-S using Torch
#'
#' This function performs the burn-in phase for training SpaTM-S using a multi-layer perceptron (MLP) model.
#'
#' @param spe A `SpatialTopicExperiment` object.
#' @param labels A vector of labels for the cells.
#' @param num_threads Number of threads to use for parallel computation (default is 1).
#' @param maxiter Maximum number of iterations for the burn-in phase (default is 100).
#' @param verbal Logical indicating whether to display progress messages (default is `TRUE`).
#' @param zero_gamma Logical indicating whether to initialize the gamma parameters to zero (default is `TRUE`).
#' @param rand_gamma Logical indicating whether to randomize the gamma parameters (default is `FALSE`).
#' @param thresh Convergence threshold for stopping the burn-in phase (default is 0.00001).
#' @param lr Learning rate for optimization (default is 0.0001).
#' @param mlp A model definition for the multi-layer perceptron (MLP).
#' @param mlp_layers A vector specifying the number of units in each layer of the MLP.
#' @param mlp_epoch The number of epochs for training the MLP during the burn-in phase.
#' @param spe_val A validation set.
#' @param burnin Maximum number of iterations to run during the STM LDA step for each sample (default is 0)
#'
#' @seealso \link{STM}
#'
#' @import torch dplyr SingleCellExperiment
#' @export
STM_Torch_burnin <- function(spe,
                             labels,
                             num_threads = 1L,
                             maxiter = 100L, verbal = TRUE,
                             zero_gamma = TRUE,
                             rand_gamma = FALSE, thresh = 0.00001,
                             lr = 0.0001,
                             mlp,
                             mlp_layers,
                             mlp_epoch = 1000,
                             burnin = 1){

  old_elbo <- -9999999
  cur_elbo <- 0
  elbo_mat <- matrix(0,maxiter,5)
  K <- ncol(alphaPrior(spe))
  mlp_parameters <- lapply(mlp$parameters,function(a){
    t(as.matrix(a))})
  for(i in 1:maxiter){
    cur_elbo <- stm_torch_estep(counts(spe),
                                spe$int_cell,
                                rowData(spe)$gene_ints,
                                alphaPrior(spe),
                                betaPrior(spe),
                                K,
                                ncol(spe),
                                zero_gamma,
                                rand_gamma,
                                labels,
                                ndk(spe),
                                nwk(spe),
                                elbo_mat,
                                mlp_layers,
                                mlp_parameters,
                                cur_elbo,
                                num_threads,
                                burnin,
                                cur_iter = i-1)
    #spe <- GTM(spe,K,ncol(spe),num_threads,maxiter,F,F,T,burnin = 100)
  }
  return(spe)
}



## Full STM Wrapper
## TODO add class imbalance to GLM
## TODO Use MLP wrapper within this function for ease of use
## TODO Resolve Torch openMP bug!


#' Spatial Topic Modeling (STM) for Single-Cell Spatial Data
#'
#' This function performs spatial topic modeling on single-cell spatial transcriptomics data. It offers two modes:
#' 1. GLM-STM (default): A generalized linear model for topic modeling.
#' 2. Torch-STM: A deep learning-based STM model using a multi-layer perceptron (MLP).
#'
#' The function supports various parameters for training, including regularization, optimization settings, and the option for a validation set to prevent overfitting in the Torch-STM mode.
#'
#' @param spe A spatial experiment object (typically of class `SummarizedExperiment`) containing spatial transcriptomics data, including counts, metadata, and gene information.
#' @param label a character representing the column in colData(spe) with the class labels
#' @param num_threads Number of threads to use for parallel computation (default is 1).
#' @param maxiter Maximum number of iterations for model training (default is 100).
#' @param verbal Logical indicating whether to display progress messages (default is `TRUE`).
#' @param zero_gamma Logical indicating whether to initialize the gamma parameters to zero (default is `TRUE`).
#' @param rand_gamma Logical indicating whether to randomize the gamma parameters (default is `FALSE`).
#' @param thresh Convergence threshold for stopping the model training (default is 0.00001).
#' @param lr Learning rate for optimization (default is 0.0001).
#' @param mlp Optional, a model definition for the multi-layer perceptron (MLP) when running Torch-STM. If `NULL`, GLM-STM is run.
#' @param mlp_layers Optional, a vector specifying the number of units in each layer of the MLP.
#' @param mlp_epoch Optional, the number of epochs for training the MLP.
#' @param spe_val A validation set (typically of class `SummarizedExperiment`) used for monitoring overfitting in Torch-STM. Must be provided for Torch-STM.
#' @param device The device for training (e.g., "cpu" or "cuda") when using Torch-STM. If `NULL`, the default device is used.
#' @param balanced_class Logical indicating whether to balance the class distribution during model training (default is `TRUE`).
#' @param dummy_topic Logical indicating whether to add a dummy topic for regularization (default is `FALSE`).
#' @param burnin_lda Number of iterations for burn-in phase in Torch-STM (default is 5).
#' @param burnin Maximum number of iterations to run during the STM LDA step for each sample (default is 0)
#' @param final_run_iter Maximum number of iterations to optimize the MLP for in the last STM epoch (default is 100)
#'
#' @return The input spatial experiment object (`spe`) with updated metadata containing the STM model results:
#'   - `STM_Weights`: A matrix of learned topic weights (GLM-STM).
#'   - `STM_MLP`: The trained MLP model (Torch-STM).
#'
#' @import torch dplyr SingleCellExperiment methods
#' @importFrom S4Vectors metadata
#' @export
STM <- function(spe,
                label = NULL,
                num_threads = 1L,
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
                      burnin_lda = 5,
                      burnin = 100,
                      final_run_iter = 100,
                      nk_init = TRUE){
  if (!is(spe,'TopicExperiment')){
    stop('This function requires a TopicExperiment object to run correctly. Please create a Spatial or SingleCell variant before using the RTM function.')
  }
  if (is.null(label)){
    stop('You need to specify the column name in the TopicExperiment that contains the classes you want predicted.')
  } else if(!is.numeric(colData(spe)[,label])){
      if(!is.factor(colData(spe)[,label]) ){
        colData(spe)[,label] <- as.factor(colData(spe)[,label])
      }
      colData(spe)[,label] <- as.numeric(colData(spe)[,label]) - 1

  }

  if (!is.null(spe_val)){
    if(!is.numeric(colData(spe_val)[,label])){
      if(!is.factor(colData(spe_val)[,label]) ){
        colData(spe_val)[,label] <- factor(colData(spe_val)[,label],levels = colData(spe)[,label])
      }
      colData(spe_val)[,label] <- as.numeric(colData(spe_val)[,label]) - 1

    }
  }

  if (is.null(mlp)){
    if (verbal){
      message('Running GLM-STM')
    }
    K <- ncol(alphaPrior(spe))
    metadata(spe)[['STM_Weights']] <- train_stm(counts(spe),
                                                spe$int_cell,
                                                rowData(spe)$gene_ints,
                                                colData(spe)[,label],
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
                                                lr,
                                                nk_init)
    return(spe)
  }
  else{
     if(verbal){
       message('Running Torch-STM')
     }
    if (!requireNamespace("torch", quietly = TRUE)) {
      stop('You do not have torch installed. Please install it before trying to use the STM-Torch workflow')
    }
    if (burnin_lda > 0){
      spe <- STM_Torch_burnin(spe,
                               colData(spe)[,label],
                               num_threads,
                               burnin_lda,
                               verbal,
                               zero_gamma,
                               rand_gamma,
                               thresh,
                               lr,
                               mlp,
                               mlp_layers,
                               mlp_epoch = 0,
                               burnin = burnin)
    }

    metadata(spe)[['STM_MLP']] <- STM_Torch(spe,
                                            colData(spe)[,label],
                                            num_threads,
                                            maxiter,
                                            verbal,
                                            zero_gamma,
                                            rand_gamma,
                                            thresh,
                                            lr,
                                            mlp,
                                            mlp_layers,
                                            mlp_epoch,
                                            spe_val,
                                            balanced_class = balanced_class,
                                            burnin = burnin,
                                            final_run_iter = final_run_iter)
    return(spe)
  }

}




#' Testthat MLP forward pass helper function
#'
#' This helper function is only used in the testing suite to test the Rcpp forward pass
#'
#' The function supports various parameters for training, including regularization, optimization settings, and the option for a validation set to prevent overfitting in the Torch-STM mode.
#'
#' @param theta_mat input feature matrix
#' @param mlp_layers number of MLP hidden layers
#' @param mlp_parameters list of MLP parameters
#' @param dummy_topic boolean for whether to use dummy topic in model or not
#'
#' @return outputs from the Rcpp forward pass implementation
#'
#' @import torch
#' @export
test_mlp <- function(theta_mat,
                     mlp_layers,
                     mlp_parameters,
                     dummy_topic = FALSE){
  test_mlp_forward(theta_mat,
                   mlp_layers,
                   mlp_parameters,
                   dummy_topic)
}

## TODO remove script and move to reproducibility repo, not necessary in the package.
##RTM Helper Functions
library(torch)
library(Rcpp)
library(RcppArmadillo)
library(utils)
library(progress)
#Build mlp
rtm_build_mlp <- function(layers = 1,d_in,d_hidden,d_out = 1,device = torch_device('cpu')){
  if (!layers %in% 0:3){
    stop('Error: MLP is set to use 3 hidden layers at most!')
  }

  else if (layers == 0){
    mlp <- nn_sequential(
      nn_linear(d_in,d_out),
      nn_sigmoid()
    )
  } else if (layers == 1){
    mlp <- nn_sequential(
      nn_linear(d_in,d_hidden),
      nn_relu(),
      nn_linear(d_hidden,d_out),
      nn_sigmoid()
    )
  } else if (layers == 2){
    mlp <- nn_sequential(
      nn_linear(d_in,d_hidden[1]),
      nn_relu(),
      nn_linear(d_hidden[1],d_hidden[2]),
      nn_relu(),
      nn_linear(d_hidden[2],d_out),
      nn_sigmoid()
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
      nn_sigmoid()
    )
  }
  return(mlp)
}





# Define fitting function
fit_rtm_mlp <- function(mlp, train_dl, lr, max_epoch) {
  opt <- optim_adam(mlp$parameters, lr = lr,weight_decay = lr/2)
  loss_fn <- nn_bce_loss()

  for (i in 1:max_epoch) {
    coro::loop(
      for (batch in train_dl) {
        x_train <- batch$x
        y_train <- batch$y
        # Train_Pred
        y_pred <- mlp(x_train) %>%
          torch_squeeze()
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

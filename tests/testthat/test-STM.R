library(testthat)
library(SpaTM)
library(SpatialExperiment)
library(Matrix)
library(torch)
library(dplyr)

#TODO need to make sure pred_mlp works on GLM and NN
#TODO need to work on how the MLP is save if desired

# Helper function to create a mock SpatialTopicExperiment object
create_test_spatial_topic_experiment <- function() {
  counts <- matrix(rpois(100 * 50, lambda = 80), nrow = 100, ncol = 50)
  rownames(counts) <- paste0("Gene", 1:100)
  colnames(counts) <- paste0("Cell", 1:50)

  spatial_coords <- matrix(runif(100, min = 0, max = 25), ncol = 2)
  colnames(spatial_coords) <- c('array_row', 'array_col')

  metadata <- data.frame(
    sample_id = sample(as.character(1:10), 50, replace = TRUE),
    label = sample(LETTERS[1:5],50,replace = T)
  )

  metadata <- cbind(metadata, as.data.frame(spatial_coords))

  spe <- SpatialExperiment(
    assays = list(counts = as(counts, 'dgCMatrix')),
    spatialCoords = spatial_coords,
    colData = metadata
  )

  SpatialTopicExperiment(spe, K = 10)
}

# Example test cases for STM function

test_that("STM function returns expected results for valid input", {
  ste <- create_test_spatial_topic_experiment()
  result <- STM(ste,label = 'label', num_threads = 1, maxiter = 5, verbal = FALSE)
  # Define expected output based on the STM function's expected behavior
  expect_true(!is.null(metadata(result)[['STM_Weights']]), info = "STM_Weights should be present in metadata")
})

test_that("STM throws error for invalid input", {
  invalid_input <- "invalid_input"
  expect_error(STM(invalid_input, num_threads = 1, maxiter = 5, verbal = FALSE))
})

# Add tests for other helper functions in STM.R

test_that("build_mlp constructs the correct model", {
  mlp <- build_mlp(layers = 2, d_in = 100, d_hidden = c(64, 32), d_out = 10)

  # Check if the model has the correct structure
  expect_true(inherits(mlp, "nn_module"), info = "MLP should be an nn_module")
  expect_equal(length(mlp$parameters), 6, info = "MLP should have 6 sets of parameters")
})


test_that("fit_mlp trains the model correctly", {
  mlp <- build_mlp(layers = 1, d_in = 100, d_hidden = 50, d_out = 10)
  x <- torch_randn(1000, 100)
  y <- torch_randint(1, 11, size = c(1000), dtype = torch_long())
  train_dl <- mlp_data(x,y) %>%
    dataloader(1000,TRUE,num_workers = 0)
  trained_mlp <- fit_mlp(mlp, train_dl, lr = 0.001, max_epoch = 5, class_prop = torch_ones(10))

  # Check if the model has been trained
  expect_true(inherits(trained_mlp, "nn_module"), info = "Trained MLP should be an nn_module")
})

test_that("pred_mlp returns correct predictions", {
  mlp <- build_mlp(layers = 1, d_in = 100, d_hidden = 50, d_out = 10)
  x <- torch_randn(1000, 100)
  y <- torch_randint(1, 11, size = c(1000), dtype = torch_long())
  train_dl <- mlp_data(x,y) %>%
    dataloader(1000,TRUE,num_workers = 0)
  trained_mlp <- fit_mlp(mlp, train_dl, lr = 0.001, max_epoch = 5, class_prop = torch_ones(10))
  predictions <- pred_mlp(trained_mlp, x)

  # Check if the predictions are correct
  expect_true(is.numeric(predictions), info = "Predictions should be numeric")
  expect_equal(length(predictions), 1000, info = "Predictions should have the same length as input data")
})

test_that("STM_Torch returns expected results for valid input", {
  ste <- create_test_spatial_topic_experiment()
  K <- ncol(alphaPrior(ste))
  classes <- length(unique(ste$label))
  mlp <- build_mlp(layers = 1, d_in = K, d_hidden = 50, d_out = classes)
  ste$label <- as.numeric(as.factor(ste$label)) -1
  result <- STM_Torch(ste, ste$label,num_threads = 1, maxiter = 5, verbal = FALSE, mlp = mlp, mlp_layers = c(50), mlp_epoch = 1, spe_val = ste, device = torch_device('cpu'))

  # Define expected output based on the STM_Torch function's expected behavior
  expect_true(inherits(result, "nn_module"), info = "Result should be an nn_module")
})



test_that("STM_Torch_burnin runs without errors", {
  spe <- create_test_spatial_topic_experiment()
  labels <- sample(0:1, ncol(spe), replace = TRUE)
  mlp <- build_mlp(layers = 1, d_in = ncol(alphaPrior(spe)), d_hidden = 50, d_out = 2)
  spe_val <- create_test_spatial_topic_experiment()
  expect_silent(STM_Torch_burnin(spe, labels, num_threads = 1, maxiter = 5, verbal = FALSE, zero_gamma = TRUE, rand_gamma = FALSE, thresh = 0.00001, lr = 0.001, mlp = mlp, mlp_layers = c(50), mlp_epoch = 1, spe_val = spe_val))
})

test_that("Rcpp forward pass provides reasonable estimator for torch MLP",{
  mlp <- build_mlp(layers = 2, d_in = 10, d_hidden = c(64,32), d_out = 10)
  mlp$to(dtype  = torch_float64())
  x <- torch_randn(3000, 10,dtype = torch_float64())
  x_mat <- as.matrix(x)
  mlp_out <- as.matrix(mlp(x))
  mlp_parameters <- lapply(mlp$parameters,function(a){
    t(as.matrix(a$cpu()))})
  mlp_out_rcpp <- mlp_forward(x_mat,
                              2,
                              mlp_parameters)
  expect_true(cor(as.vector(mlp_out),as.vector(mlp_out_rcpp)) > 0.5)
})


test_that("progress_bar prints the correct output", {
  # Capture the console output
  output <- capture.output(progress_bar(100, 100, -1234.56, 85.0, 90.0))

  expected_output <- "[==================================================] 100% || Iter: 100 || Train: 85 || Val: 90 || ELBO: -1234.56 \r"
  # Compare the outputs
  expect_equal(output, expected_output)
})

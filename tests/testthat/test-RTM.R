suppressPackageStartupMessages({
  library(testthat)
  library(SpatialExperiment)
  library(SpaTM)
  library(bluster)
})
#TODO need tests for the helper functions

# Helper function to create a mock SpatialTopicExperiment object
create_test_spatial_topic_experiment <- function() {
  counts <- matrix(rpois(100 * 50, lambda = 80), nrow = 100, ncol = 50)
  rownames(counts) <- paste0("Gene", 1:100)
  colnames(counts) <- paste0("Cell", 1:50)

  spatial_coords <- matrix(runif(100, min = 0, max = 25),
                           ncol = 2)
  colnames(spatial_coords) <- c('array_row','array_col')
  # TODO need SpatialTopicExperiment to check for this and rename accordingly
  metadata <- data.frame(
    sample_id = sample(as.character(1:10), 50, replace = TRUE)
  )

  metadata <- cbind(metadata,
                    as.data.frame(spatial_coords))
  #TODO address this in SpTE object as well.
  #TODO Can check for array in colData and switch to spatialCoords if unavailable. Include message to user.
  spe <- SpatialExperiment(
    assays = list(counts = as(counts, 'dgCMatrix')),
    spatialCoords = spatial_coords,
    colData = metadata
  )

  SpatialTopicExperiment(spe,K = 10)
}

test_that("RTM function returns expected results for valid input", {
  ste <- create_test_spatial_topic_experiment()
  nbr_list <- get_nbrs(ste, samples = "sample_id", cell_ids = "int_cell", dist = 5)
  K <- ncol(alphaPrior(ste))
  result <- RTM(ste, K = K, nbr_list = nbr_list, num_threads = 1, maxiter = 5, verbal = FALSE)

  # Define expected output based on the RTM function's expected behavior
  expect_true(!is.null(metadata(result)[['RTM_weights']]), info = "RTM_weights should be present in metadata")
})

test_that("RTM function returns expected results for valid input using positive-example-only loss function", {
  ste <- create_test_spatial_topic_experiment()
  nbr_list <- get_nbrs(ste, samples = "sample_id", cell_ids = "int_cell", dist = 5,loss_fun = 0)
  K <- ncol(alphaPrior(ste))
  result <- RTM(ste, K = K, nbr_list = nbr_list, num_threads = 1, maxiter = 5, verbal = FALSE,loss_fun = 0)

  # Define expected output based on the RTM function's expected behavior
  expect_true(!is.null(metadata(result)[['RTM_weights']]), info = "RTM_weights should be present in metadata")
})

test_that("RTM function returns expected results for valid input using distance prediction loss_fun", {
  ste <- create_test_spatial_topic_experiment()
  nbr_list <- get_nbrs(ste, samples = "sample_id", cell_ids = "int_cell", dist = 5,loss_fun = 2)
  K <- ncol(alphaPrior(ste))
  result <- RTM(ste, K = K, nbr_list = nbr_list, num_threads = 1, maxiter = 5, verbal = FALSE,loss_fun = 2)

  # Define expected output based on the RTM function's expected behavior
  expect_true(!is.null(metadata(result)[['RTM_weights']]), info = "RTM_weights should be present in metadata")
})



test_that("RTM returns error if no neighbors included", {
  ste <- create_test_spatial_topic_experiment()
  #Building empty nbr_list
  nbr_list <- lapply(1:ncol(ste),function(a){
    matrix(0,nrow = 0,ncol= 2)
  })
  K <- ncol(alphaPrior(ste))
  expect_error(RTM(ste, K = K, nbr_list = nbr_list, num_threads = 1, maxiter = 5, verbal = FALSE))
})




test_that("RTM throws error for invalid input", {
  invalid_input <- "invalid_input"
  nbr_list <- lapply(1:20,function(a){
    matrix(0,nrow = 5,ncol= 2)
  })
  expect_error(RTM(invalid_input, K = 5, nbr_list = nbr_list, num_threads = 1, maxiter = 5, verbal = FALSE))
})

#### new tests



test_that("get_all_pred returns correct predictions", {
  ste <- create_test_spatial_topic_experiment()
  nbr_list <- get_nbrs(ste, samples = "sample_id", cell_ids = "int_cell", dist = 5)
  ste <- RTM(ste, K = ncol(alphaPrior(ste)), nbr_list = nbr_list, num_threads = 1, maxiter = 5, verbal = FALSE)
  result <- get_all_pred(ste, loss_fun = 1)

  # Define expected output based on the get_all_pred function's expected behavior
  expect_true(is.matrix(result), info = "Result should be a matrix")
  expect_true(ncol(result) == nrow(result) & ncol(result) == ncol(ste))
})


#TODO check an expected output as well
test_that("rtm_smooth returns correct labels", {
  ste <- create_test_spatial_topic_experiment()
  labels <- sample(c("Type1", "Type2"), 50, replace = TRUE)
  nbr_list <- get_nbrs(ste, samples = "sample_id", cell_ids = "int_cell", dist = 5)
  result <- rtm_smooth(ste, labels, nbr_list)
  # Define expected output based on the rtm_smooth function's expected behavior
  expect_equal(length(result), length(labels), info = "Result length should match labels length")
  result <- rtm_smooth(ste, labels, k = 5)
  expect_equal(length(result), length(labels), info = "Result length should match labels length")


  expect_error(rtm_smooth(ste, labels),info = 'Needs k or nbr_list argument')
  expect_error(rtm_smooth(ste, k = 5),info = 'Needs a vector of labels to update')
})


test_that("clust_optim returns correct clustering", {
  ste <- create_test_spatial_topic_experiment()
  reducedDim(ste,'PCA') <- matrix(rnorm(5*ncol(ste)),nrow = ncol(ste),ncol = 5)
  result <- clust_optim(ste, target_clusters = 2, max_iterations = 10, X = 'PCA', k = 4, alg = 'louvain')

  # Define expected output based on the clust_optim function's expected behavior
  expect_true(is.list(result), info = "Result should be a list")
  expect_true("clustering" %in% names(result), info = "Result should contain 'clustering'")
  expect_true("resolution" %in% names(result), info = "Result should contain 'resolution'")
})

test_that("minmax_norm normalizes matrix correctly", {
  adj_mat <- matrix(runif(100,10,200), nrow = 10, ncol = 10)
  result <- minmax_norm(adj_mat)

  # Define expected output based on the minmax_norm function's expected behavior
  expect_true(is.matrix(result), info = "Result should be a matrix")
  expect_true(all(result >= 0 & result <= 1), info = "Result values should be between 0 and 1")
})

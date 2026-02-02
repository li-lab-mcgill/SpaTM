suppressPackageStartupMessages({
  library(testthat)
  library(SpatialExperiment)
  library(SpaTM)
})

# Helper function to create a fixed mock SpatialExperiment object
mock_spe <- function() {
  colData <- DataFrame(
    array_row = c(1, 2, 2, 3, 3),
    array_col = c(1, 1, 2, 2, 3),
    sample_id = rep("sample1", 5),
    group = c("A", "A", "B", "B", "A"),
    cell_id = 1:5
  )

  spe <- SpatialExperiment(colData = colData)
  return(spe)
}

# Unit tests for get_nbrs

test_that("Basic functionality without grouping", {
  spe <- mock_spe()
  nbrs <- get_nbrs(spe, samples = "sample_id", cell_ids = "cell_id", dist = 2)
  expect_type(nbrs, "list")
  expect_true(length(nbrs) > 0)
})


## TODO address this test
test_that("Basic functionality with grouping", {
  spe <- mock_spe()
  nbrs <- get_nbrs(spe, samples = "sample_id", cell_ids = "cell_id", group_by = "group")
  expect_type(nbrs, "list")
  expect_true(do.call('sum',lapply(nbrs,nrow)) > 0)
})


test_that("Handles missing required columns", {
  spe <- SpatialExperiment(colData = DataFrame(sample_id = 1:5, cell_id = 1:5))
  expect_error(get_nbrs(spe, samples = "sample_id", cell_ids = "cell_id"))
})

test_that("Loss function 0 (no loss) works", {
  spe <- mock_spe()
  nbrs <- get_nbrs(spe, samples = "sample_id", cell_ids = "cell_id", loss_fun = 0)
  expect_type(nbrs, "list")
  all_nbrs <- do.call('rbind',nbrs)
  expect_true(sum(all_nbrs[,2]) == nrow(all_nbrs),
              info = "Should only contain positive examples.")

  expect_true(ncol(all_nbrs) == 2,info = 'Should only have nbr index and binary indicator.')
})

test_that("Loss function 1 (negative sampling) works", {
  spe <- mock_spe()
  nbrs <- get_nbrs(spe, samples = "sample_id", cell_ids = "cell_id", loss_fun = 1)
  expect_type(nbrs, "list")
  all_nbrs <- do.call('rbind',nbrs)
  expect_true(all(c(0,1) %in% all_nbrs[,2]),info = "Should contain + and - examples.")
  expect_true(ncol(all_nbrs) == 2,info = 'Should only have nbr index and binary indicator.')
})

test_that("Loss function 2 (Euclidean loss) works", {
  spe <- mock_spe()
  nbrs <- get_nbrs(spe, samples = "sample_id", cell_ids = "cell_id", loss_fun = 2)
  expect_type(nbrs, "list")
  all_nbrs <- do.call('rbind',nbrs)
  expect_true(all(c(0,1) %in% all_nbrs[,2]),info = "Should contain + and - examples.")
  expect_true(ncol(all_nbrs) == 3,info = 'Should nbr index binary indicator and distance.')
})

test_that("Correct neighbors are returned", {
  spe <- mock_spe()
  nbrs <- get_nbrs(spe, samples = "sample_id", cell_ids = "cell_id", dist = 1,loss_fun = 0)
  expected_neighbors <- list(
    matrix(c(2, 1), ncol = 2),
    matrix(c(1, 3, 1, 1), ncol = 2),
    matrix(c(2, 4, 1, 1), ncol = 2),
    matrix(c(3, 5, 1, 1), ncol = 2),
    matrix(c(4, 1), ncol = 2)
  )
  expected_neighbors <- lapply(expected_neighbors,
                               function(a){
                                 a[,1] <- a[,1]-1
                                 a})
  expect_equal(lapply(nbrs, unname), lapply(expected_neighbors, unname))
})



test_that("Get warning if no neighbors detected", {
  spe <- mock_spe()
  expect_message(get_nbrs(spe, samples = "sample_id", cell_ids = "cell_id", dist = 0),
                 'Warning. Your current parametrization returned 0 neighbors for all samples.')
})


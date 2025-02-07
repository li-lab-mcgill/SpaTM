library(testthat)
library(SpatialExperiment)
library(Matrix)
library(SpaTM)  # Ensure package is loaded

test_that("SpatialTopicExperiment class is correctly defined", {
  library(SpaTM)  # Ensure package is loaded

  expect_true(!is.null(getClassDef("SpatialTopicExperiment")))

  cls_def <- getClassDef("SpatialTopicExperiment")
  expect_true("SpatialExperiment" %in% names(cls_def@contains))
})

test_that("SpatialTopicExperiment constructor initializes correctly", {
  # Create a simple SpatialExperiment object
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:10)
  colnames(counts_matrix) <- paste0("Cell", 1:10)

  spe <- SpatialExperiment(assays = list(counts = counts_matrix))

  # Initialize SpatialTopicExperiment
  spte <- SpatialTopicExperiment(spe, K = 5, verbal = FALSE)

  # Check class
  expect_s4_class(spte, "SpatialTopicExperiment")

  # Check if required slots exist
  expect_true(!is.null(spte@theta))
  expect_true(!is.null(spte@phi))
  expect_true(!is.null(spte@alphaPrior))
  expect_true(!is.null(spte@betaPrior))
  expect_true(!is.null(spte@nwk))
  expect_true(!is.null(spte@ndk))

  # Check dimensions
  expect_equal(dim(spte@theta), c(10, 5))  # 10 cells, 5 topics
  expect_equal(dim(spte@phi), c(10, 5))    # 10 genes, 5 topics
})

test_that("SpatialTopicExperiment handles HVG filtering", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  spe <- SpatialExperiment(assays = list(counts = counts_matrix))

  spte <- SpatialTopicExperiment(spe, hvg = 5, verbal = FALSE)
  expect_true(nrow(spte) <= 5) # Should be reduced to 5 HVGs or less if below threshold
})

test_that("Guided mode requires labels", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  spe <- SpatialExperiment(assays = list(counts = counts_matrix))

  expect_error(SpatialTopicExperiment(spe, guided = TRUE),
               "labels argument must be defined if guided parameter is set to TRUE")
})

test_that("Guided mode assigns correct priors", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  coldata <- data.frame(Celltype = factor(rep(1:2, each = 5)))
  spe <- SpatialExperiment(assays = list(counts = counts_matrix), colData = coldata)

  spte <- SpatialTopicExperiment(spe, guided = TRUE, labels = "Celltype", verbal = FALSE)

  expect_equal(ncol(spte@alphaPrior), length(unique(coldata$Celltype)))  # Should match number of cell types
})

test_that("Default K is set when missing", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  spe <- SpatialExperiment(assays = list(counts = counts_matrix))

  expect_message(SpatialTopicExperiment(spe), "Setting topics to K = 10")
})

test_that("Initialized matrices have correct dimensions", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  spe <- SpatialExperiment(assays = list(counts = counts_matrix))
  spte <- SpatialTopicExperiment(spe, K = 5, verbal = FALSE)

  expect_equal(dim(spte@theta), c(10, 5))
  expect_equal(dim(spte@phi), c(10, 5))
  expect_equal(dim(spte@alphaPrior), c(10, 5))
  expect_equal(dim(spte@betaPrior), c(10, 5))
  expect_equal(dim(spte@nwk), c(10, 5))
  expect_equal(dim(spte@ndk), c(10, 5))
})

test_that("Row and column names are correctly assigned", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:10)
  colnames(counts_matrix) <- paste0("Cell", 1:10)
  spe <- SpatialExperiment(assays = list(counts = counts_matrix))

  spte <- SpatialTopicExperiment(spe, K = 5, verbal = FALSE)

  expect_equal(rownames(spte@phi), rownames(counts_matrix))
  expect_equal(rownames(spte@theta), colnames(counts_matrix))
})



test_that("Count matrix is converted to Sparse matrix if it is not already", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  spe <- SpatialExperiment(assays = list(counts = counts_matrix))
  expect_message(SpatialTopicExperiment(spe, K = 5),
                 regexp = "Converting counts to sparse matrix \\(dgCmatrix\\)",
                 fixed = FALSE)
  spe <- SpatialTopicExperiment(spe)
  expect_true(is(counts(spe),'dgCMatrix'))
})


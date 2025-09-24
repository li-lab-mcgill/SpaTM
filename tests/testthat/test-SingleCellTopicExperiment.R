library(testthat)
library(SingleCellExperiment)
library(Matrix)
library(SpaTM)  # Ensure package is loaded

# TODO proofread tests to ensure correctness


test_that("SingleCellTopicExperiment class is correctly defined", {
  library(SpaTM)  # Ensure package is loaded

  expect_true(!is.null(getClassDef("SingleCellTopicExperiment")))

  cls_def <- getClassDef("SingleCellTopicExperiment")
  expect_true("SingleCellExperiment" %in% names(cls_def@contains))
  expect_true("TopicExperiment" %in% names(cls_def@contains))
})

test_that("SingleCellTopicExperiment constructor initializes correctly", {
  # Create a small SingleCellExperiment object
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:10)
  colnames(counts_matrix) <- paste0("Cell", 1:10)
  counts_matrix <- as(counts_matrix,'dgCMatrix')
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

  # Initialize SingleCellTopicExperiment
  scte <- SingleCellTopicExperiment(sce, K = 5, verbal = FALSE)

  # Check class
  expect_s4_class(scte, "SingleCellTopicExperiment")

  # Check if required slots exist
  expect_true(!is.null(scte@theta))
  expect_true(!is.null(scte@phi))
  expect_true(!is.null(scte@alphaPrior))
  expect_true(!is.null(scte@betaPrior))
  expect_true(!is.null(scte@nwk))
  expect_true(!is.null(scte@ndk))

  # Check dimensions
  expect_equal(dim(scte@theta), c(10, 5))  # 10 cells, 5 topics
  expect_equal(dim(scte@phi), c(10, 5))    # 10 genes, 5 topics
})



test_that("SingleCellTopicExperiment handles HVG filtering", {
  counts_matrix <- matrix(rpois(5000, lambda = 100), nrow = 50, ncol = 100)
  counts_matrix <- as(counts_matrix,'dgCMatrix')
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

  scte <- SingleCellTopicExperiment(sce,K = 5, hvg = 5, verbal = FALSE)
  expect_true(nrow(scte) <= 5) # Should be reduced to 5 HVGs or less if below threshold
})

test_that("Guided mode requires labels", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  counts_matrix <- as(counts_matrix,'dgCMatrix')
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

  expect_error(SingleCellTopicExperiment(sce, guided = TRUE),
               "labels argument must be defined if guided parameter is set to TRUE")
})

test_that("Guided mode assigns correct priors", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  counts_matrix <- as(counts_matrix,'dgCMatrix')
  coldata <- data.frame(Celltype = factor(rep(1:2, each = 5)))
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix), colData = coldata)

  scte <- SingleCellTopicExperiment(sce, guided = TRUE, labels = "Celltype", verbal = FALSE)

  expect_equal(ncol(scte@alphaPrior), length(unique(coldata$Celltype)))  # Should match the number of cell types
})

test_that("Default K is set when missing", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  counts_matrix <- as(counts_matrix,'dgCMatrix')
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

  expect_message(SingleCellTopicExperiment(sce), "Setting topics to K = 10")
})

test_that("Initialized matrices have correct dimensions", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  counts_matrix <- as(counts_matrix,'dgCMatrix')
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
  scte <- SingleCellTopicExperiment(sce, K = 5, verbal = FALSE)

  expect_equal(dim(scte@theta), c(10, 5))
  expect_equal(dim(scte@phi), c(10, 5))
  expect_equal(dim(scte@alphaPrior), c(10, 5))
  expect_equal(dim(scte@betaPrior), c(10, 5))
  expect_equal(dim(scte@nwk), c(10, 5))
  expect_equal(dim(scte@ndk), c(10, 5))
})

test_that("Row and column names are correctly assigned", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  counts_matrix <- as(counts_matrix,'dgCMatrix')
  rownames(counts_matrix) <- paste0("Gene", 1:10)
  colnames(counts_matrix) <- paste0("Cell", 1:10)
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

  scte <- SingleCellTopicExperiment(sce, K = 5, verbal = FALSE)

  expect_equal(rownames(scte@phi), rownames(counts_matrix))
  expect_equal(rownames(scte@theta), colnames(counts_matrix))
})




test_that("Count matrix is converted to Sparse matrix if it is not already", {
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
  expect_message(sce <- SingleCellTopicExperiment(sce, K = 5),
                 regexp = "Converting counts to sparse matrix \\(dgCmatrix\\)",
                 fixed = FALSE)
  #sce <- SingleCellTopicExperiment(sce)
  expect_true(is(counts(sce),'dgCMatrix'))
})



test_that("Class Imbalance correction for alphaPrior is correctly handled",{
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  counts_matrix <- as(counts_matrix,'dgCMatrix')
  coldata <- data.frame(Celltype = factor(c(rep('A',8),rep('B',2))))
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix), colData = coldata)

  scte <- SingleCellTopicExperiment(sce, guided = TRUE, labels = "Celltype", verbal = FALSE, balanced = TRUE)

  expected_alpha_prior <- matrix(0.01, nrow = ncol(sce), ncol = length(unique(coldata$Celltype)))
  expected_alpha_prior[1:8,1] <- (ncol(sce) / 8)
  expected_alpha_prior[9:10,2] <- (ncol(sce) / 2)

  expect_true(all(abs(scte@alphaPrior - expected_alpha_prior) < 1e-6),
              info = "alphaPrior should be adjusted for class imbalance")
})

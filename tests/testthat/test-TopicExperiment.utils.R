library(testthat)
library(SpatialExperiment)
library(SpaTM)
# Mock TopicExperiment objects
test_spe <- SpatialExperiment(
  assays = list(counts = matrix(rpois(100 * 50, lambda = 10), nrow = 100, ncol = 50))
)

test_sce <- SingleCellExperiment(
  assays = list(counts = matrix(rpois(100 * 50, lambda = 10), nrow = 100, ncol = 50))
)

# Unit tests for TopicExperiment

test_that("TopicExperiments initialize correctly", {
  spte <- SpatialTopicExperiment(test_spe, K = 5)
  scte <- SingleCellTopicExperiment(test_sce, K = 5)
  expect_s4_class(spte, "SpatialTopicExperiment")
  expect_equal(ncol(theta(spte)), 5)  # Ensure K topics are set correctly
  expect_equal(ncol(phi(spte)), 5)

  expect_s4_class(scte, "SingleCellTopicExperiment")
  expect_equal(ncol(theta(scte)), 5)  # Ensure K topics are set correctly
  expect_equal(ncol(phi(scte)), 5)


})



# Test Show Method
test_that("show method prints correct output", {
  spte <- SpatialTopicExperiment(test_spe, K = 5)
  expect_output(print(spte), "TopicPriors")
  expect_output(print(spte), "SuffStats")
  expect_output(print(spte), "ExpTopics")
  scte <- SingleCellTopicExperiment(test_sce, K = 5)
  expect_output(print(scte), "TopicPriors")
  expect_output(print(scte), "SuffStats")
  expect_output(print(scte), "ExpTopics")
})

# Test Subsetting Method
test_that("subsetting maintains topic matrices", {
  spte <- SpatialTopicExperiment(test_spe, K = 5)
  scte <- SingleCellTopicExperiment(test_sce, K = 5)

  sub_spte <- spte[1:10, 1:20]
  sub_scte <- scte[1:10,1:20]
  expect_equal(nrow(sub_spte@betaPrior), 10)
  expect_equal(nrow(sub_spte@nwk), 10)
  expect_equal(nrow(sub_spte@phi), 10)
  expect_equal(nrow(sub_spte@alphaPrior), 20)
  expect_equal(nrow(sub_spte@ndk), 20)
  expect_equal(nrow(sub_spte@theta), 20)


  expect_equal(nrow(sub_scte@betaPrior), 10)
  expect_equal(nrow(sub_scte@nwk), 10)
  expect_equal(nrow(sub_scte@phi), 10)
  expect_equal(nrow(sub_scte@alphaPrior), 20)
  expect_equal(nrow(sub_scte@ndk), 20)
  expect_equal(nrow(sub_scte@theta), 20)
})

# Test Dimnames Setter
test_that("dimnames are updated correctly", {
  spte <- SpatialTopicExperiment(test_spe, K = 5)
  scte <- SingleCellTopicExperiment(test_sce, K = 5)

  new_row_names <- paste0("Gene_New", 1:nrow(spte))
  new_col_names <- paste0("Cell_New", 1:ncol(spte))
  dimnames(spte) <- list(new_row_names, new_col_names)
  dimnames(scte) <- list(new_row_names, new_col_names)

  expect_equal(rownames(spte@betaPrior), new_row_names)
  expect_equal(rownames(spte@nwk), new_row_names)
  expect_equal(rownames(spte@phi), new_row_names)
  expect_equal(rownames(spte@alphaPrior), new_col_names)
  expect_equal(rownames(spte@ndk), new_col_names)
  expect_equal(rownames(spte@theta), new_col_names)


  expect_equal(rownames(scte@betaPrior), new_row_names)
  expect_equal(rownames(scte@nwk), new_row_names)
  expect_equal(rownames(scte@phi), new_row_names)
  expect_equal(rownames(scte@alphaPrior), new_col_names)
  expect_equal(rownames(scte@ndk), new_col_names)
  expect_equal(rownames(scte@theta), new_col_names)
})


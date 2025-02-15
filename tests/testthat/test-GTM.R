library(testthat)
library(SingleCellExperiment)
library(Matrix)
library(SpaTM)


##Logical Function to check that outputs are 'nearly' identical
#accomodates float point imprecisions
near_match <- function(target, test, tol = 1e-8) {
  all(abs(target - test) < tol)
}

test_sce <- SingleCellExperiment(
  assays = list(counts = matrix(rpois(100 * 50, lambda = 80), nrow = 100, ncol = 50))
)
counts(test_sce) <- as(counts(test_sce),'dgCMatrix')
# Test 1: Check that ndk and nwk are no longer zero matrices after training
test_that("ndk and nwk are updated correctly after training", {
  scte <- SingleCellTopicExperiment(test_sce, K = 5)

  # Train the model using GTM
  scte_trained <- GTM(scte, K = 5, D = ncol(scte), num_threads = 1, maxiter = 5,verbal = F)

  expect_false(all(ndk(scte_trained) == 0), info = "ndk should not be a zero matrix")

  row_sums_ndk <- rowSums(ndk(scte_trained))
  colsums_counts <- Matrix::colSums(counts(scte_trained))
  expect_true(near_match(row_sums_ndk,colsums_counts), info = "Row sums of ndk should match row sums of counts matrix")

  # Test that nwk is no longer a zero matrix and sums match the total counts
  expect_false(all( nwk(scte_trained) == 0), info = "nwk should not be a zero matrix")

  total_nwk_sum <- sum(nwk(scte_trained))
  total_counts_sum <- sum(counts(scte_trained))
  expect_equal(total_nwk_sum, total_counts_sum, info = "Sum of nwk should match the total counts")


  scte_trained <- buildTheta(scte_trained)
  scte_trained <- buildPhi(scte_trained)
  expect_true(near_match(1,rowSums(theta(scte_trained))), info = "All row sums of theta should be equal to 1")
  expect_true(near_match(1,colSums(phi(scte_trained))), info = "All col sums of phi should be equal to 1")

})

# Test 2: Check that phi and theta have colSums and rowSums equal to 1
test_that("inferTopics should only update ndk and theta matrices with correct sums", {
  scte <- SingleCellTopicExperiment(test_sce, K = 5)
  phi <- matrix(runif(500), ncol = 5)  # Replace with actual phi matrix from a trained model
  phi <- sweep(phi,2,colSums(phi),FUN = "/")


  # Perform topic inference
  scte_inferred <- inferTopics(scte, num_threads = 1, maxiter = 5, verbal = FALSE, phi = phi)
  expect_true(all( nwk(scte_inferred) == 0), info = "nwk should be a zero matrix")
  expect_false(all( ndk(scte_inferred) == 0), info = "ndk should not be a zero matrix")


  row_sums_ndk <- rowSums(ndk(scte_inferred))
  colsums_counts <- Matrix::colSums(counts(scte_inferred))
  expect_true(near_match(row_sums_ndk,colsums_counts), info = "Row sums of ndk should match col sums of counts matrix")

  scte_inferred <- buildTheta(scte_inferred)
  expect_true(near_match(1,rowSums(theta(scte_inferred))), info = "All row sums of theta should be equal to 1")



})


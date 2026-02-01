library(testthat)
library(SingleCellExperiment)
library(Matrix)
library(SpaTM)

## Logical helper to check approximate equality
near_match <- function(target, test, tol = 1e-8) {
  all(abs(target - test) < tol)
}

set.seed(42)

test_sce_sgtm <- SingleCellExperiment(
  assays = list(counts = matrix(rpois(80 * 40, lambda = 60), nrow = 80, ncol = 40))
)
counts(test_sce_sgtm) <- as(counts(test_sce_sgtm),'dgCMatrix')

test_that("sGTM updates ndk and produces valid nwk", {
  scte <- SingleCellTopicExperiment(test_sce_sgtm, K = 4)

  scte_trained <- sGTM(scte, K = 4, D = ncol(scte), batch_size = 8,
                       num_threads = 1, maxiter = 3, verbal = FALSE,
                       lr = 0.5, shuffle = FALSE)

  expect_false(all(ndk(scte_trained) == 0), info = "ndk should not be a zero matrix")

  row_sums_ndk <- rowSums(ndk(scte_trained))
  colsums_counts <- Matrix::colSums(counts(scte_trained))
  expect_true(near_match(row_sums_ndk,colsums_counts), info = "Row sums of ndk should match col sums of counts matrix")

  expect_false(all(nwk(scte_trained) == 0), info = "nwk should not be a zero matrix")
  expect_true(all(is.finite(nwk(scte_trained))), info = "nwk should contain finite values")
})

test_that("sGTM handles full-batch setting", {
  scte <- SingleCellTopicExperiment(test_sce_sgtm, K = 3)

  scte_trained <- sGTM(scte, K = 3, D = ncol(scte), batch_size = ncol(scte),
                       num_threads = 1, maxiter = 2, verbal = FALSE,
                       lr = 0.3, shuffle = FALSE)

  row_sums_ndk <- rowSums(ndk(scte_trained))
  colsums_counts <- Matrix::colSums(counts(scte_trained))
  expect_true(near_match(row_sums_ndk,colsums_counts), info = "Full batch ndk sums should match counts")
})

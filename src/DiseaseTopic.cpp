#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

/*
 * Script for alpha prior helper function for disease-informed modelling
 */




//' @title Build Disease-informed alpha priors
//' @description updates the initialized label-guided alpha prior matrix to include disease information. We assume that most cells from healthy patients are healthy and that there is an equal chance of a cell from a diseased patient representing a healthy or diseased state.
//' @param a A numeric matrix representing the original alpha priors.
//' @param disease An unsigned integer vector indicating disease presence (1) or absence (0) for each row.
//' @return A numeric matrix with adjusted alpha values. Rows are modified based on the corresponding disease indicator:
//' - If `disease[cur] == 1`, the row is element-wise multiplied by {0.5, 0.5}.
//' - If `disease[cur] == 0`, the row is element-wise multiplied by {0.9, 0.1}.
//' @examples
//' \dontrun{
//' a <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
//' disease <- c(1, 0)
//' BuildDiseaseAlpha(a, disease)
//' }
//' @export
// [[Rcpp::export]]
arma::mat BuildDiseaseAlpha(const arma::mat& a, const arma::uvec& disease) {
  arma::mat result = a;

  arma::rowvec disease_vec = arma::rowvec(a.n_cols, arma::fill::ones) * 0.5;

  arma::rowvec ctrl_vec(a.n_cols);
  for (size_t i = 0; i < a.n_cols; ++i) {
    ctrl_vec[i] = (i % 2 == 0) ? 0.9 : 0.1;
  }

  for (size_t cur = 0; cur < disease.n_elem; ++cur) {
    if (disease[cur]) {
      result.row(cur) = a.row(cur) % disease_vec; // Element-wise multiplication
    } else {
      result.row(cur) = a.row(cur) % ctrl_vec; // Element-wise multiplication
    }
  }

  return result;
}

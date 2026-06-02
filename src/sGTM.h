#ifndef SGTM_H_
#define SGTM_H_

#include "CellMap.h"

#include <RcppArmadillo.h>
#include <cstdlib>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

void train_sgtm(arma::sp_mat& counts,
                     arma::vec& celltypes,
                     arma::vec& genes,
                     const arma::mat& alpha,
                     const arma::mat& beta,
                     int K,int D,
                     arma::mat& n_dk,
                     arma::mat& n_wk,
                     int batch_size,
                     int num_threads,
                     int maxiter,
                     bool verbal,
                     bool zero_gamma,
                     bool rand_gamma,
                     double thresh,
                     int burnin,
                     double lr,
                     bool shuffle);

#endif

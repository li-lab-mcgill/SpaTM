#ifndef GTM_H_
#define GTM_H_

#include "CellMap.h"

#include <RcppArmadillo.h>
#include <cstdlib>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;


//==========================
/*
 infer sufficient statistics
 */

//COMPLETE
void build_nwk(arma::mat& n_wk,
               std::unordered_map<int,Cell>& CellMap,
               int D, int K);
//COMPLETE
void build_ndk(arma::mat& n_dk,
               std::unordered_map<int,Cell>& CellMap, int D, int K);


// Complete
arma::mat get_theta(arma::mat& n_dk,
                    arma::mat& alpha);

//Complete
arma::mat get_phi(arma::mat& n_wk,
                  arma::mat& beta);
//==========================


double get_ELBO(std::unordered_map<int,Cell>& CellMap,
                const arma::mat& alpha,
                const arma::mat& beta,
                int M,
                int D,
                int K,
                arma::mat& n_dk,
                arma::mat& n_wk,
                int C);

//COMPLETE
void run_epoch(std::unordered_map<int,Cell>& CellMap,
               const arma::mat& alpha,
               const arma::mat& beta,
               int K,int D,
               arma::mat& n_dk,
               arma::mat& n_wk,
               int num_threads);
//COMPLETE
void train_gtm(arma::sp_mat& counts,
                     arma::vec& celltypes,
                     arma::vec& genes,
                     const arma::mat& alpha,
                     const arma::mat& beta,
                     int K,int D,
                     arma::mat& n_dk,
                     arma::mat& n_wk,
                     int num_threads,
                     int maxiter,
                     bool verbal,
                     bool zero_gamma,
                     bool rand_gamma,
                     double thresh);




// Prediction
void predict_epoch(std::unordered_map<int,Cell>& CellMap,const double& alpha,
                   int K,int D,int M,arma::mat& n_dk,
                   const arma::mat& phi, int num_threads);

void infer_topics_cpp(arma::sp_mat& counts,
                       arma::vec& celltypes,
                       arma::vec& genes,
                       double& alpha,
                       int K,int D,int M,arma::mat& n_dk, const arma::mat& phi,
                       int num_threads,int maxiter,
                       bool verbal);

#endif

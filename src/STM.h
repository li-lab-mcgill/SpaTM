#ifndef STM_H_
#define STM_H_
#include "CellMap.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using namespace std;


// Supervised Cell Map
class STM_Cell: public Cell {
  public:
    int cell_label;
    
    STM_Cell();
    
    STM_Cell(int cell_id,arma::mat cell_mtx,
         const arma::mat& alpha,
         const arma::mat& beta,
         int K,
         bool zero_gamma,
         bool rand_gamma,
         int label);
    
    STM_Cell(int cell_id,arma::mat cell_mtx,
         int K,
         int label);
    
    void print();
};


std::unordered_map<int,STM_Cell> build_STM_Cell_Map(arma::sp_mat& counts,
                                            arma::vec& celltypes, 
                                            arma::vec& genes,
                                            const arma::mat& alpha,
                                            const arma::mat& beta,
                                            int D, int K,
                                            bool zero_gamma,
                                            bool rand_gamma,
                                            arma::vec& labels);


// Suff Stat Helpers
void build_nwk_stm(arma::mat& n_wk,
                   std::unordered_map<int,STM_Cell>& STM_CellMap,
                   int D, int K);

void build_ndk_stm(arma::mat& n_dk,
                   std::unordered_map<int,STM_Cell>& STM_CellMap,
                   int D, int K);


double get_stm_ELBO(std::unordered_map<int,STM_Cell>& STM_CellMap,
                    const arma::mat& alpha,
                    const arma::mat& beta,
                    int M,
                    int D,
                    int K,
                    arma::mat& n_dk,
                    arma::mat& n_wk,
                    int C,
                    arma::mat& elbo,
                    int cur_iter);

arma::rowvec get_label_prob(arma::rowvec ndk_row,
                            const arma::rowvec& gamma_weight,
                            arma::mat& model_weights,
                            int gene_count,int K, int label);
#endif /* STM_H_ */

#ifndef RTM_H_
#define RTM_H_

#include "CellMap.h"



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using namespace std;


// Supervised Cell Map
class RTM_Cell: public Cell {
  public:
    arma::uvec cell_nbrs;
    arma::vec bool_nbrs;
    arma::vec nbr_dist;
    RTM_Cell();
    
    RTM_Cell(int cell_id,arma::mat cell_mtx,
         const arma::mat& alpha,
         const arma::mat& beta,
         int K,
         const arma::uvec& nbrs,
         bool zero_gamma,
         bool rand_gamma);
    
    RTM_Cell(int cell_id,arma::mat cell_mtx,
             const arma::mat& alpha,
             const arma::mat& beta,
             int K,
             const arma::uvec& nbrs,
             const arma::vec& nbr_bool,
             bool zero_gamma,
             bool rand_gamma);
    
    RTM_Cell(int cell_id,arma::mat cell_mtx,
             const arma::mat& alpha,
             const arma::mat& beta,
             int K,
             const arma::uvec& nbrs,
             const arma::vec& nbr_bool,
             const arma::vec& nbr_dist,
             bool zero_gamma,
             bool rand_gamma);
    
    RTM_Cell(int cell_id,arma::mat cell_mtx,
         int K,
         const arma::uvec& nbrs);
    
    RTM_Cell(int cell_id,arma::mat cell_mtx,
             int K,
             const arma::uvec& nbrs,
             const arma::vec& nbr_bool);
    
    void print();
};


std::unordered_map<int,RTM_Cell> build_RTM_Cell_Map(arma::sp_mat& counts,
                                            arma::vec& celltypes, 
                                            arma::vec& genes,
                                            const arma::mat& alpha,
                                            const arma::mat& beta,
                                            int D, int K,
                                            bool zero_gamma,
                                            bool rand_gamma,
                                            Rcpp::List& nbr_list,
                                            int loss_fun = 0);


// Suff Stat Helpers
void build_nwk_RTM(arma::mat& n_wk,
                   std::unordered_map<int,RTM_Cell>& RTM_CellMap,
                   int D, int K);

void build_ndk_RTM(arma::mat& n_dk,
                   std::unordered_map<int,RTM_Cell>& RTM_CellMap,
                   int D, int K);


double get_RTM_ELBO(std::unordered_map<int,RTM_Cell>& RTM_CellMap,
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


arma::rowvec get_nbr_prob(arma::rowvec ndk_row,
                          const arma::rowvec& gamma_weight,
                          arma::vec& model_weights,
                          int model_bias,
                          arma::mat& nbr_ndk,
                          int gene_count,int K,
                          int loss_fun,
                          arma::vec& nbr_bool,
                          arma::vec& nbr_dist);


#endif /* RTM_H_ */

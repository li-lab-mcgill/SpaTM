#ifndef CELLMAP_H_
#define CELLMAP_H_


#include <RcppArmadillo.h>
#include <cstdlib>
#include <omp.h>


class Cell {
public:
  int cell_id;
  //arma::rowvec ndk;
  arma::mat cell_mtx;
  arma::mat cell_gamma;
  Cell();

  Cell(int cell_id,arma::mat cell_mtx,
       const arma::mat& alpha,
       const arma::mat& beta,
       int K,
       bool zero_gamma,
       bool rand_gamma);

  Cell(int cell_id,arma::mat cell_mtx,
       int K);
  //~Cell();
  void print();


};




std::unordered_map<int,Cell> build_Cell_Map(arma::sp_mat& counts,
                                            arma::vec& celltypes,
                                            arma::vec& genes,
                                            const arma::mat& alpha,
                                            const arma::mat& beta,
                                            int D, int K,
                                            bool zero_gamma,
                                            bool rand_gamma,
                                            bool verbal,
                                            int num_threads);

std::unordered_map<int,Cell> build_Predict_Cell_Map(arma::sp_mat& counts,
                                                    arma::vec& celltypes,
                                                    arma::vec& genes,
                                                    int D, int K);

void test_cell_map(arma::sp_mat& counts,
                   arma::vec& celltypes,
                   arma::vec& genes,
                   const arma::mat& alpha,
                   const arma::mat& beta,
                   int D, int K,
                   bool zero_gamma,
                   bool rand_gamma,
                   int cellid);




#endif /* CELLMAP_H_ */


#include "CellMap.h"

#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

/*
 * Script containing helper functions for generating list of Cells
 */


Cell::Cell(){
  cell_id = 0;
  cell_mtx = arma::mat();
  cell_gamma = arma::mat();
}

Cell::Cell(int id,arma::mat cellmtx,
           const arma::mat& alpha, const arma::mat& beta,int K,
           bool zero_gamma = false,
           bool rand_gamma = true){
  cell_id = id;
  cell_mtx = cellmtx;
  int tokens = cell_mtx.n_rows;
  cell_gamma.resize(tokens,K);
  arma::rowvec beta_sum = arma::sum(beta,0);
  int gene,cell;
  arma::rowvec gamma_k = arma::zeros<arma::rowvec>(K);

  //Initialize posteriors
  for (int i = 0; i < tokens; i++){
    if (zero_gamma){
      gamma_k = arma::zeros<arma::rowvec>(K);
    }
    else if (rand_gamma){
      //gamma_k = arma::zeros<arma::rowvec>(K);
      gamma_k = arma::randu<arma::rowvec>(K);
      gamma_k = gamma_k/sum(gamma_k);
    }
    else{
      cell = cell_mtx(i,1)-1;
      gene = cell_mtx(i,2)-1;
      gamma_k = alpha.row(cell) % beta.row(gene) / beta_sum;
      gamma_k = gamma_k/sum(gamma_k);
    }


    //sequential par
    cell_gamma.row(i) = gamma_k;
  }
}


Cell::Cell(int id,arma::mat cellmtx,
           int K){
  cell_id = id;
  cell_mtx = cellmtx;
  int tokens = cell_mtx.n_rows;
  cell_gamma.resize(tokens,K);
  // arma::rowvec gamma_k = arma::zeros<arma::rowvec>(k);
  double gamma_init = 1.0/K;
  // Rcout << gamma_init << std::endl;
  // Rcout << 1.0/k << std::endl;
  // Rcout << k << std::endl;
  cell_gamma.fill(gamma_init);
  //Initialize posteriors
  // for (int i = 0; i < tokens; i++){
  //   gamma_k =  arma::rowvec(k,(alpha * beta / G_beta));
  //   gamma_k = gamma_k/sum(gamma_k);
  //   //sequential par
  //   cell_gamma.row(i) = gamma_k;
  // }
}

void Cell::print() {
  Rcout << "cell id: " << cell_id << std::endl;
  Rcout << "mtx: \n" << cell_mtx << std::endl;
  Rcout << "gamma: " << cell_gamma << std::endl;
}




// Cell::~Cell(){
//   //cell_id.clear();
//   cell_mtx.clear();
//   cell_gamma.clear();
// }
//build Cell hash list
std::unordered_map<int,Cell> build_Cell_Map(arma::sp_mat& counts,
                                            arma::vec& celltypes,
                                            arma::vec& genes,
                                            const arma::mat& alpha,
                                            const arma::mat& beta,
                                            int D, int K,
                                            bool zero_gamma,
                                            bool rand_gamma = true,
                                            bool verbal = false){
  std::unordered_map<int,Cell> CellMap;
  //int cell_start,cell_end;
  for (int i = 0; i < D; i++){
    arma::uvec geneidx = arma::find(counts.col(i));
    arma::mat mtx(geneidx.n_elem,3);

    mtx.col(0) = arma::nonzeros(counts.col(i));//counts(nonzero_counts);
    mtx.col(1).fill(celltypes(i));
    mtx.col(2) = genes(geneidx);
    if (verbal){
      if ( i == 10000){
        Rcout << i << " Cells Processed" << std::endl;
      }
      else if ( i == 50000){
        Rcout << i << " Cells Processed" << std::endl;
      }
      else if (i == 100000){
        Rcout << i << " Cells Processed" << std::endl;
      }
      else if (i == 150000){
        Rcout << i << " Cells Processed" << std::endl;
      }
    }
    CellMap[i+1] = Cell(i+1,mtx,alpha,beta,K,
                        zero_gamma,
                        rand_gamma);
  }

  return CellMap;
}

std::unordered_map<int,Cell> build_Predict_Cell_Map(arma::sp_mat& counts,
                                                    arma::vec& celltypes,
                                                    arma::vec& genes,
                                                    int D, int K){
  std::unordered_map<int,Cell> CellMap;

  for (int i = 0; i < D; i++){
    arma::uvec geneidx = arma::find(counts.col(i));
    arma::mat mtx(geneidx.n_elem,3);

    mtx.col(0) = arma::nonzeros(counts.col(i));//counts(nonzero_counts);
    mtx.col(1).fill(celltypes(i));
    mtx.col(2) = genes(geneidx);
    CellMap[i+1] = Cell(i+1,mtx,K);
  }

  return CellMap;
}

// [[Rcpp::export]]
void test_cell_map(arma::sp_mat& counts,
                   arma::vec& celltypes,
                   arma::vec& genes,
                   const arma::mat& alpha,
                   const arma::mat& beta,
                   int D, int K,
                   bool zero_gamma = false,
                   bool rand_gamma = true,
                   int cellid = 1){
  std::unordered_map<int,Cell> temp = build_Cell_Map(counts,
                                                     celltypes,
                                                     genes,
                                                     alpha,
                                                     beta,
                                                     D,K,
                                                     zero_gamma,
                                                     rand_gamma);
  temp[cellid].print();

  Rcout << "function ran" << std::endl;
}





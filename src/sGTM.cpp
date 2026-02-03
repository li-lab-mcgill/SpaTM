#include "CellMap.h"
#include "sGTM.h"
#include <omp.h>
#include <algorithm>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

static void build_nwk_batch(arma::mat& n_wk_batch,
                     std::unordered_map<int,Cell>& CellMap,
                     const arma::uvec& batch_ids){
  int tokens;
  n_wk_batch.zeros();
  for (arma::uword b = 0; b < batch_ids.n_elem; b++){
    int cell_id = batch_ids(b);
    tokens = CellMap[cell_id].cell_mtx.n_rows;
    for (int i = 0;i < tokens; i++){
      int gene = CellMap[cell_id].cell_mtx(i,2)-1;
      n_wk_batch.row(gene) += CellMap[cell_id].cell_gamma.row(i) * CellMap[cell_id].cell_mtx(i,0);
    }
  }
}

static void run_epoch_batch(std::unordered_map<int,Cell>& CellMap,const arma::mat& alpha,
               const arma::mat& beta,int K,arma::mat& n_dk,arma::mat& n_wk,
               const arma::uvec& batch_ids,
               int num_threads = 1,
               int burnin = 1){
  omp_set_num_threads(num_threads);

  const double eps = 1e-12;
  arma::rowvec beta_sum = arma::sum(beta,0);
  arma::rowvec nwk_sum = arma::sum(n_wk,0);
  #pragma omp parallel for
  for (arma::uword b = 0; b < batch_ids.n_elem; b++){
    int cell_id = batch_ids(b);
    arma::mat cur_gamma = CellMap[cell_id].cell_gamma;
    arma::mat m = CellMap[cell_id].cell_mtx;
    int ndk_id = cell_id - 1;
    int tokens = m.n_rows;
    double diff = 1;
    int max_update = burnin;
    int u = 0;
    arma::rowvec cur_ndk = n_dk.row(ndk_id);
    while (diff > 0.01 && u < max_update){
      u++;
      for (int token = 0; token < tokens; token++){
        int counts = m(token,0);
        int gene = m(token,2)-1;
        arma::rowvec cur_counts = cur_gamma.row(token)*counts;
        arma::rowvec denom = beta_sum + nwk_sum - cur_counts;
        denom = arma::clamp(denom, eps, arma::datum::inf);
        arma::rowvec gamma_k = (alpha.row(ndk_id) + n_dk.row(ndk_id) - cur_counts) %
          (beta.row(gene)+n_wk.row(gene)- cur_counts) /
            denom;
        double gamma_sum = arma::accu(gamma_k);
        if (!gamma_k.is_finite() || gamma_sum <= eps || arma::any(gamma_k < 0)){
          gamma_k = (alpha.row(ndk_id) + n_dk.row(ndk_id)) % (beta.row(gene)+
            n_wk.row(gene)) /
              (beta_sum+nwk_sum);
          gamma_sum = arma::accu(gamma_k);
          if (!gamma_k.is_finite() || gamma_sum <= eps){
            gamma_k.fill(1.0 / static_cast<double>(K));
          }
        }
        gamma_k = gamma_k/arma::accu(gamma_k);
        cur_gamma.row(token) = gamma_k;
      }

      n_dk.row(ndk_id).zeros();
      for (int token = 0; token < tokens; token++){
        n_dk.row(ndk_id) += cur_gamma.row(token)*m(token,0);
      }

      diff = abs(arma::accu(cur_ndk - n_dk.row(ndk_id)));
      cur_ndk = n_dk.row(ndk_id);
    }
    CellMap[cell_id].cell_gamma = cur_gamma;
  }
}

// [[Rcpp::export]]
void train_sgtm(arma::sp_mat& counts,
                     arma::vec& celltypes,
                     arma::vec& genes,
                     const arma::mat& alpha, const arma::mat& beta,
                     int K,int D,arma::mat& n_dk,arma::mat& n_wk,
                     int batch_size = 1024,
                     int num_threads = 1,int maxiter = 100,
                     bool verbal = true,
                     bool zero_gamma = false,
                     bool rand_gamma = true,
                     double thresh = 0.00001,
                     int burnin = 1,
                     double lr = 0.1,
                     bool shuffle = true){
  if (batch_size <= 0){
    batch_size = D;
  }
  if (batch_size > D){
    batch_size = D;
  }

  int M = genes.n_elem;
  int num_batches = (D + batch_size - 1) / batch_size;
  double prog = 0;
  int prog_width = 50;
  int pos = 0;
  arma::mat old_nwk = n_wk;
  const double tau0 = 1.0;
  const double kappa = 0.7;
  double step = 0.0;

  for (int iter = 0; iter < maxiter; iter++){
    arma::uvec order = shuffle ? (arma::randperm(D) + 1) : arma::regspace<arma::uvec>(1,D);

    for (int b = 0; b < num_batches; b++){
      arma::uword start = static_cast<arma::uword>(b * batch_size);
      arma::uword end = static_cast<arma::uword>(std::min(D, static_cast<int>(start) + batch_size));
      arma::uvec batch_ids = order.subvec(start, end - 1);

      std::unordered_map<int,Cell> CellMap = build_Cell_Map_batch(counts,
                                                                  celltypes,
                                                                  genes,
                                                                  alpha,
                                                                  beta,
                                                                  batch_ids,
                                                                  K,
                                                                  zero_gamma,
                                                                  rand_gamma,
                                                                  num_threads);

      run_epoch_batch(CellMap,alpha,beta,K,n_dk,n_wk,
                      batch_ids,num_threads,burnin);

      arma::mat batch_nwk = arma::zeros<arma::mat>(M,K);
      build_nwk_batch(batch_nwk,CellMap,batch_ids);

      double scale = static_cast<double>(D) / static_cast<double>(batch_ids.n_elem);
      double lr_t = lr * std::pow(tau0 + step, -kappa);
      lr_t = std::min(1.0, std::max(lr_t, 1e-6));
      n_wk = (1.0 - lr_t) * n_wk + lr_t * (scale * batch_nwk);
      step += 1.0;
      CellMap.clear();

      if (verbal){
        double batch_prog = (static_cast<double>(iter) +
          (static_cast<double>(b + 1) / static_cast<double>(num_batches))) /
          static_cast<double>(maxiter);
        Rcout.flush();
        prog = batch_prog;
        Rcout << "[";
        pos = int(prog_width * prog);
        for(int p = 0; p < prog_width; p++){
          if (p < pos){ Rcout << "=" ;}
          else { Rcout << " " ;}
        }
        double diff = arma::max(arma::max(arma::abs(n_wk - old_nwk)));
        Rcout << "] " << int(prog * 100.0) << "% || Iter: " << iter
              << " || Batches: " << (b + 1) << "/" << num_batches
              << " || Delta: " << diff <<" \r";
      }
    }

    double diff = arma::max(arma::max(arma::abs(n_wk - old_nwk)));
    if (verbal){
      prog = (static_cast<double>(iter + 1)) / static_cast<double>(maxiter);
    }

    if (diff < thresh){
      if (verbal){
        Rcout.flush();
        Rcout << "[";
        for(int p = 0; p < prog_width; p++){
          Rcout << "=" ;
        }
          Rcout << "] " << "100% || Iter: " << iter
            << " || Batches: " << num_batches << "/" << num_batches
            << " || Delta: " << diff << std::endl;
        Rcout << "Model Converged at iteration : " << iter << std::endl;
      }
      return;
    }

    old_nwk = n_wk;
  }

  if (verbal){
    Rcout.flush();
    Rcout << "[";
    for(int p = 0; p < prog_width; p++){
      Rcout << "=" ;
    }
        Rcout << "] " << "100% || Iter: " << maxiter
          << " || Batches: " << num_batches << "/" << num_batches
          << " || Delta: " << arma::max(arma::max(arma::abs(n_wk - old_nwk))) << std::endl;
    Rcout << "Max Iteration Reached" << std::endl;
  }
  return;

}

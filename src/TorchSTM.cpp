#include "CellMap.h"
#include "STM.h"
#include <cstdlib>
#include <omp.h>
#include <fstream>

using namespace Rcpp;
using namespace std;
// Script for Torch-STM Implementation
//TODO
//use previous n_dk for updates
//update after parsing all tokens
//compare to original
//rowvec old_ndk = ndk.row(ndk_id);

// Progress bar for Torch-based model
// [[Rcpp::export]]
void progress_bar(double iter,int maxiter,double elbo, double train_acc,double val_acc){
  int prog_width = 50;
  int pos = 0;
  Rcout.flush();
  double prog = iter/maxiter;
  Rcout << "[";
  pos = int(prog_width * prog);
  for(int p = 0; p < prog_width; p++){
    if (p < pos){ Rcout << "=" ;}
    else { Rcout << " " ;}
  }
  Rcout << "] " << int(prog * 100.0) << "% || Iter: " << iter << " || Train: "<< train_acc << " || Val: "<< val_acc  << " || ELBO: "<< elbo <<" \r";

}




//forward MLP
// [[Rcpp::export]]
arma::mat mlp_forward(const arma::mat& X,
                      int layers,
                      Rcpp::List& weights,
                      bool dummy_topic = false) {


  arma::mat cur_val = X;
  //Added dummy topic feature
  if (dummy_topic){
    cur_val.shed_col(cur_val.n_cols - 1);
  }
  ////
  int i  = 0;
  // Calculate hidden layer output using ReLU activation
  arma::mat weight;
  arma::rowvec bias;
  while (i < weights.length()){
    weight = Rcpp::as<arma::mat>(weights[i]);
    bias = Rcpp::as<arma::rowvec>(weights[i+1]);

    cur_val = cur_val*weight;
    cur_val.each_row() += bias;
    //Rcout << size(cur_val);
    if (layers > 0 && i != layers){cur_val = cur_val % (cur_val > 0);}
    i += 2;
  }

  // Calculate output
  // if (cur_val.has_nan()){
  //   Rcout << cur_val << std::endl;
  //   stop("Label Out has Nans");
  // }
  arma::mat output = arma::exp(cur_val);
  // if (output.has_nan()){
  //   Rcout << output << std::endl;
  //   stop("Exp-Label Out has NaNs");
  // }
  arma::vec sum_exp = arma::sum(output,1);
  // if (sum_exp.has_nan()){
  //   Rcout << sum_exp << std::endl;
  //   stop("Sum Exp-Label Out has NaNs");
  // }
  output = output.each_col()/sum_exp;
  // if (output.has_nan()){
  //   Rcout << output << std::endl;
  //   stop("Softmax Out has NaNs");
  //   Rcout << sum_exp << std::endl;
  // }
  return output;
}


// [[Rcpp::export]]
arma::rowvec get_label_prob(arma::rowvec ndk_row,
                            const arma::rowvec& gamma_weight,
                            int gene_count,int K,
                            int label,int layers,
                            Rcpp::List& weights,
                            bool dummy_topic = false){

  arma::mat label_prob = arma::zeros<arma::mat>(K,K);


  ndk_row -= gamma_weight*gene_count; //remove cur value


  label_prob.each_row() = ndk_row;

  label_prob.diag() += gene_count;

  label_prob.each_col() /= arma::sum(label_prob,0).t();

  //////

  arma::mat mlp_out = mlp_forward(label_prob,
                                  layers,
                                  weights,
                                  dummy_topic);
  arma::rowvec final_prob = mlp_out.col(label).t();
  // concatenate 1/K to final_prob if dummy topic
  // Added dummy topic feature
  if (dummy_topic){
    final_prob = arma::join_horiz(final_prob,arma::rowvec(1/K));
  }
  ///
  // if (final_prob.has_nan()){
  //   Rcout << final_prob << std::endl;
  //   Rcout << std::endl;
  //   Rcout << mlp_out << std::endl;
  //   stop("Final Prob has NaNs (get_label_prob)");
  // }
  // get weights
  return final_prob;
}

//
// ////////////////////
// // Epoch train
//
// [[Rcpp::export]]
double stm_torch_estep(arma::sp_mat& spe_counts,
                     arma::vec& celltypes,
                     arma::vec& genes,
                     const arma::mat& alpha,
                     const arma::mat& beta,
                     int K, int D,
                     bool zero_gamma,
                     bool rand_gamma,
                     arma::vec& labels,
                     arma::mat& n_dk,
                     arma::mat& n_wk,
                     arma::mat& elbo_mat,
                     int layers,
                     Rcpp::List& weights,
                     double& elbo,
                     int num_threads = 1,
                     int iter = 0,
                     bool dummy_topic = false){
  omp_set_num_threads(num_threads);
  int C = arma::accu(spe_counts);
  int M = genes.n_elem;

  std::unordered_map<int,STM_Cell> STM_CellMap = build_STM_Cell_Map(spe_counts,
                                                                    celltypes,
                                                                    genes,
                                                                    alpha,
                                                                    beta,
                                                                    D,K,
                                                                    zero_gamma,
                                                                    rand_gamma,
                                                                    labels);

  arma::rowvec beta_sum = arma::sum(beta,0);
  arma::rowvec nwk_sum = arma::sum(n_wk,0);
  #pragma omp parallel for private(weights) shared(STM_CellMap, n_dk, n_wk, beta, beta_sum, nwk_sum)
  for (int i = 1; i <= D; i++){
    arma::rowvec gamma_k = arma::zeros<arma::rowvec>(K);
    arma::rowvec cur_counts = arma::zeros<arma::rowvec>(K);
    arma::rowvec cur_ndk = arma::zeros<arma::rowvec>(K);
    arma::rowvec label_post = arma::zeros<arma::rowvec>(K);
    double d_diff = 999;
    int iter = 100;

    arma::mat cur_gamma = STM_CellMap[i].cell_gamma; //TODO
    arma::mat m = STM_CellMap[i].cell_mtx; //TODO

    while((d_diff > 0.01) & (iter > 0)){ //consider OR statement
      iter--;

      int ndk_id = i-1; //To accomodate row indices for matrix
      int tokens = m.n_rows;
      cur_ndk = n_dk.row(ndk_id);

      for (int token = 0; token < tokens; token++){

        int counts = m(token,0);
        int gene = m(token,2)-1;
        cur_counts = cur_gamma.row(token)*counts;


        label_post = get_label_prob(cur_ndk,cur_gamma.row(token),
                                    counts,K,
                                    STM_CellMap[i].cell_label,
                                    layers,
                                    weights,
                                    dummy_topic);
        if (label_post.has_nan()){
          label_post = arma::ones<arma::rowvec>(K)/K;
          //Rcout << "Encountered Nan in Label Probs" << std::endl;
        }
        gamma_k = (alpha.row(ndk_id) + cur_ndk - cur_counts) %
          ((beta.row(gene)+n_wk.row(gene)- cur_counts) /
            (beta_sum+nwk_sum- cur_counts)) %
              label_post;


        gamma_k = gamma_k/sum(gamma_k);
        //sequential par
        if (arma::any(gamma_k < 0)){
          gamma_k = (alpha.row(ndk_id) + cur_ndk) % (beta.row(gene)+
            n_wk.row(gene)) /
              (beta_sum+nwk_sum) %
                label_post;
        }
        if (gamma_k.has_nan()){
          Rcout << "cell:  " << ndk_id << " gene: " << gene << " gamma: " << gamma_k << std::endl;
          Rcout << "n_wk:" << std::endl;
          Rcout << n_wk.row(gene) << std::endl;
          Rcout << "n_dk" << std::endl;
          Rcout << n_dk.row(ndk_id) << std::endl;
          Rcout << "P(Label|Z)" << std::endl;
          Rcout << label_post << std::endl;
          stop("Variational Estimates contain NAN");
        }

        cur_gamma.row(token) = gamma_k;

      }

      n_dk.row(ndk_id).zeros();

      for (int token = 0; token < tokens; token++){ //TODO vectorise
        n_dk.row(ndk_id) += cur_gamma.row(token)*m(token,0);
      }

      d_diff = std::abs(arma::accu(cur_ndk - n_dk.row(ndk_id)));

    }
    STM_CellMap[i].cell_gamma = cur_gamma;
  }
  //Rcout << "finished iteration" << std::endl;
  build_nwk_stm(n_wk,STM_CellMap,D,K);
  build_ndk_stm(n_dk,STM_CellMap,D,K);

  elbo =  get_stm_ELBO(STM_CellMap,
                            alpha,
                            beta,
                            M,D,K,
                            n_dk,
                            n_wk,
                            C,
                            elbo_mat,
                            iter);
  STM_CellMap.clear();

  return elbo;
}

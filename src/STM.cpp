#include "CellMap.h"
#include "STM.h"
#include <cstdlib>
#include <omp.h>
#include <fstream>

using namespace Rcpp;
using namespace std;

// Script for STM Implementation

// Supervised Cell Map
STM_Cell::STM_Cell() : Cell(){
  cell_label = 0;
}

STM_Cell::STM_Cell(int id,arma::mat cellmtx,
                   const arma::mat& alpha,
                   const arma::mat& beta,int K,
                   bool zero_gamma = false,
                   bool rand_gamma = true,
                   int label = 0) : Cell(id,cellmtx,
                   alpha,
                   beta,
                   K,
                   zero_gamma,
                   rand_gamma){
  cell_label = label;
}

STM_Cell::STM_Cell(int id,arma::mat mtx,
                   int K,
                   int label): Cell::Cell(id,mtx,
                   K) {
                     cell_label = label;
                   }


void STM_Cell::print(){
  Rcout << "cell id: " << cell_id << std::endl;
  Rcout << "Label: " << cell_label << std::endl;
  Rcout << "mtx: \n" << cell_mtx << std::endl;
  Rcout << "gamma: " << cell_gamma << std::endl;
}


// [[Rcpp::export]]
void test_stm_cell(){
  STM_Cell test_cell = STM_Cell();
  test_cell.print();
}


// build STM Cell Map
std::unordered_map<int,STM_Cell> build_STM_Cell_Map(arma::sp_mat& counts,
                                                    arma::vec& celltypes,
                                                    arma::vec& genes,
                                                    const arma::mat& alpha,
                                                    const arma::mat& beta,
                                                    int D, int K,
                                                    bool zero_gamma,
                                                    bool rand_gamma,
                                                    arma::vec& labels){
  std::unordered_map<int,STM_Cell> STM_CellMap;

  for (int i = 0; i < D; i++){
    arma::uvec geneidx = arma::find(counts.col(i));
    arma::mat mtx(geneidx.n_elem,3);
    int cur_label = labels(i);
    mtx.col(0) = arma::nonzeros(counts.col(i));//counts(nonzero_counts);
    mtx.col(1).fill(celltypes(i));
    mtx.col(2) = genes(geneidx);

    STM_CellMap[i+1] = STM_Cell(i+1,mtx,alpha,beta,K,
                                zero_gamma,
                                rand_gamma,
                                cur_label);
  }

  return STM_CellMap;
}

// Test STM Cell Map
// [[Rcpp::export]]
void test_stm_cellmap(arma::sp_mat& counts,
                      arma::vec& celltypes,
                      arma::vec& genes,
                      const arma::mat& alpha,
                      const arma::mat& beta,
                      int D, int K,
                      bool zero_gamma,
                      bool rand_gamma,
                      arma::vec& labels,
                      int print_cell = 5){
  std::unordered_map<int,STM_Cell> test_map = build_STM_Cell_Map(counts,
                                                                 celltypes,
                                                                 genes,
                                                                 alpha,
                                                                 beta,
                                                                 D,K,
                                                                 zero_gamma,
                                                                 rand_gamma,
                                                                 labels);
  test_map[print_cell].print();
  Rcout << "function ran" << std::endl;
}

// STM  Helpers (softmax)




arma::rowvec get_label_prob(arma::rowvec ndk_row,
                            const arma::rowvec& gamma_weight,
                            arma::mat& model_weights,
                            int gene_count,int K, int label){

  arma::mat label_prob = arma::zeros<arma::mat>(K,K);

  ndk_row -= gamma_weight*gene_count; //remove cur value
  arma::rowvec cur_ndk;

  /// TODO
  for (int i = 0; i < K; i++){
    //adjust sufficient statistics

    cur_ndk = ndk_row;

    cur_ndk(i) += gene_count;
    cur_ndk /= arma::sum(cur_ndk);


    label_prob.row(i) = cur_ndk;

  }
  //////

  if (label_prob.has_nan()){
    Rcout << "Label Prob Fun : \n"<< label_prob << std::endl;
    stop("label probs contain NAN");
  }

  arma::mat model_out = label_prob * model_weights; // K x L
  if (model_out.has_nan()){
    Rcout << "Model out Fun : \n"<< model_out << std::endl;
    stop("model out contains NAN");
  }
  arma::rowvec final_prob = arma::zeros<arma::rowvec>(K); // K x 1 (rowvec)
  for (int i = 0; i < K; i++){
    //adjust sufficient statistics

    arma::rowvec softmax_val = arma::exp(model_out.row(i)); //1 x L

    if (softmax_val.has_nan()){
      Rcout << "softmax val 1 : \n"<< softmax_val << std::endl;
      stop("softmax contains NAN");
    }
    if (softmax_val.has_inf()){
      Rcout << "Inf softmax val 1 : \n"<< softmax_val << std::endl;
      stop("Infinite softmax val");
    }
    final_prob(i) = softmax_val(label)/arma::sum(softmax_val); // only get val for label of interest

    if (final_prob.has_nan()){
      Rcout << "final Prob Fun Partial : \n"<< final_prob << std::endl;
      stop("Intermediate probability contains NAN");
    }
  }
  final_prob /= arma::sum(final_prob); // normalize across topic assignments //1 x K
  if (final_prob.has_nan()){
    Rcout << "final Prob Fun full : \n"<< final_prob << std::endl;
    stop("Final probability contains NAN");
  }
  // get weights
  return final_prob;
}


////////////////////
// Epoch train

void run_stm_epoch(std::unordered_map<int,STM_Cell>& STM_CellMap,
                   arma::mat& model_weights,
                   const arma::mat& alpha,
                   const arma::mat& beta,int K,int D,arma::mat& n_dk,arma::mat& n_wk,
                   int num_threads = 1){
  omp_set_num_threads(num_threads);


  arma::rowvec beta_sum = arma::sum(beta,0);
  arma::rowvec nwk_sum = arma::sum(n_wk,0);
  int gene,counts,ndk_id;
  arma::rowvec gamma_k = arma::zeros<arma::rowvec>(K);
  arma::rowvec cur_counts = arma::zeros<arma::rowvec>(K);
  arma::rowvec cur_ndk = arma::zeros<arma::rowvec>(K);
  arma::rowvec label_post = arma::zeros<arma::rowvec>(K);
#pragma omp parallel for private(label_post,cur_ndk,gamma_k,cur_counts,gene,counts,ndk_id)
  for (int i = 1; i <= D; i++){

    arma::mat cur_gamma = STM_CellMap[i].cell_gamma; //TODO
    arma::mat m = STM_CellMap[i].cell_mtx; //TODO

    ndk_id = i-1; //To accomodate row indices for matrix
    int tokens = m.n_rows;
    cur_ndk = n_dk.row(ndk_id);

    for (int token = 0; token < tokens; token++){


      counts = m(token,0);
      gene = m(token,2)-1;
      cur_counts = cur_gamma.row(token)*counts; //TODO


      label_post = get_label_prob(cur_ndk,cur_gamma.row(token),
                                  model_weights,counts,K,
                                  STM_CellMap[i].cell_label);

      gamma_k = (alpha.row(ndk_id) + n_dk.row(ndk_id) - cur_counts) %
        ((beta.row(gene)+n_wk.row(gene)- cur_counts) /
          (beta_sum+nwk_sum- cur_counts)) %
            label_post;


      gamma_k = gamma_k/sum(gamma_k);
      //sequential par
      if (arma::any(gamma_k < 0)){
        gamma_k = (alpha.row(ndk_id) + n_dk.row(ndk_id)) % (beta.row(gene)+
          n_wk.row(gene)) /
            (beta_sum+nwk_sum) %
              label_post;
      }

      if (gamma_k.has_nan()){
        Rcout << "Cell: " << i  << std::endl;
        Rcout << "token : " << token << std::endl;
        Rcout << "current gamma: " << gamma_k << std::endl;
        stop("Variational Estimates contain NAN");
      }

      cur_gamma.row(token) = gamma_k; //TODO

      //update params

      // n_wk.row(gene) += (gamma_k * counts); // NOTE: this line is moved inside the critical block
      //
      // n_dk.row(ndk_id) += (gamma_k * counts) - cur_counts ;
      // nwk_sum += (gamma_k * counts);

      if (n_wk.has_nan()){
        Rcout << "Cell: " << i  << std::endl;
        Rcout << "token : " << token << std::endl;
        Rcout << "NWK: " << n_wk.row(gene) << std::endl;
        stop("NWK Estimates contain NAN");
      }

      if (n_dk.has_nan()){
        Rcout << "Cell: " << i  << std::endl;
        Rcout << "token : " << token << std::endl;
        Rcout << "NDK: " << n_dk.row(ndk_id) << std::endl;
        stop("NDK Estimates contain NAN");
      }

    }
    STM_CellMap[i].cell_gamma = cur_gamma; //TODO


  }



}

// Suff Stat Helpers
void build_nwk_stm(arma::mat& n_wk,
                   std::unordered_map<int,STM_Cell>& STM_CellMap,
                   int D, int K){
  int tokens;
  n_wk.zeros();
  for (int j = 1; j <= D; j++){ //Note: indexing starts at one to match cell ids
    tokens = STM_CellMap[j].cell_mtx.n_rows;
    //For each token, update the corresponding gene value
    for (int i = 0;i < tokens; i++){
      int gene = STM_CellMap[j].cell_mtx(i,2)-1;
      n_wk.row(gene) += STM_CellMap[j].cell_gamma.row(i) * STM_CellMap[j].cell_mtx(i,0);
    }
  }

}

void build_ndk_stm(arma::mat& n_dk,
                   std::unordered_map<int,STM_Cell>& STM_CellMap,
                   int D, int K){
  int tokens;
  n_dk.zeros();
  for(int j = 1; j <= D; j++){ //Note: indexing starts at one to match cell ids
    tokens  = STM_CellMap[j].cell_mtx.n_rows;
    //For each token, update the corresponding cell value
    for (int i = 0; i < tokens; i++){
      int cell = STM_CellMap[j].cell_mtx(i,1)-1;
      n_dk.row(cell) += STM_CellMap[j].cell_gamma.row(i) * STM_CellMap[j].cell_mtx(i,0);
    }
  }
}



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
                int cur_iter
                ){
  double ELBO_theta = 0;
  double ELBO_phi;
  double ELBO_phi_prior = 0;
  double ELBO_phi_like = 0;
  double ELBO_phi_norm = 0;
  double ELBO_gamma = 0;
  int tokens,theta_denom;
  // ELBO theta
  ELBO_theta += (lgamma(arma::accu(alpha)) - arma::accu(lgamma(alpha)));
  for (int i = 0; i < D; i++){
    theta_denom = lgamma(arma::sum(alpha.row(i)) + arma::sum(n_dk.row(i)));
    for (int j = 0; j < K; j++){
      ELBO_theta += (lgamma(alpha(i,j) + n_dk(i,j))) - theta_denom;
    }
  }

  // ELBO phi
  ELBO_phi_prior += (lgamma(arma::accu(beta)) - arma::accu(lgamma(beta)));

  for (int j = 0; j < M; j++){
    ELBO_phi_like += arma::accu(lgamma(beta.row(j) + n_wk.row(j)));
  }
  for (int j = 0; j < M; j++){
    //for (int i = 0; i < K; i++){
      ELBO_phi_norm += lgamma(arma::sum(beta.row(j)) + arma::sum(n_wk.row(j)));
    //}
  }

  ELBO_phi = (K*ELBO_phi_prior + ELBO_phi_like - ELBO_phi_norm);

  // ELBO gamma
  for (int cell = 1; cell <= D; cell++){
    tokens = STM_CellMap[cell].cell_mtx.n_rows;
    for (int token = 0; token < tokens; token++){
      ELBO_gamma += arma::sum(STM_CellMap[cell].cell_gamma.row(token) %
        arma::log(STM_CellMap[cell].cell_gamma.row(token) + 0.0000001)*STM_CellMap[cell].cell_mtx(token,0));
    }
  }
  //Rcout << ELBO_theta/C << ELBO_phi/C << ELBO_gamma/C << std::endl;
  arma::rowvec cur_val = {ELBO_theta/C,ELBO_phi/C,ELBO_gamma/C, (ELBO_theta + ELBO_phi -ELBO_gamma)/C};
  elbo.row(cur_iter).subvec(0,3) = cur_val;
  return (ELBO_theta + ELBO_phi + -ELBO_gamma)/C;
}



// Initialize + Iterate through weights
//////////////////////////////////////////

// [[Rcpp::export]]
arma::mat train_stm(arma::sp_mat& counts,
                          arma::vec& celltypes,
                          arma::vec& genes,
                          arma::vec& labels,
                          const arma::mat& alpha, const arma::mat& beta,
                          int K,int D,arma::mat& n_dk,arma::mat& n_wk,
                          int num_threads = 1,int maxiter = 100,
                          bool verbal = true,
                          bool zero_gamma = false,
                          bool rand_gamma = true,
                          double thresh = 0.00001,
                          double lr = 0.0001){
  double old_elbo = 0;
  double cur_elbo = 0;
  arma::mat elbo_val = arma::zeros<arma::mat>(maxiter,5);
  int M = genes.n_elem;
  int C = arma::accu(counts);
  double prog = 0;
  int prog_width = 50;
  int pos = 0;
  std::unordered_map<int,STM_Cell> STM_CellMap = build_STM_Cell_Map(counts,
                                                                    celltypes,
                                                                    genes,
                                                                    alpha,
                                                                    beta,
                                                                    D,K,
                                                                    zero_gamma,
                                                                    rand_gamma,
                                                                    labels);
  double ce_loss = 0;
  if (verbal){

  Rcout << "STM_CellMap is built" << std::endl;
  }

  //Complete // Consider making the unordered_map its own object
  build_nwk_stm(n_wk,STM_CellMap,D,K);
  //double lr = 0.0001;
  //TComplete // Consider making the unordered_map its own object
  build_ndk_stm(n_dk,STM_CellMap,D,K);

  // Initialize Weights
  arma::vec unique_label = arma::unique(labels);
  int L = unique_label.n_elem;
  arma::mat model_weights = arma::randn(K,L);
  arma::mat model_grad = arma::zeros(K,L);
  // store prob for labels
  arma::mat label_prob = arma::zeros(D,L);
  arma::mat ce_temp = arma::zeros(D,L);
  arma::mat label_true = arma::zeros(D,L);

  arma::mat ndk_norm = arma::zeros(D,K);

  for (int i = 0; i < D; i++){
    label_true.row(i)(labels(i)) = 1;
  }


  for (int i = 0; i < maxiter; i++){
    run_stm_epoch(STM_CellMap,
                  model_weights,
                  alpha,beta,K,D,n_dk,n_wk,
                  num_threads);
    build_nwk_stm(n_wk,STM_CellMap,D,K);
    build_ndk_stm(n_dk,STM_CellMap,D,K);

    ndk_norm = n_dk;

    for (int d = 0; d < D; d++){
      ndk_norm.row(d) /= arma::sum(ndk_norm.row(d));
    }

    int m_step_iter = 0;
    double cur_ce = 99999;
    double ce_diff = 999;
    while ((m_step_iter < 100) & (ce_diff > 0.00001)){
      label_prob = ndk_norm * model_weights;
      label_prob = arma::exp(label_prob);
      for (int d = 0; d < D; d++){
        label_prob.row(d) /= arma::sum(label_prob.row(d));
      }
      model_grad = trans(ndk_norm)*(label_prob-label_true);
      model_weights -= lr*model_grad;
      ce_temp = label_true % arma::log(label_prob);
      ce_loss = -(1.0/D)*arma::accu(ce_temp);
      ce_diff = abs(cur_ce-ce_loss);
      cur_ce = ce_loss;
      m_step_iter++;
    }


    cur_elbo = get_stm_ELBO(STM_CellMap,
                        alpha,
                        beta,
                        M,D,K,
                        n_dk,
                        n_wk,
                        C,
                        elbo_val,
                        i);//(i+1);
    cur_elbo -= ce_loss;
    elbo_val.row(i)[4] = -ce_loss;
    if (i == 0){old_elbo = cur_elbo;}

    old_elbo = cur_elbo;
    if (verbal){
      //Calc cross entropy
      Rcout.flush();

      prog += 1.0/maxiter;
      Rcout << "[";
      pos = int(prog_width * prog);
      for(int p = 0; p < prog_width; p++){
        if (p < pos){ Rcout << "=" ;}
        else { Rcout << " " ;}
      }
      Rcout << "] " << int(prog * 100.0) << "% || Iter: " << i << " || CE: "<< ce_loss << " || ELBO: "<< cur_elbo <<" \r";
    }

  }


  if (verbal){
    Rcout.flush();
    Rcout << "[";
    for(int p = 0; p < prog_width; p++){
      Rcout << "=" ;
    }
    Rcout << "] " << "100% || Iter: " << maxiter << " || CE: "<< ce_loss << " || ELBO: "<< cur_elbo << std::endl;
    Rcout << "Max Iteration Reached" << std::endl;
  }
  return model_weights;
}

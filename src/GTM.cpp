#include "CellMap.h"
#include "GTM.h"
#include <omp.h>

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;



//==========================
/*
 infer sufficient statistics
 */


void build_nwk(arma::mat& n_wk,
               std::unordered_map<int,Cell>& CellMap,
               int D, int K){
  int tokens;
  n_wk.zeros();
  for (int j = 1; j <= D; j++){ //Note: indexing starts at one to match cell ids
    tokens = CellMap[j].cell_mtx.n_rows;
    //For each token, update the corresponding gene value
    for (int i = 0;i < tokens; i++){
      int gene = CellMap[j].cell_mtx(i,2)-1;
      n_wk.row(gene) += CellMap[j].cell_gamma.row(i) * CellMap[j].cell_mtx(i,0);
    }
  }

}

void build_ndk(arma::mat& n_dk,
               std::unordered_map<int,Cell>& CellMap,
               int D, int K){
  int tokens;
  n_dk.zeros();
  for(int j = 1; j <= D; j++){ //Note: indexing starts at one to match cell ids
    tokens  = CellMap[j].cell_mtx.n_rows;
    //For each token, update the corresponding cell value
    for (int i = 0; i < tokens; i++){
      int cell = CellMap[j].cell_mtx(i,1)-1;
      n_dk.row(cell) += CellMap[j].cell_gamma.row(i) * CellMap[j].cell_mtx(i,0);
    }
  }
}


// [[Rcpp::export]]
arma::mat get_theta(const arma::mat& n_dk,
                    const arma::mat& alpha){
  //Rcout << "1" << std::endl;
  arma::mat theta = arma::zeros<arma::mat>(alpha.n_rows, alpha.n_cols);
  //Rcout << "2" << std::endl;
  int cells = n_dk.n_rows;
  for (int i = 0; i < cells; i++){
    theta.row(i) = (alpha.row(i) + n_dk.row(i))/(arma::sum(alpha.row(i)) +
      arma::sum(n_dk.row(i)));
  }
  return theta;
}

// [[Rcpp::export]]
arma::mat get_phi(const arma::mat& n_wk,
                  const arma::mat& beta){
  arma::mat phi = arma::zeros<arma::mat>(beta.n_rows,beta.n_cols);
  arma::rowvec phi_denom = arma::sum(beta,0) + arma::sum(n_wk,0);
  int genes = n_wk.n_rows;
  for (int i = 0; i < genes; i++){
    phi.row(i) = (beta.row(i) + n_wk.row(i))/phi_denom;
  }
  return phi;
}
//==========================


double get_ELBO(std::unordered_map<int,Cell>& CellMap,
                const arma::mat& alpha,
                const arma::mat& beta,
                int M,
                int D,
                int K,
                arma::mat& n_dk,
                arma::mat& n_wk,
                int C,
                arma::mat& elbo,
                int cur_iter){
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
  //theta_denom = lgamma(arma::accu(alpha) + arma::accu(n_dk));
  //ELBO_theta -= theta_denom;

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
    tokens = CellMap[cell].cell_mtx.n_rows;
    for (int token = 0; token < tokens; token++){
      ELBO_gamma += arma::sum(CellMap[cell].cell_gamma.row(token) %
        arma::log(CellMap[cell].cell_gamma.row(token) + 0.0000001)*CellMap[cell].cell_mtx(token,0));
    }
  }
  //Rcout << ELBO_theta/C << ELBO_phi/C << ELBO_gamma/C << std::endl;
  arma::rowvec cur_val = {ELBO_theta/C,ELBO_phi/C,ELBO_gamma/C, (ELBO_theta + ELBO_phi -ELBO_gamma)/C};
  elbo.row(cur_iter).subvec(0,3) = cur_val;
  return (ELBO_theta + ELBO_phi + -ELBO_gamma)/C;
}


void run_epoch(std::unordered_map<int,Cell>& CellMap,const arma::mat& alpha,
               const arma::mat& beta,int K,int D,arma::mat& n_dk,arma::mat& n_wk,
               int num_threads = 1){
  //omp_set_num_threads(num_threads);


  arma::rowvec beta_sum = arma::sum(beta,0);
  arma::rowvec nwk_sum = arma::sum(n_wk,0);
  //int gene,counts,ndk_id;
  arma::rowvec gamma_k = arma::zeros<arma::rowvec>(K);
  arma::rowvec cur_counts = arma::zeros<arma::rowvec>(K);
  //#pragma omp parallel for private(cur_counts,gamma_k) reduction(+:n_wk) reduction(+:nwk_sum)

  for (int i = 1; i <= D; i++){

    arma::mat cur_gamma = CellMap[i].cell_gamma; //TODO
    arma::mat m = CellMap[i].cell_mtx; //TODO
    int ndk_id = i-1; //To accomodate row indices for matrix
    int tokens = m.n_rows;
    double diff = 1;
    int max_update = 1;
    //int u = 0;
    //arma::rowvec cur_ndk = n_dk.row(ndk_id);
    //Gamma burn-in
    //while (diff > 0.01 && u < max_update){
    //  u++;
    for (int token = 0; token < tokens; token++){

      //Rcout << "1" << std::endl;
      int counts = m(token,0);
      int gene = m(token,2)-1;
      arma::rowvec cur_counts = cur_gamma.row(token)*counts; //TODO
      // Rcout << "2" << std::endl;

      // E-step
      arma::rowvec gamma_k = (alpha.row(ndk_id) + n_dk.row(ndk_id) - cur_counts) %
        (beta.row(gene)+n_wk.row(gene)- cur_counts) /
          (beta_sum+nwk_sum- cur_counts);
      gamma_k = gamma_k/sum(gamma_k);
      //sequential par
      if (arma::any(gamma_k < 0)){
        gamma_k = (alpha.row(ndk_id) + n_dk.row(ndk_id)) % (beta.row(gene)+
          n_wk.row(gene)) /
            (beta_sum+nwk_sum);
      }
      //Rcout << "3" << std::endl;
      cur_gamma.row(token) = gamma_k; //TODO

      //M-Step
      n_dk.row(ndk_id) += (gamma_k * counts) - cur_counts;
      n_wk.row(gene) += (gamma_k * counts);// - cur_counts; // NOTE: this line is moved inside the critical block
      nwk_sum += (gamma_k * counts); //- cur_counts;

    }
    //  diff = abs(arma::accu(cur_ndk - n_dk.row(ndk_id)));
    //  cur_ndk = n_dk.row(ndk_id);
    //}
    CellMap[i].cell_gamma = cur_gamma; //TODO
  }



}


/* Full model

 */


// [[Rcpp::export]]
void train_gtm(arma::sp_mat& counts,
                     arma::vec& celltypes,
                     arma::vec& genes,
                     const arma::mat& alpha, const arma::mat& beta,
                     int K,int D,arma::mat& n_dk,arma::mat& n_wk,
                     int num_threads = 1,int maxiter = 100,
                     bool verbal = true,
                     bool zero_gamma = false,
                     bool rand_gamma = true,
                     double thresh = 0){
  double old_elbo;
  double cur_elbo = 0;
  int M = genes.n_elem;
  int C = arma::accu(counts);
  double prog = 0;
  int prog_width = 50;
  int pos = 0;
  std::unordered_map<int,Cell> CellMap = build_Cell_Map(counts,
                                                        celltypes,
                                                        genes,
                                                        alpha,
                                                        beta,
                                                        D,K,
                                                        zero_gamma,
                                                        rand_gamma);

  Rcout << "CellMap is built" << std::endl;
  //Complete
  build_nwk(n_wk,CellMap,D,K);

  //TComplete
  build_ndk(n_dk,CellMap,D,K);


  arma::mat elbo_val = arma::zeros<arma::mat>(maxiter,4);
  for (int i = 0; i < maxiter; i++){
    run_epoch(CellMap,alpha,beta,K,D,n_dk,n_wk,
              num_threads);
    //build_nwk(n_wk,CellMap,D,K);



    cur_elbo = get_ELBO(CellMap,
                        alpha,
                        beta,
                        M,D,K,
                        n_dk,
                        n_wk,
                        C,
                        elbo_val,
                        i);//(i+1);
    if (verbal){

      Rcout.flush();

      prog += 1.0/maxiter;
      Rcout << "[";
      pos = int(prog_width * prog);
      for(int p = 0; p < prog_width; p++){
        if (p < pos){ Rcout << "=" ;}
        else { Rcout << " " ;}
      }
      Rcout << "] " << int(prog * 100.0) << "% || Iter: " << i << " || ELBO: "<< cur_elbo <<" \r";
    }

    if (i > 0){
      if (abs(cur_elbo - old_elbo) < thresh){
        Rcout.flush();
        Rcout << "[";
        for(int p = 0; p < prog_width; p++){
          Rcout << "=" ;
        }
        Rcout << "] " << "100% || Iter: " << i << " || ELBO: "<< cur_elbo << std::endl;
        Rcout << "Model Converged at iteration : " << i << std::endl;

        return;
      }
    }

    old_elbo = cur_elbo;

  }
  Rcout.flush();
  Rcout << "[";
  for(int p = 0; p < prog_width; p++){
    Rcout << "=" ;
  }
  Rcout << "] " << "100% || Iter: " << maxiter << " || ELBO: "<< cur_elbo << std::endl;
  Rcout << "Max Iteration Reached" << std::endl;

  return;

}


/// PREDICTION
void predict_epoch(std::unordered_map<int,Cell>& CellMap,const double& alpha,
                   int K,int D,int M,arma::mat& n_dk,
                   const arma::mat& phi, int num_threads = 1){
  //omp_set_num_threads(num_threads);
  int gene,counts,ndk_id;
  arma::rowvec gamma_k = arma::zeros<arma::rowvec>(K);
  arma::rowvec cur_counts = arma::zeros<arma::rowvec>(K);
  //#pragma omp parallel for
  for (int i = 1; i <= D; i++){

    arma::mat cur_gamma = CellMap[i].cell_gamma; //TODO
    arma::mat m = CellMap[i].cell_mtx; //TODO
    ndk_id = i-1; //To accomodate row indices for matrix
    int tokens = m.n_rows;
    for (int token = 0; token < tokens; token++){
      //Rcout << "1" << std::endl;
      counts = m(token,0);
      gene = m(token,2)-1;
      cur_counts = cur_gamma.row(token)*counts; //TODO
      // Rcout << "2" << std::endl;
      gamma_k = (alpha + n_dk.row(ndk_id) - cur_counts) % phi.row(gene);

      //sequential par
      if (arma::any(gamma_k < 0)){
        gamma_k = (alpha + n_dk.row(ndk_id)) % phi.row(gene);
      }
      gamma_k = gamma_k/sum(gamma_k);
      //Rcout << "3" << std::endl;
      cur_gamma.row(token) = gamma_k; //TODO

      //update params
      //Rcout << "4" << std::endl;

      n_dk.row(ndk_id) += (gamma_k * counts) - cur_counts ;
    }
    CellMap[i].cell_gamma = cur_gamma; //TODO


  }
}


// [[Rcpp::export]]
void infer_topics_cpp(arma::sp_mat& counts,
                       arma::vec& celltypes,
                       arma::vec& genes,
                       double& alpha,
                       int K,int D,int M,arma::mat& n_dk, const arma::mat& phi,
                       int num_threads = 1,int maxiter = 100,
                       bool verbal = true){
  std::unordered_map<int,Cell> CellMap = build_Predict_Cell_Map(counts,
                                                                celltypes,
                                                                genes,
                                                                D,K);
  //Rcout << "Prediction CellMap is built" << std::endl;
  // CellMap[1].print();
  // CellMap[29].print();
  build_ndk(n_dk,CellMap,D,K);
  //Rcout << "Initialized NDK" << std::endl;

  double prog = 0;
  int prog_width = 50;
  int pos = 0;

  for (int i = 0; i < maxiter; i++){
    predict_epoch(CellMap,alpha,K,D,M,n_dk,phi,num_threads);
    if (verbal){

      Rcout.flush();

      prog += 1.0/maxiter;
      Rcout << "[";
      pos = int(prog_width * prog);
      for(int p = 0; p < prog_width; p++){
        if (p < pos){ Rcout << "=" ;}
        else { Rcout << " " ;}
      }
      Rcout << "] " << int(prog * 100.0) << "% || Iter: " << i << " \r";
    }
  }
  if (verbal){
    Rcout.flush();
    Rcout << "[";
    for(int p = 0; p < prog_width; p++){
      Rcout << "=" ;
    }
    Rcout << "] " << "100% || Iter: " << maxiter << std::endl;
    Rcout << "Max Iteration Reached" << std::endl;
  }
  return;
}

#include "CellMap.h"
#include "RTM.h"
#include <cstdlib>
#include <omp.h>
#include <fstream>

using namespace Rcpp;
using namespace std;
// Script for RTM Implementation

// Supervised Cell Map
RTM_Cell::RTM_Cell() : Cell(){
}

RTM_Cell::RTM_Cell(int id,arma::mat cellmtx,
                   const arma::mat& alpha,
                   const arma::mat& beta,int K,
                   const arma::uvec& nbrs,
                   bool zero_gamma = false,
                   bool rand_gamma = true) : Cell(id,
                   cellmtx,
                   alpha,
                   beta,
                   K,
                   zero_gamma,
                   rand_gamma){
  //Rcout << "In function" << std::endl;
  cell_nbrs = nbrs;
  //Rcout << "Correct behaviour" << std::endl;
}
RTM_Cell::RTM_Cell(int id,arma::mat cellmtx,
                   const arma::mat& alpha,
                   const arma::mat& beta,int K,
                   const arma::uvec& nbrs,
                   const arma::vec& nbr_bool,
                   bool zero_gamma = false,
                   bool rand_gamma = true) : Cell(id,
                   cellmtx,
                   alpha,
                   beta,
                   K,
                   zero_gamma,
                   rand_gamma){
  //Rcout << "In function" << std::endl;
  cell_nbrs = nbrs;
  bool_nbrs = nbr_bool;
  //Rcout << "Correct behaviour" << std::endl;
}
RTM_Cell::RTM_Cell(int id,arma::mat cellmtx,
                   const arma::mat& alpha,
                   const arma::mat& beta,int K,
                   const arma::uvec& nbrs,
                   const arma::vec& nbr_bool,
                   const arma::vec& cell_dist,
                   bool zero_gamma = false,
                   bool rand_gamma = true) : Cell(id,
                   cellmtx,
                   alpha,
                   beta,
                   K,
                   zero_gamma,
                   rand_gamma){
  //Rcout << "In function" << std::endl;
  cell_nbrs = nbrs;
  bool_nbrs = nbr_bool;
  nbr_dist = cell_dist;
  //Rcout << "Correct behaviour" << std::endl;
}


RTM_Cell::RTM_Cell(int id,arma::mat mtx,
                   int K,
                   const arma::uvec& nbrs): Cell::Cell(id,
                   mtx,
                   K) {
                     cell_nbrs = nbrs;
                   }

RTM_Cell::RTM_Cell(int id,arma::mat mtx,
                   int K,
                   const arma::uvec& nbrs,
                   const arma::vec& nbr_bool): Cell::Cell(id,
                   mtx,
                   K) {
                     cell_nbrs = nbrs;
                     bool_nbrs = nbr_bool;
                   }



void RTM_Cell::print(){
  Rcout << "cell id: " << cell_id << std::endl;
  Rcout << "Neighbors: " << cell_nbrs << std::endl;
  Rcout << "Neighbor Distance: " << nbr_dist << std::endl;
  Rcout << "mtx: \n" << cell_mtx << std::endl;
  Rcout << "gamma: " << cell_gamma << std::endl;
  Rcout << "nbr_bool" << bool_nbrs << std::endl;
}


// [[Rcpp::export]]
void test_RTM_cell(){
  RTM_Cell test_cell = RTM_Cell();
  test_cell.print();
}


// build RTM Cell Map
std::unordered_map<int,RTM_Cell> build_RTM_Cell_Map(arma::sp_mat& counts,
                                                    arma::vec& celltypes,
                                                    arma::vec& genes,
                                                    const arma::mat& alpha,
                                                    const arma::mat& beta,
                                                    int D, int K,
                                                    bool zero_gamma,
                                                    bool rand_gamma,
                                                    Rcpp::List& nbr_list,
                                                    int loss_fun){
  //Rcout << "here 3" << std::endl;
  std::unordered_map<int,RTM_Cell> RTM_CellMap;
  //Rcout << "here 4" << std::endl;
  for (int i = 0; i < D; i++){
    //Rcout << "here 5" << std::endl;
    arma::uvec geneidx = arma::find(counts.col(i));
    arma::mat mtx(geneidx.n_elem,3);
    //Rcout << "here 6" << std::endl;
    mtx.col(0) = arma::nonzeros(counts.col(i));//counts(nonzero_counts);
    mtx.col(1).fill(celltypes(i));
    mtx.col(2) = genes(geneidx);
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
    //Rcout << "here 7" << std::endl;
    if (loss_fun == 0){
      arma::mat nbr_mat = Rcpp::as<arma::mat>(nbr_list[i]);
      arma::vec temp_nbr = nbr_mat.col(0);
      arma::uvec nbr_ids = arma::conv_to<arma::uvec>::from(temp_nbr);
      RTM_CellMap[i+1] = RTM_Cell(i+1,mtx,alpha,beta,K,
                                  nbr_ids,
                                  zero_gamma,
                                  rand_gamma);
    }
    else {
      //nbr_list is a list of matrices, need first column
      arma::mat nbr_mat = Rcpp::as<arma::mat>(nbr_list[i]);
      arma::vec temp_nbr = nbr_mat.col(0);
      arma::uvec nbr_ids = arma::conv_to<arma::uvec>::from(temp_nbr);
      arma::vec nbr_bool = arma::conv_to<arma::vec>::from(nbr_mat.col(1));
      if (loss_fun == 1){
        RTM_CellMap[i+1] = RTM_Cell(i+1,mtx,alpha,beta,K,
                                    nbr_ids,
                                    nbr_bool,
                                    zero_gamma,
                                    rand_gamma);
      }
      else {
        //Rcout << "here 1" << std::endl;
        arma::vec nbr_dist = arma::conv_to<arma::vec>::from(nbr_mat.col(2));
        //Rcout << "here 2" << std::endl;
        RTM_CellMap[i+1] = RTM_Cell(i+1,mtx,alpha,beta,K,
                                    nbr_ids,
                                    nbr_bool,
                                    nbr_dist,
                                    zero_gamma,
                                    rand_gamma);
      }
    }


   // Rcout << "After function: " << i+1 << std::endl;
  }
  //Rcout << "Complete" << std::endl;
  return RTM_CellMap;
}

// Test RTM Cell Map
// [[Rcpp::export]]
void test_RTM_cellmap(arma::sp_mat& counts,
                      arma::vec& celltypes,
                      arma::vec& genes,
                      const arma::mat& alpha,
                      const arma::mat& beta,
                      int D, int K,
                      bool zero_gamma,
                      bool rand_gamma,
                      Rcpp::List& nbr_list,
                      int print_cell = 5,
                      int loss_fun = 0){
  //Rcout << "here 1" << std::endl;
  std::unordered_map<int,RTM_Cell> test_map = build_RTM_Cell_Map(counts,
                                                                 celltypes,
                                                                 genes,
                                                                 alpha,
                                                                 beta,
                                                                 D,K,
                                                                 zero_gamma,
                                                                 rand_gamma,
                                                                 nbr_list,
                                                                 loss_fun);
  //Rcout << "here 2" << std::endl;
  test_map[print_cell].print();
  Rcout << "function ran" << std::endl;
}

// RTM  Helpers (softmax)
arma::rowvec get_nbr_prob(arma::rowvec ndk_row,
                                  const arma::rowvec& gamma_weight,
                                  arma::vec& model_weights,
                                  int model_bias,
                                  arma::mat& nbr_ndk,
                                  int gene_count,int K,
                                  int loss_fun,
                                  arma::vec& nbr_bool,
                                  const arma::vec* nbr_dist = nullptr){
  arma::mat label_prob = arma::zeros<arma::mat>(K,K);
  ndk_row -= gamma_weight*gene_count; //remove cur value
  arma::rowvec cur_ndk = arma::zeros<arma::rowvec>(K);
  arma::rowvec final_prob = arma::zeros<arma::rowvec>(K); // K x 1 (rowvec)
  arma::rowvec nbr_prob = arma::zeros<arma::rowvec>(K);


  double sum_ndk = arma::sum(ndk_row);
  for (int i = 0; i < K; i++){
    //adjust sufficient statistics
    cur_ndk = ndk_row;
    cur_ndk(i) += gene_count;
    cur_ndk /= arma::sum(cur_ndk);
    label_prob.row(i) = cur_ndk;
  }
  //Rcout << label_prob << std::endl;
  //Rcout << nbr_ndk << std::endl;
  //Rcout << "label loop start" << std::endl;
  for (arma::uword i = 0; i < nbr_ndk.n_rows; i++){
    //Rcout << "label-topic loop start" << std::endl;
    for (int j = 0; j < K; j++){
      if (loss_fun < 2){
        nbr_prob(j) = arma::dot((label_prob.row(j) % nbr_ndk.row(i)),model_weights);
        nbr_prob(j) += model_bias;
      } else{
        //Rcout << "B1" << std::endl;
        nbr_prob(j) = arma::dot((label_prob.row(j) % nbr_ndk.row(i)),model_weights);
        nbr_prob(j) += model_bias;

        nbr_prob(j) = ((model_weights(j) * (*nbr_dist)(i))/ sum_ndk)
          - (1.0 / (2.0 * sum_ndk * sum_ndk)) *
            (2.0 * model_weights(j) * nbr_prob(j) + nbr_prob(j) * nbr_prob(j));
        //Rcout << "B3" << std::endl;
     }
    }
    //Rcout << "label-topic loop end" << std::endl;
    //final_prob = final_prob % nbr_prob;
    if (loss_fun == 1){
      //Want to minimize negative examples and maximize positives
      final_prob += nbr_bool.at(i)*nbr_prob -(1-nbr_bool.at(i))*nbr_prob;
      //final_prob += nbr_bool.at(i)*nbr_prob +(1-nbr_bool.at(i))*(1-nbr_prob);
    }
    else {

      //Rcout << nbr_bool << std::endl;
      //Rcout << i << std::endl;
      final_prob += nbr_prob;
    }

  }
  //Rcout << "label loop end" << std::endl;
  final_prob = arma::exp(final_prob);
  final_prob /= arma::sum(final_prob);

  return final_prob;
}


//TODO////////////////////////////////
void update_weights(std::unordered_map<int,RTM_Cell>& RTM_CellMap,
                    int K,
                    int D, int rho, int num_links,
                    arma::mat& n_dk,
                    double& model_bias,
                    arma::vec& model_weights,
                    int loss_fun = 0,
                    double lr = 0.0001){
  arma::mat ndk_norm = n_dk;
  ndk_norm.each_col() /= arma::sum(ndk_norm,1);
 if (loss_fun == 0){
   arma::rowvec nbr_sums = arma::zeros<arma::rowvec>(K);
   for (int j = 0; j < D; j++){
     arma::uvec nbr_ids = RTM_CellMap[j+1].cell_nbrs;


     arma::mat sub_mat = ndk_norm.rows(nbr_ids);

     nbr_sums += arma::sum(sub_mat.each_row() %
       ndk_norm.row(j) ,0);

   }

   model_bias = (std::log(num_links -arma::accu(nbr_sums)) -
     std::log(rho*((K-1.0)/K) + num_links -arma::accu(nbr_sums)));

   model_weights = arma::log(nbr_sums.t()) -
     arma::log(nbr_sums.t() + rho/(K*K)) -
     model_bias;
 } else if (loss_fun == 1) {
   //Logistic Regression update for bias and weights
   //double lr = 0.0001;
   double old_bias = model_bias;
   arma::vec old_weights = model_weights;
   //Gradient descent update
   arma::vec  weight_grad = arma::zeros<arma::vec>(K);
   double bias_grad = 0.0;

   for (int z = 0; z < 100; z++){
     for (int j = 0; j < D; j++){
       arma::uvec nbr_ids = RTM_CellMap[j+1].cell_nbrs;
       arma::vec nbr_bool = RTM_CellMap[j+1].bool_nbrs;
       if (nbr_ids.n_elem == 0){
         continue;
       }
       arma::mat sub_mat = ndk_norm.rows(nbr_ids);
       //Rcout << "here 1" << std::endl;
       sub_mat = sub_mat.each_row() %
         ndk_norm.row(j);
       //Rcout << "here 2" << std::endl;
       //Rcout << size(sub_mat) << std::endl;
       //Rcout << size(model_weights) << std::endl;
       arma::rowvec model_output = (sub_mat*model_weights).t() + model_bias;
       //Rcout << "here 3" << std::endl;
       model_output = 1/(1+arma::exp(-model_output));
       //Rcout << "here 4" << std::endl;
       weight_grad += ((nbr_bool.t() - model_output)*sub_mat).t();
       //Rcout << "here 5" << std::endl;
       bias_grad += arma::sum(nbr_bool.t() - model_output);
       //Rcout << "here 6" << std::endl;

     }
     model_bias -= lr*bias_grad/num_links;
     model_weights -= lr*weight_grad/num_links;
   }
 } else if (loss_fun == 2){
   //TODO////////////////////////////////
   // Loop for gradient computation and weight updates
   arma::vec  weight_grad = arma::zeros<arma::vec>(K);
   double bias_grad = 0.0;
   for (int z = 0; z < 100; z++) {
     // Reset gradients at the start of each iteration
     weight_grad.zeros();
     bias_grad = 0.0;

     for (int j = 0; j < D; j++) {
       arma::uvec nbr_ids = RTM_CellMap[j+1].cell_nbrs;
       arma::vec nbr_dist = RTM_CellMap[j+1].nbr_dist;
       if (nbr_ids.n_elem == 0) {
         continue; // Skip if no neighbors
       }

       // Compute the sub-matrix and model output
       arma::mat sub_mat = ndk_norm.rows(nbr_ids);
       sub_mat = sub_mat.each_row() % ndk_norm.row(j);
       arma::rowvec model_output = (sub_mat * model_weights).t() + model_bias;

       // Calculate prediction error
       arma::vec error = model_output.t() - nbr_dist;

       // Compute gradients
       weight_grad += sub_mat.t() * error; // Gradient w.r.t. weights
       bias_grad += arma::sum(error);     // Gradient w.r.t. bias
     }

     // Apply gradient updates
     model_weights -= lr * weight_grad / num_links; // Normalize by D
     model_bias -= lr * bias_grad / num_links;      // Normalize by D
   }
 }
}
////////////////////
// Epoch train
//TODO FIX HERE
void run_RTM_epoch(std::unordered_map<int,RTM_Cell>& RTM_CellMap,
                   arma::vec& model_weights,
                   int model_bias,
                   const arma::mat& alpha,
                   const arma::mat& beta,int K,int D,arma::mat& n_dk,arma::mat& n_wk,
                   int num_threads = 1,
                   int loss_fun = 0,
                   int burnin = 1){
  omp_set_num_threads(num_threads);

  //Rcout << "here 4" << std::endl;
  arma::rowvec beta_sum = arma::sum(beta,0);
  arma::rowvec nwk_sum = arma::sum(n_wk,0);
  int gene,counts,ndk_id;
  arma::rowvec gamma_k = arma::zeros<arma::rowvec>(K);
  arma::rowvec cur_counts = arma::zeros<arma::rowvec>(K);
  arma::rowvec cur_ndk = arma::zeros<arma::rowvec>(K);
  arma::rowvec label_post = arma::zeros<arma::rowvec>(K);
#pragma omp parallel for private(label_post,cur_ndk,gamma_k,cur_counts,gene,counts,ndk_id)
  for (int i = 1; i <= D; i++){
    //Rcout << "here 5" << std::endl;
    arma::mat cur_gamma = RTM_CellMap[i].cell_gamma;
    arma::mat m = RTM_CellMap[i].cell_mtx;
    arma::uvec nbr_ids = RTM_CellMap[i].cell_nbrs;
    arma::mat nbr_ndk = n_dk.rows(nbr_ids);
    arma::vec nbr_bool = RTM_CellMap[i].bool_nbrs;

    //normalize each row
    nbr_ndk.each_col() /= arma::sum(nbr_ndk,1);

    ndk_id = i-1; //To accommodate row indices for matrix
    int tokens = m.n_rows;
    cur_ndk = n_dk.row(ndk_id);

    double diff = 1;
    int max_update = burnin;
    int u = 0;
    while (diff > 0.01 && u < max_update){
      u++;
      for (int token = 0; token < tokens; token++){


        counts = m(token,0);
        gene = m(token,2)-1;
        cur_counts = cur_gamma.row(token)*counts; //TODO

        if (nbr_ids.n_elem == 0){
          label_post = arma::ones<arma::rowvec>(K);
          label_post /= arma::sum(label_post);
        }
        else {

          if (loss_fun == 2){
            label_post = get_nbr_prob(cur_ndk,cur_gamma.row(token),
                                      model_weights,model_bias,
                                      nbr_ndk,counts,K,loss_fun,nbr_bool,&RTM_CellMap[i].nbr_dist);
          } else {
            label_post = get_nbr_prob(cur_ndk,cur_gamma.row(token),
                                      model_weights,model_bias,
                                      nbr_ndk,counts,K,loss_fun,nbr_bool);
          }
          //Rcout << "here C" << std::endl;
          //Rcout << "End label prob" << std::endl;
        }

        //Rcout << "here 6" << std::endl;
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

        cur_gamma.row(token) = gamma_k;


        if (n_wk.has_nan()){
          Rcout << "Cell: " << i  << std::endl;
          Rcout << "token : " << token << std::endl;
          Rcout << "NWK: " << n_wk.row(gene) << std::endl;
          stop("NWK matrix contains NAN");
        }

        if (n_dk.has_nan()){
          Rcout << "Cell: " << i  << std::endl;
          Rcout << "token : " << token << std::endl;
          Rcout << "NDK: " << n_dk.row(ndk_id) << std::endl;
          stop("NDK matrix contains NAN");
        }

      }
      n_dk.row(ndk_id).zeros();
      for (int token = 0; token < tokens; token++){ //TODO vectorise
        n_dk.row(ndk_id) += cur_gamma.row(token)*m(token,0);
      }

      diff = abs(arma::accu(cur_ndk - n_dk.row(ndk_id)));
      cur_ndk = n_dk.row(ndk_id);
    }
    RTM_CellMap[i].cell_gamma = cur_gamma;


  }



}

// Suff Stat Helpers
void build_nwk_RTM(arma::mat& n_wk,
                   std::unordered_map<int,RTM_Cell>& RTM_CellMap,
                   int D, int K){
  int tokens;
  n_wk.zeros();
  for (int j = 1; j <= D; j++){ //Note: indexing starts at one to match cell ids
    tokens = RTM_CellMap[j].cell_mtx.n_rows;
    //For each token, update the corresponding gene value
    for (int i = 0;i < tokens; i++){
      int gene = RTM_CellMap[j].cell_mtx(i,2)-1;
      n_wk.row(gene) += RTM_CellMap[j].cell_gamma.row(i) * RTM_CellMap[j].cell_mtx(i,0);
    }
  }

}

void build_ndk_RTM(arma::mat& n_dk,
                   std::unordered_map<int,RTM_Cell>& RTM_CellMap,
                   int D, int K){
  int tokens;
  n_dk.zeros();
  for(int j = 1; j <= D; j++){ //Note: indexing starts at one to match cell ids
    tokens  = RTM_CellMap[j].cell_mtx.n_rows;
    //For each token, update the corresponding cell value
    for (int i = 0; i < tokens; i++){
      int cell = RTM_CellMap[j].cell_mtx(i,1)-1;
      n_dk.row(cell) += RTM_CellMap[j].cell_gamma.row(i) * RTM_CellMap[j].cell_mtx(i,0);
    }
  }
}


//TODO FIX HERE
//ALSO USE THE UPDATED ELBO CALCULATION
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
    ELBO_theta += arma::sum(lgamma(alpha.row(i) + n_dk.row(i))) - theta_denom;
  }


  // ELBO phi
  ELBO_phi_prior += (lgamma(arma::accu(beta)) - arma::accu(lgamma(beta)));

  for (int j = 0; j < K; j++){
    ELBO_phi_like += arma::accu(lgamma(beta.col(j) + n_wk.col(j)));
    ELBO_phi_norm += lgamma(arma::sum(beta.col(j)) + arma::sum(n_wk.col(j)));
  }

  ELBO_phi = (K*ELBO_phi_prior + ELBO_phi_like - ELBO_phi_norm);

  // ELBO gamma
  for (int cell = 1; cell <= D; cell++){
    tokens = RTM_CellMap[cell].cell_mtx.n_rows;
    for (int token = 0; token < tokens; token++){
      ELBO_gamma += arma::sum(RTM_CellMap[cell].cell_gamma.row(token) %
        arma::log(RTM_CellMap[cell].cell_gamma.row(token) + 0.0000001)*RTM_CellMap[cell].cell_mtx(token,0));
    }
  }
  //Rcout << ELBO_theta/C << ELBO_phi/C << ELBO_gamma/C << std::endl;
  arma::rowvec cur_val = {ELBO_theta/C,ELBO_phi/C,ELBO_gamma/C, (ELBO_theta + ELBO_phi -ELBO_gamma)/C};
  elbo.row(cur_iter).subvec(0,3) = cur_val;
  return (ELBO_theta + ELBO_phi + -ELBO_gamma)/C;
}



// Initialize + Iterate through weights
//////////////////////////////////////////


// int loss_fun 0,1,2 - exp_mean,cross_ent, softmax cross_ent
// [[Rcpp::export]]
arma::vec train_RTM(arma::sp_mat& counts,
                          arma::vec& celltypes,
                          arma::vec& genes,
                          Rcpp::List& nbr_list,
                          const arma::mat& alpha, const arma::mat& beta,
                          int K,int D,arma::mat& n_dk,arma::mat& n_wk,
                          int num_threads = 1,int maxiter = 100,
                          bool verbal = true,
                          bool zero_gamma = false,
                          bool rand_gamma = true,
                          double thresh = 0.00001,
                          double lr = 0.001,
                          double rho = 1000.0,
                          int loss_fun = 0,
                          bool m_update = true,
                          int burnin = 1){
  double old_elbo = 0;
  double cur_elbo = 0;
  arma::mat elbo_val = arma::zeros<arma::mat>(maxiter,5);
  int M = genes.n_elem;
  int C = arma::accu(counts);
  double prog = 0;
  int prog_width = 50;
  int pos = 0;
  arma::mat ndk_norm = arma::zeros<arma::mat>(D,K);
 if (verbal){
   Rcout << "Building RTM_CellMap" << std::endl;
 }

  std::unordered_map<int,RTM_Cell> RTM_CellMap = build_RTM_Cell_Map(counts,
                                                               celltypes,
                                                               genes,
                                                               alpha,
                                                               beta,
                                                               D,K,
                                                               zero_gamma,
                                                               rand_gamma,
                                                               nbr_list,
                                                               loss_fun);
  if (verbal){
    Rcout << "RTM_CellMap is built" << std::endl;
  }

  //Complete // Consider making the unordered_map its own object
  build_nwk_RTM(n_wk,RTM_CellMap,D,K);
  //Rcout << "here nwk" << std::endl;
  //double lr = 0.0001;
  //TComplete // Consider making the unordered_map its own object
  build_ndk_RTM(n_dk,RTM_CellMap,D,K);
  //Rcout << "here ndk" << std::endl;
  // Initialize Weights

  arma::vec model_weights = arma::ones(K);
  double model_bias = 0.0;
  double num_links = 0.0;
  for (int i = 1; i <= D; i++){
    num_links += RTM_CellMap[i].cell_nbrs.n_elem;
  }
  num_links /= 2;
  //TODO FIX HERE
  for (int i = 0; i < maxiter; i++){
    //E-Step
    //Rcout << "Start E-step" << std::endl;

    run_RTM_epoch(RTM_CellMap,
                  model_weights,
                  model_bias,
                  alpha,beta,K,D,n_dk,n_wk,
                  num_threads,
                  loss_fun,
                  burnin);
   //Rcout << "End E-step" << std::endl;
    //M-Step
    //Rcout << "here 8" << std::endl;
    build_nwk_RTM(n_wk,RTM_CellMap,D,K);
    //Rcout << "here 9" << std::endl;
    build_ndk_RTM(n_dk,RTM_CellMap,D,K);
    //Rcout << "here 10" << std::endl;
    // M STEP GOES HERE
    //Rcout << "Start M-step" << std::endl;
    if (m_update){
      update_weights(RTM_CellMap,
                     K,
                     D,
                     rho,
                     num_links,
                     n_dk,
                     model_bias,
                     model_weights,
                     loss_fun,
                     lr);
    }
    //Rcout << "End M-step" << std::endl;
    /////////
    //TODO ENDS here
    //Rcout << model_weights << std::endl;
    //Rcout << "ELBO Start" << std::endl;
    cur_elbo = get_RTM_ELBO(RTM_CellMap,
                        alpha,
                        beta,
                        M,D,K,
                        n_dk,
                        n_wk,
                        C,
                        elbo_val,
                        i);//(i+1);

    //Rcout << "ELBO Mid" << std::endl;
    //Add RTM ELBO Term -TODO As well
    double rtm_elbo = 0.0;
    //Rcout << "ELBO" << std::endl;
    if (loss_fun == 0){
      for (int j = 0; j < D; j++){
        arma::uvec nbr_ids = RTM_CellMap[j+1].cell_nbrs;
        if (nbr_ids.n_elem == 0){
          continue;
        }
        arma::mat sub_mat = ndk_norm.rows(nbr_ids);
        // Rcout << size(sub_mat) << std::endl;
        // Rcout << size((ndk_norm.row(j) % model_weights.t())) << std::endl;
        rtm_elbo += arma::accu(arma::sum(sub_mat.each_row() %
          (ndk_norm.row(j) % model_weights.t()),0) + model_bias);
      }
      //TODO ENDS here
      //Rcout << "here 2" << std::endl;
    } else{
      // use the nbr_pred function to make this faster
      //create an adjacency matrix based on inputs

      for (int j = 0; j < D; j++){
        arma::uvec nbr_ids = RTM_CellMap[j+1].cell_nbrs;
        if (nbr_ids.n_elem == 0){
          continue;
        }

        arma::mat sub_mat = ndk_norm.rows(nbr_ids);
        arma::vec cur_out = sub_mat*model_weights + model_bias;
        if (loss_fun == 1){
          arma::vec nbr_bool = RTM_CellMap[j+1].bool_nbrs;
          cur_out = 1/(1+arma::exp(-cur_out));
          rtm_elbo += arma::accu(nbr_bool % arma::log1p(cur_out) + (1-nbr_bool) % arma::log1p(1-cur_out));
        } else if ( loss_fun == 2) {
          arma::vec nbr_dist = RTM_CellMap[j+1].nbr_dist;
          arma::vec pdf_values = arma::exp(-arma::square(nbr_dist - cur_out) / 2.0) / std::sqrt(2.0 * M_PI);
          rtm_elbo += arma::accu(arma::log1p(pdf_values));                 // Sum log1p results
        }
      }
    }
    rtm_elbo /= (2*num_links);
    elbo_val(i,4) = rtm_elbo;
    //Rcout << "ELBO End" << std::endl;
    cur_elbo += rtm_elbo;
    //Rcout << "ELBO Done" << std::endl;
    old_elbo = cur_elbo;
    //TODO ENDS here

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

  }


  if (verbal){
    Rcout.flush();
    Rcout << "[";
    for(int p = 0; p < prog_width; p++){
      Rcout << "=" ;
    }
    Rcout << "] " << "100% || Iter: " << maxiter <<  " || ELBO: "<< cur_elbo << std::endl;
    Rcout << "Max Iteration Reached" << std::endl;
  }
  //elbo_val.save("elbo_out.csv",arma::csv_ascii);
  //Rcout << model_bias << std::endl;
  //Rcout << model_weights << std::endl;
  //Combine model_bias and model_weights
  arma::vec model_weights_out = arma::zeros<arma::vec>(K+1);
  model_weights_out(0) = model_bias;
  model_weights_out.subvec(1,K) = model_weights;
  RTM_CellMap.clear();
  return model_weights_out;
}



// [[Rcpp::export]]
arma::mat nbr_pred(arma::mat theta,arma::rowvec weights,double max_val = 1,int loss_fun =0){
  arma::mat pred(theta.n_rows,theta.n_rows);
  for (int i = 0; i < theta.n_rows; i++){
    for (int j = 0; j < theta.n_rows; j++){
      arma::rowvec cur_ndk = theta.row(i);
      arma::rowvec nbr_ndk = theta.row(j);
      cur_ndk = cur_ndk % nbr_ndk;
      cur_ndk.insert_cols(0,1); //Adding bias term
      if (loss_fun == 0){
        double pred_val = arma::as_scalar(arma::exp(cur_ndk*weights.t()));
        pred(i,j) = pred_val/max_val;
      }
      else if (loss_fun == 1){
        double pred_val = arma::as_scalar(arma::exp(-cur_ndk*weights.t()));
        pred(i,j) = 1/(1+pred_val);
      } else if (loss_fun == 2){
        pred(i,j) = arma::as_scalar(cur_ndk*weights.t());
      }
    }
  }
  return(pred);
}



// [[Rcpp::export]]
NumericMatrix get_dist_cpp(NumericMatrix x) {
  int n = x.nrow();
  NumericMatrix out(n,n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      out(i,j) = sqrt(pow(x(i,0) - x(j,0),2) + pow(x(i,1) - x(j,1),2));
    }
  }
  return out;
}

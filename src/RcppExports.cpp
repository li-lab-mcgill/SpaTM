// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// test_cell_map
void test_cell_map(arma::sp_mat& counts, arma::vec& celltypes, arma::vec& genes, const arma::mat& alpha, const arma::mat& beta, int D, int K, bool zero_gamma, bool rand_gamma, int cellid);
RcppExport SEXP _SpaTM_test_cell_map(SEXP countsSEXP, SEXP celltypesSEXP, SEXP genesSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP DSEXP, SEXP KSEXP, SEXP zero_gammaSEXP, SEXP rand_gammaSEXP, SEXP cellidSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type celltypes(celltypesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< bool >::type zero_gamma(zero_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type rand_gamma(rand_gammaSEXP);
    Rcpp::traits::input_parameter< int >::type cellid(cellidSEXP);
    test_cell_map(counts, celltypes, genes, alpha, beta, D, K, zero_gamma, rand_gamma, cellid);
    return R_NilValue;
END_RCPP
}
// test_RTM_cell
void test_RTM_cell();
RcppExport SEXP _SpaTM_test_RTM_cell() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test_RTM_cell();
    return R_NilValue;
END_RCPP
}
// test_RTM_cellmap
void test_RTM_cellmap(arma::sp_mat& counts, arma::vec& celltypes, arma::vec& genes, const arma::mat& alpha, const arma::mat& beta, int D, int K, bool zero_gamma, bool rand_gamma, Rcpp::List& nbr_list, int print_cell, int loss_fun);
RcppExport SEXP _SpaTM_test_RTM_cellmap(SEXP countsSEXP, SEXP celltypesSEXP, SEXP genesSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP DSEXP, SEXP KSEXP, SEXP zero_gammaSEXP, SEXP rand_gammaSEXP, SEXP nbr_listSEXP, SEXP print_cellSEXP, SEXP loss_funSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type celltypes(celltypesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< bool >::type zero_gamma(zero_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type rand_gamma(rand_gammaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type nbr_list(nbr_listSEXP);
    Rcpp::traits::input_parameter< int >::type print_cell(print_cellSEXP);
    Rcpp::traits::input_parameter< int >::type loss_fun(loss_funSEXP);
    test_RTM_cellmap(counts, celltypes, genes, alpha, beta, D, K, zero_gamma, rand_gamma, nbr_list, print_cell, loss_fun);
    return R_NilValue;
END_RCPP
}
// train_RTM
arma::vec train_RTM(arma::sp_mat& counts, arma::vec& celltypes, arma::vec& genes, Rcpp::List& nbr_list, const arma::mat& alpha, const arma::mat& beta, int K, int D, arma::mat& n_dk, arma::mat& n_wk, int num_threads, int maxiter, bool verbal, bool zero_gamma, bool rand_gamma, double thresh, double lr, double rho, int loss_fun, bool m_update);
RcppExport SEXP _SpaTM_train_RTM(SEXP countsSEXP, SEXP celltypesSEXP, SEXP genesSEXP, SEXP nbr_listSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP KSEXP, SEXP DSEXP, SEXP n_dkSEXP, SEXP n_wkSEXP, SEXP num_threadsSEXP, SEXP maxiterSEXP, SEXP verbalSEXP, SEXP zero_gammaSEXP, SEXP rand_gammaSEXP, SEXP threshSEXP, SEXP lrSEXP, SEXP rhoSEXP, SEXP loss_funSEXP, SEXP m_updateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type celltypes(celltypesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type nbr_list(nbr_listSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type n_dk(n_dkSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type n_wk(n_wkSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbal(verbalSEXP);
    Rcpp::traits::input_parameter< bool >::type zero_gamma(zero_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type rand_gamma(rand_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< double >::type lr(lrSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type loss_fun(loss_funSEXP);
    Rcpp::traits::input_parameter< bool >::type m_update(m_updateSEXP);
    rcpp_result_gen = Rcpp::wrap(train_RTM(counts, celltypes, genes, nbr_list, alpha, beta, K, D, n_dk, n_wk, num_threads, maxiter, verbal, zero_gamma, rand_gamma, thresh, lr, rho, loss_fun, m_update));
    return rcpp_result_gen;
END_RCPP
}
// nbr_pred
arma::mat nbr_pred(arma::mat theta, arma::rowvec weights, double max_val, int loss_fun);
RcppExport SEXP _SpaTM_nbr_pred(SEXP thetaSEXP, SEXP weightsSEXP, SEXP max_valSEXP, SEXP loss_funSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type max_val(max_valSEXP);
    Rcpp::traits::input_parameter< int >::type loss_fun(loss_funSEXP);
    rcpp_result_gen = Rcpp::wrap(nbr_pred(theta, weights, max_val, loss_fun));
    return rcpp_result_gen;
END_RCPP
}
// get_CE
double get_CE(arma::mat pred_adj, arma::mat grd_adj);
RcppExport SEXP _SpaTM_get_CE(SEXP pred_adjSEXP, SEXP grd_adjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type pred_adj(pred_adjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type grd_adj(grd_adjSEXP);
    rcpp_result_gen = Rcpp::wrap(get_CE(pred_adj, grd_adj));
    return rcpp_result_gen;
END_RCPP
}
// get_acc
double get_acc(arma::mat pred_adj, arma::mat grd_adj, double thresh);
RcppExport SEXP _SpaTM_get_acc(SEXP pred_adjSEXP, SEXP grd_adjSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type pred_adj(pred_adjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type grd_adj(grd_adjSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(get_acc(pred_adj, grd_adj, thresh));
    return rcpp_result_gen;
END_RCPP
}
// get_dist_cpp
NumericMatrix get_dist_cpp(NumericMatrix x);
RcppExport SEXP _SpaTM_get_dist_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(get_dist_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// rtm_mlp_forward
arma::rowvec rtm_mlp_forward(const arma::mat& X, int layers, Rcpp::List& weights);
RcppExport SEXP _SpaTM_rtm_mlp_forward(SEXP XSEXP, SEXP layersSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type layers(layersSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(rtm_mlp_forward(X, layers, weights));
    return rcpp_result_gen;
END_RCPP
}
// mlp_nbr_pred
arma::mat mlp_nbr_pred(const arma::mat& X, int layers, Rcpp::List weights);
RcppExport SEXP _SpaTM_mlp_nbr_pred(SEXP XSEXP, SEXP layersSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type layers(layersSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(mlp_nbr_pred(X, layers, weights));
    return rcpp_result_gen;
END_RCPP
}
// test_stm_cell
void test_stm_cell();
RcppExport SEXP _SpaTM_test_stm_cell() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test_stm_cell();
    return R_NilValue;
END_RCPP
}
// test_stm_cellmap
void test_stm_cellmap(arma::sp_mat& counts, arma::vec& celltypes, arma::vec& genes, const arma::mat& alpha, const arma::mat& beta, int D, int K, bool zero_gamma, bool rand_gamma, arma::vec& labels, int print_cell);
RcppExport SEXP _SpaTM_test_stm_cellmap(SEXP countsSEXP, SEXP celltypesSEXP, SEXP genesSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP DSEXP, SEXP KSEXP, SEXP zero_gammaSEXP, SEXP rand_gammaSEXP, SEXP labelsSEXP, SEXP print_cellSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type celltypes(celltypesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< bool >::type zero_gamma(zero_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type rand_gamma(rand_gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< int >::type print_cell(print_cellSEXP);
    test_stm_cellmap(counts, celltypes, genes, alpha, beta, D, K, zero_gamma, rand_gamma, labels, print_cell);
    return R_NilValue;
END_RCPP
}
// train_stm_scgtm
arma::mat train_stm_scgtm(arma::sp_mat& counts, arma::vec& celltypes, arma::vec& genes, arma::vec& labels, const arma::mat& alpha, const arma::mat& beta, int K, int D, arma::mat& n_dk, arma::mat& n_wk, int num_threads, int maxiter, bool verbal, bool zero_gamma, bool rand_gamma, double thresh, double lr);
RcppExport SEXP _SpaTM_train_stm_scgtm(SEXP countsSEXP, SEXP celltypesSEXP, SEXP genesSEXP, SEXP labelsSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP KSEXP, SEXP DSEXP, SEXP n_dkSEXP, SEXP n_wkSEXP, SEXP num_threadsSEXP, SEXP maxiterSEXP, SEXP verbalSEXP, SEXP zero_gammaSEXP, SEXP rand_gammaSEXP, SEXP threshSEXP, SEXP lrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type celltypes(celltypesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type n_dk(n_dkSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type n_wk(n_wkSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbal(verbalSEXP);
    Rcpp::traits::input_parameter< bool >::type zero_gamma(zero_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type rand_gamma(rand_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< double >::type lr(lrSEXP);
    rcpp_result_gen = Rcpp::wrap(train_stm_scgtm(counts, celltypes, genes, labels, alpha, beta, K, D, n_dk, n_wk, num_threads, maxiter, verbal, zero_gamma, rand_gamma, thresh, lr));
    return rcpp_result_gen;
END_RCPP
}
// progress_bar
void progress_bar(double iter, int maxiter, double elbo, double train_acc, double val_acc);
RcppExport SEXP _SpaTM_progress_bar(SEXP iterSEXP, SEXP maxiterSEXP, SEXP elboSEXP, SEXP train_accSEXP, SEXP val_accSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type elbo(elboSEXP);
    Rcpp::traits::input_parameter< double >::type train_acc(train_accSEXP);
    Rcpp::traits::input_parameter< double >::type val_acc(val_accSEXP);
    progress_bar(iter, maxiter, elbo, train_acc, val_acc);
    return R_NilValue;
END_RCPP
}
// mlp_forward
arma::mat mlp_forward(const arma::mat& X, int layers, Rcpp::List& weights, bool dummy_topic);
RcppExport SEXP _SpaTM_mlp_forward(SEXP XSEXP, SEXP layersSEXP, SEXP weightsSEXP, SEXP dummy_topicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type layers(layersSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type dummy_topic(dummy_topicSEXP);
    rcpp_result_gen = Rcpp::wrap(mlp_forward(X, layers, weights, dummy_topic));
    return rcpp_result_gen;
END_RCPP
}
// get_label_prob
arma::rowvec get_label_prob(arma::rowvec ndk_row, const arma::rowvec& gamma_weight, int gene_count, int K, int label, int layers, Rcpp::List& weights, bool dummy_topic);
RcppExport SEXP _SpaTM_get_label_prob(SEXP ndk_rowSEXP, SEXP gamma_weightSEXP, SEXP gene_countSEXP, SEXP KSEXP, SEXP labelSEXP, SEXP layersSEXP, SEXP weightsSEXP, SEXP dummy_topicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type ndk_row(ndk_rowSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type gamma_weight(gamma_weightSEXP);
    Rcpp::traits::input_parameter< int >::type gene_count(gene_countSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type label(labelSEXP);
    Rcpp::traits::input_parameter< int >::type layers(layersSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type dummy_topic(dummy_topicSEXP);
    rcpp_result_gen = Rcpp::wrap(get_label_prob(ndk_row, gamma_weight, gene_count, K, label, layers, weights, dummy_topic));
    return rcpp_result_gen;
END_RCPP
}
// stm_torch_estep
double stm_torch_estep(arma::sp_mat& spe_counts, arma::vec& celltypes, arma::vec& genes, const arma::mat& alpha, const arma::mat& beta, int K, int D, bool zero_gamma, bool rand_gamma, arma::vec& labels, arma::mat& n_dk, arma::mat& n_wk, arma::mat& elbo_mat, int layers, Rcpp::List& weights, double& elbo, int num_threads, int iter, bool dummy_topic);
RcppExport SEXP _SpaTM_stm_torch_estep(SEXP spe_countsSEXP, SEXP celltypesSEXP, SEXP genesSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP KSEXP, SEXP DSEXP, SEXP zero_gammaSEXP, SEXP rand_gammaSEXP, SEXP labelsSEXP, SEXP n_dkSEXP, SEXP n_wkSEXP, SEXP elbo_matSEXP, SEXP layersSEXP, SEXP weightsSEXP, SEXP elboSEXP, SEXP num_threadsSEXP, SEXP iterSEXP, SEXP dummy_topicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type spe_counts(spe_countsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type celltypes(celltypesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< bool >::type zero_gamma(zero_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type rand_gamma(rand_gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type n_dk(n_dkSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type n_wk(n_wkSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type elbo_mat(elbo_matSEXP);
    Rcpp::traits::input_parameter< int >::type layers(layersSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double& >::type elbo(elboSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< bool >::type dummy_topic(dummy_topicSEXP);
    rcpp_result_gen = Rcpp::wrap(stm_torch_estep(spe_counts, celltypes, genes, alpha, beta, K, D, zero_gamma, rand_gamma, labels, n_dk, n_wk, elbo_mat, layers, weights, elbo, num_threads, iter, dummy_topic));
    return rcpp_result_gen;
END_RCPP
}
// get_theta
arma::mat get_theta(const arma::mat& n_dk, const arma::mat& alpha);
RcppExport SEXP _SpaTM_get_theta(SEXP n_dkSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type n_dk(n_dkSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(get_theta(n_dk, alpha));
    return rcpp_result_gen;
END_RCPP
}
// get_phi
arma::mat get_phi(const arma::mat& n_wk, const arma::mat& beta);
RcppExport SEXP _SpaTM_get_phi(SEXP n_wkSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type n_wk(n_wkSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(get_phi(n_wk, beta));
    return rcpp_result_gen;
END_RCPP
}
// train_scgtm
void train_scgtm(arma::sp_mat& counts, arma::vec& celltypes, arma::vec& genes, const arma::mat& alpha, const arma::mat& beta, int K, int D, arma::mat& n_dk, arma::mat& n_wk, int num_threads, int maxiter, bool verbal, bool zero_gamma, bool rand_gamma, double thresh);
RcppExport SEXP _SpaTM_train_scgtm(SEXP countsSEXP, SEXP celltypesSEXP, SEXP genesSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP KSEXP, SEXP DSEXP, SEXP n_dkSEXP, SEXP n_wkSEXP, SEXP num_threadsSEXP, SEXP maxiterSEXP, SEXP verbalSEXP, SEXP zero_gammaSEXP, SEXP rand_gammaSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type celltypes(celltypesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type n_dk(n_dkSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type n_wk(n_wkSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbal(verbalSEXP);
    Rcpp::traits::input_parameter< bool >::type zero_gamma(zero_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type rand_gamma(rand_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    train_scgtm(counts, celltypes, genes, alpha, beta, K, D, n_dk, n_wk, num_threads, maxiter, verbal, zero_gamma, rand_gamma, thresh);
    return R_NilValue;
END_RCPP
}
// predict_scgtm
void predict_scgtm(arma::sp_mat& counts, arma::vec& celltypes, arma::vec& genes, double& alpha, int K, int D, int M, arma::mat& n_dk, const arma::mat& phi, int num_threads, int maxiter, bool verbal);
RcppExport SEXP _SpaTM_predict_scgtm(SEXP countsSEXP, SEXP celltypesSEXP, SEXP genesSEXP, SEXP alphaSEXP, SEXP KSEXP, SEXP DSEXP, SEXP MSEXP, SEXP n_dkSEXP, SEXP phiSEXP, SEXP num_threadsSEXP, SEXP maxiterSEXP, SEXP verbalSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type celltypes(celltypesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type n_dk(n_dkSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbal(verbalSEXP);
    predict_scgtm(counts, celltypes, genes, alpha, K, D, M, n_dk, phi, num_threads, maxiter, verbal);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SpaTM_test_cell_map", (DL_FUNC) &_SpaTM_test_cell_map, 10},
    {"_SpaTM_test_RTM_cell", (DL_FUNC) &_SpaTM_test_RTM_cell, 0},
    {"_SpaTM_test_RTM_cellmap", (DL_FUNC) &_SpaTM_test_RTM_cellmap, 12},
    {"_SpaTM_train_RTM", (DL_FUNC) &_SpaTM_train_RTM, 20},
    {"_SpaTM_nbr_pred", (DL_FUNC) &_SpaTM_nbr_pred, 4},
    {"_SpaTM_get_CE", (DL_FUNC) &_SpaTM_get_CE, 2},
    {"_SpaTM_get_acc", (DL_FUNC) &_SpaTM_get_acc, 3},
    {"_SpaTM_get_dist_cpp", (DL_FUNC) &_SpaTM_get_dist_cpp, 1},
    {"_SpaTM_rtm_mlp_forward", (DL_FUNC) &_SpaTM_rtm_mlp_forward, 3},
    {"_SpaTM_mlp_nbr_pred", (DL_FUNC) &_SpaTM_mlp_nbr_pred, 3},
    {"_SpaTM_test_stm_cell", (DL_FUNC) &_SpaTM_test_stm_cell, 0},
    {"_SpaTM_test_stm_cellmap", (DL_FUNC) &_SpaTM_test_stm_cellmap, 11},
    {"_SpaTM_train_stm_scgtm", (DL_FUNC) &_SpaTM_train_stm_scgtm, 17},
    {"_SpaTM_progress_bar", (DL_FUNC) &_SpaTM_progress_bar, 5},
    {"_SpaTM_mlp_forward", (DL_FUNC) &_SpaTM_mlp_forward, 4},
    {"_SpaTM_get_label_prob", (DL_FUNC) &_SpaTM_get_label_prob, 8},
    {"_SpaTM_stm_torch_estep", (DL_FUNC) &_SpaTM_stm_torch_estep, 19},
    {"_SpaTM_get_theta", (DL_FUNC) &_SpaTM_get_theta, 2},
    {"_SpaTM_get_phi", (DL_FUNC) &_SpaTM_get_phi, 2},
    {"_SpaTM_train_scgtm", (DL_FUNC) &_SpaTM_train_scgtm, 15},
    {"_SpaTM_predict_scgtm", (DL_FUNC) &_SpaTM_predict_scgtm, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_SpaTM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

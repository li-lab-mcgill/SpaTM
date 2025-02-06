#' Spatial Topic Model (SpaTM)
#'
#' @description
#' A topic modelling framework for Spatial and Single cell analyses.
#' SpaTM employs three topic models:
#'
#' 1. The Guided Topic Model for anchoring topics to cell types or other features
#'
#' 2. A Supervised Topic Model for predicting metadata features
#'
#' 3. A Relational Topic Model for predicting neighboring cells/spots
#'
#' @author Adrien Osakwe
#'
#' @name SpaTM-pkg
#' @useDynLib SpaTM, .registration=TRUE
#' @import Rcpp RcppArmadillo
NULL

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STM.R
\name{STM}
\alias{STM}
\title{Spatial Topic Modeling (STM) for Single-Cell Spatial Data}
\usage{
STM(
  spe,
  label = NULL,
  num_threads = 1L,
  maxiter = 100L,
  verbal = TRUE,
  zero_gamma = TRUE,
  rand_gamma = FALSE,
  thresh = 1e-05,
  lr = 1e-04,
  mlp = NULL,
  mlp_layers = NULL,
  mlp_epoch = NULL,
  spe_val = NULL,
  device = NULL,
  balanced_class = TRUE,
  dummy_topic = FALSE,
  burnin = 5
)
}
\arguments{
\item{spe}{A spatial experiment object (typically of class `SummarizedExperiment`) containing spatial transcriptomics data, including counts, metadata, and gene information.}

\item{num_threads}{Number of threads to use for parallel computation (default is 1).}

\item{maxiter}{Maximum number of iterations for model training (default is 100).}

\item{verbal}{Logical indicating whether to display progress messages (default is `TRUE`).}

\item{zero_gamma}{Logical indicating whether to initialize the gamma parameters to zero (default is `TRUE`).}

\item{rand_gamma}{Logical indicating whether to randomize the gamma parameters (default is `FALSE`).}

\item{thresh}{Convergence threshold for stopping the model training (default is 0.00001).}

\item{lr}{Learning rate for optimization (default is 0.0001).}

\item{mlp}{Optional, a model definition for the multi-layer perceptron (MLP) when running Torch-STM. If `NULL`, GLM-STM is run.}

\item{mlp_layers}{Optional, a vector specifying the number of units in each layer of the MLP.}

\item{mlp_epoch}{Optional, the number of epochs for training the MLP.}

\item{spe_val}{A validation set (typically of class `SummarizedExperiment`) used for monitoring overfitting in Torch-STM. Must be provided for Torch-STM.}

\item{device}{The device for training (e.g., "cpu" or "cuda") when using Torch-STM. If `NULL`, the default device is used.}

\item{balanced_class}{Logical indicating whether to balance the class distribution during model training (default is `TRUE`).}

\item{dummy_topic}{Logical indicating whether to add a dummy topic for regularization (default is `FALSE`).}

\item{burnin}{Number of iterations for burn-in phase in Torch-STM (default is 5).}
}
\value{
The input spatial experiment object (`spe`) with updated metadata containing the STM model results:
  - `STM_Weights`: A matrix of learned topic weights (GLM-STM).
  - `STM_MLP`: The trained MLP model (Torch-STM).
}
\description{
This function performs spatial topic modeling on single-cell spatial transcriptomics data. It offers two modes:
1. GLM-STM (default): A generalized linear model for topic modeling.
2. Torch-STM: A deep learning-based STM model using a multi-layer perceptron (MLP).
}
\details{
The function supports various parameters for training, including regularization, optimization settings, and the option for a validation set to prevent overfitting in the Torch-STM mode.
}
\examples{
# Example 1: Running GLM-STM with default settings
spe_result <- STM(spe = spe_object)

# Example 2: Running Torch-STM with MLP model
spe_result <- STM(spe = spe_object, mlp = my_mlp_model, mlp_layers = c(64, 32), spe_val = spe_validation)

}

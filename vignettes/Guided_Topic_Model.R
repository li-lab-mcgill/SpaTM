## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE,warning=F-------------------------------------------
library(SpaTM)
library(SingleCellExperiment)

## ----simulate_data------------------------------------------------------------
set.seed(42)

# Parameters for simulation
n_genes <- 100
n_cells_per_type <- 50
n_cell_types <- 4
total_cells <- n_cells_per_type * n_cell_types

# Create base expression patterns for each cell type
base_patterns <- matrix(
  rgamma(n_genes * n_cell_types, shape = 2, scale = 1),
  nrow = n_genes,
  ncol = n_cell_types
)

# Create empty count matrix
counts <- matrix(0, nrow = n_genes, ncol = total_cells)
colnames(counts) <- paste0("cell_", 1:total_cells)
rownames(counts) <- paste0("gene_", 1:n_genes)

# Generate cell type labels
cell_types <- rep(paste0("celltype_", 1:n_cell_types), each = n_cells_per_type)

# Fill count matrix with simulated data
for (i in 1:n_cell_types) {
  start_idx <- ((i-1) * n_cells_per_type) + 1
  end_idx <- i * n_cells_per_type
  
  # Add noise and generate counts
  counts[, start_idx:end_idx] <- sapply(1:n_cells_per_type, function(x) {
    rpois(n_genes, lambda = base_patterns[, i] * runif(n_genes, 0.8, 1.2))
  })
}

# Create SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = DataFrame(cell_type = cell_types)
)

# Display dimensions and structure
sce
as.matrix(table(sce$cell_type))/ncol(sce)

## -----------------------------------------------------------------------------
#Converting to a SingleCellTopicExperiment
sce <- SingleCellTopicExperiment(sce,
                                 guided = T,
                                 labels = 'cell_type')
K <- ncol(alphaPrior(sce))
D <- ncol(sce)
sce

## -----------------------------------------------------------------------------
sce <- GTM(sce,K,D,num_threads = 1,maxiter = 50,
           verbal = F,zero_gamma = F,rand_gamma = F)

#Calculating parameter estimates
sce <- buildTheta(sce)
sce <- buildPhi(sce)

#Check training performance
sce$Pred <- levels(sce$cell_type)[apply(theta(sce),1,which.max)]
table(sce$cell_type,sce$Pred)

## -----------------------------------------------------------------------------

## Creating pseudobulk of training dataset
bulk_mat <- as.matrix(rowSums(counts))
bulk_se <- SingleCellExperiment(assays = list(counts = bulk_mat))
bulk_se <- SingleCellTopicExperiment(bulk_se,
                                     K = K)
## Running Topic Inference
bulk_se <- inferTopics(bulk_se,maxiter = 100,phi = phi(sce),verbal = F)
bulk_se <- buildTheta(bulk_se)
colnames(theta(bulk_se)) <- colnames(phi(sce))

##Printing Deconvolution Results
theta(bulk_se)

## -----------------------------------------------------------------------------
sessionInfo()


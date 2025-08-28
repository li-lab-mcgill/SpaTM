# (Archive) A Spatial Topic Model (SpaTM) for Inferring Spatially Informed Transcriptional Programs

**NOTE: This branch was made to archive the exact version of SpaTM used to generate the results in our current preprint as of July 29 2025.** An updated version of the code that has been adapted to improve ease-of-use is being updated into the **main branch**.

The corresponding version of the manuscript [can be found here](https://www.biorxiv.org/content/10.1101/2025.01.24.634726v2).

The corresponding analysis scripts [can be found here](https://github.com/aosakwe/SpaTM_Analysis).


## Details
This repository contains the R package for the Spatial Topic Model we have developed for flexible and interpretable analysis of spatial transcriptomics data. The method integrates multiple topic modelling frameworks to enable:
- Inferrence of transcriptional programs that can predict histology-based annotations
- Identification of spatially informed clusters through neighbor prediction tasks
- Cell type Deconvolution (with or without a reference single cell dataset)
- Impute spatial information in single cell atlases

## Installation
```
devtools::install_github("li-lab-mcgill/SpaTM", ref = "archive/v0.3")
```

## Issues
Please use the [issues](https://github.com/li-lab-mcgill/SpaTM/issues) tab to contact us if you encounter any issues/bugs or have suggestions for our software!


## SpaTM Workflow
![SpaTM Workflow](https://github.com/aosakwe/SpaTM_Analysis/blob/main/SpaTM.png)

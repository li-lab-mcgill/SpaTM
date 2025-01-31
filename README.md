# A Spatial Topic Model (SpaTM) for Inferring Spatially Informed Transcriptional Programs

This repository contains the R package for the Spatial Topic Model we have developed for flexible and interpretable analysis of spatial transcriptomics data. The method integrates multiple topic modelling frameworks to enable:

- Inferrence of transcriptional programs that can predict histology-based annotations
- Identification of spatially informed clusters through neighbor prediction tasks
- Cell type Deconvolution (with or without a reference single cell dataset)
- Impute spatial information in single cell atlases

Our manuscript pre-print detailing the method's performance and analyses [can be found here](https://www.biorxiv.org/content/10.1101/2025.01.24.634726v1)

## Installation
```
devtools::install_github("li-lab-mcgill/SpaTM")
```
## Issues
Please use the [issues](https://github.com/li-lab-mcgill/SpaTM/issues) tab to contact us if you encounter any issues/bugs or have suggestions for our software!


## SpaTM Workflow
![SpaTM Workflow](https://github.com/aosakwe/SpaTM_Analysis/blob/main/SpaTM.png)

# A Spatial Topic Model (SpaTM) for Inferring Spatially Informed Transcriptional Programs

This repository contains the R package for the Spatial Topic Model we have developed for flexible and interpretable analysis of spatial transcriptomics data. The method integrates multiple topic modelling frameworks to enable:

- Inferrence of transcriptional programs that can predict histology-based annotations
- Identification of spatially informed clusters through neighbor prediction tasks
- Cell type Deconvolution (with or without a reference single cell dataset)
- Impute spatial information in single cell atlases

Our manuscript (now published in _Briefings in Bioinformatics_) detailing the method's performance and analyses [can be found here](https://academic.oup.com/bib/article/26/6/bbaf657/8374204). The analysis scripts used to generate the results found in the manuscript (using [version 0.3 of SpaTM](https://github.com/li-lab-mcgill/SpaTM/tree/archive/v0.3)) [can be found here](https://github.com/aosakwe/SpaTM_Analysis).

## Installation
```
devtools::install_github("li-lab-mcgill/SpaTM")
```
## Issues
Please use the [issues](https://github.com/li-lab-mcgill/SpaTM/issues) tab to contact us if you encounter any issues/bugs or have suggestions for our software!


## SpaTM Workflow
![SpaTM Workflow](https://github.com/aosakwe/SpaTM_Analysis/blob/main/SpaTM.png)

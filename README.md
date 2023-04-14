# DESpace_manuscript
This repository contains the code used for the analyses presented in the manuscript "DESpace: spatially variable gene detection via differential expression testing of spatial clusters".

The scripts for all analysis are organized in *./Analysis* folders:

- `01_preprocessing`: to do quality control and filtering;
- `02_run_methods`: to split mouse cerebellum *sce* to count matrix and meta data and scripts to run spatially variable gene detection methods in python;
- `03_simulations`: to get simulated data sets and R results used in the paper;
- `04_multiple_samples`: to get multiple samples simulated data sets and R results used in the paper;
- `05_individual_clusters`: to get individual clusters results;
- `06_real_data`: to run all R SVGs detection method based on real data.

The scripts for all plots are organized in *./Figures/Scripts* folders:

- `Scripts` contains the code to make all Figures available in the manuscript and its Supplementary;
- `Figures` contains all Figures shown in the manuscript and its Supplementary.

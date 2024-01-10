# DESpace_manuscript
This repository contains the code used for the analyses presented in the manuscript “DESpace: spatially variable gene detection via differential expression testing of spatial clusters”.

Scripts in *./Analyses* folder are organized as follows:

- `01_preprocessing`: quality control and filtering;
- `02_run_methods`: transform mouse cerebellum *sce* into count matrix and meta data, and run spatially variable gene detection methods in python;
- `03_simulations`: generate simulated data sets and corresponding analyses in R;
- `04_multiple_samples`: generate the multiple sample simulated data sets and corresponding analyses in R;
- `05_individual_clusters`: obtain results for the individual cluster test;
- `06_real_data`: run all SVG methods on real data in R.

Files in *./Figures* are arranged as follows:

- `Scripts` contains the code to make all Figures available in the manuscript and its Supplementary;
- `Figures` contains all Figures shown in the manuscript and its Supplementary.

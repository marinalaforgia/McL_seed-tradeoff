# McL_seed-tradeoff

Repository for data and code associated with McLaughlin HMM paper

-   **Data**: all data and output from models
    -   **McL-community-PA-data_annuals_2019.csv**: presence absence of annual species from 2000-2019 in 400 plots at McLaughlin Natural Reserve
    -   **HMM-meta-mcl.csv**: metadata on McLaughlin community data
    -   **McL_seed-traits.csv**: seed trait data
    -   **seed-trait-metadata.csv**: metadata on columns found in McL_seed-traits.csv
    -   **McL-species-boot.RDS**: HMM output from *01_HMM-McL-boot.R*
    -   **McL-sim-species-validation-50.RDS**: output of *02_HMM-McL-sim.R* using 50 patches (min number of patches a species can occur in)
    -   **McL-sim-154-correlation.RDS**: output of of *03_HMM-McL-cor.R* using simulated data from 154 patches (average number of patches a species occurred in)
-   **Scripts**: all scripts
    -   **00_HMM_Functions.R**: Creation of HMM functions
    -   **00_HMM_Simulation.R**: Creation of function to simulate time series data for a single species
    -   **01_HMM-McL-boot.R**: Code to run HMM functions on McLaughlin Data
    -   **02_HMM-McL-sim.R**: Code for checking accuracy of model estimates per species (produces Fig S3, corresponds to Pluntz et al. Fig A4)
    -   **03_HMM-McL-cor.R**: Code for checking correlation of data vs correlation from estimation procedure (produces Fig S2, corresponds to Pluntz et al. Fig A6)
    -   **04_HMM-McL-analysis.R**: Analyses and figures for paper

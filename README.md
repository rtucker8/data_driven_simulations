Simulation Code for Gonzalez, Dempsey, and Boonstra (2026+)
================

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Reproducibility](#reproducibility)
- [Reproducing the Simulations](#reproducing-the-simulations)
- [Reproducing the Figures](#reproducing-the-figures)
- [Computational Requirements](#computational-requirements)
- [Citation](#citation)

# Overview

This repository contains the R code and supporting files used to
generate the simulation results presented in:

> Gonzalez, Dempsey, and Boonstra. (2026+). *Estimation of State
> Occupation Probabilities for the Illness-Death Model via Multiple
> Imputation*. \[in preparation\].

The simulations evaluate the MSMI-KM and MSMI-Cox methods proposed in
the paper. Specifically, we examine the bias and coverage probabilities
of the proposed state occupation probability estimators and compare
their performance to the classical Aalen-Johansen estimator.

The repository is organized to facilitate reproduction of the
simulations, figures, and tables reported in the paper.

For questions about the code or simulations, please contact Rachel
Gonzalez, rbtucker@umich.edu.

# Repository Structure

``` text
.
├── R/                  # R functions used by the simulations
├── Code/               # Code for running simulations and analyses
├── Data/               # data from TCGA-CDR and random seeds
├── Results/            # Simulation output
  ├── lihc/               # Simulation setting based on lihc TCGA-CDR data
  ├── meso/               # Simulation setting based on meso TCGA-CDR data
  ├── prad/               # Simulation setting based on prad TCGA-CDR data
├── results_precomputed/# Precomputed simulation output
  ├── lihc/               # Simulation setting based on lihc TCGA-CDR data
  ├── meso/               # Simulation setting based on meso TCGA-CDR data
  ├── prad/               # Simulation setting based on prad TCGA-CDR data
├── figures/            # Figures generated from simulation results
├── README.Rmd          # Source file for this README
├── README.md           # GitHub README
└── renv.lock           # R package environment
```

# Requirements

The required R packages can be installed with:

``` r
install.packages(c(
  #data wrangling
  "tidyverse",
  "readxl",
  #statistical packages
  "mstate",
  "PWEXP",
  "survival",
  "sp",
  #visualization
  "patchwork",
  "ggtern",
  #other
  "devtools"
))

#proposed method
install_github("rtucker8/msmi")
```

# Reproducibility

All analyses reported in the paper were conducted using the code
contained in this repository.

To facilitate reproducibility:

- Random seeds are explicitly specified
- Simulation parameters are derived from the publicly available TCGA-CDR
  data using a process documented in the repository
- R package versions are recorded using `renv`
- Figures and tables are generated directly from simulation output

To reproduce the R package environment used for the simulations:

``` r
install.packages("renv")
renv::restore()
```

For more information about the publicly available data that is used for
this simulation study, see:

Liu, J., Lichtenberg, T., Hoadley, K. A., Poisson, L. M., Lazar, A. J.,
Cherniack, A. D., Kovatich, A. J., Benz, C. C., Levine, D. A., Lee, A.
V., Omberg, L., Wolf, D. M., Shriver, C. D., Thorsson, V., Cancer Genome
Atlas Research Network, & Hu, H. (2018). An Integrated TCGA Pan-Cancer
Clinical Data Resource to Drive High-Quality Survival Outcome Analytics.
Cell, 173(2), 400–416.e11.

The data published alongside the above manuscript is in this
repository’s Data folder.

# Reproducing the Simulations

After cloning this repository to your local machine, the main simulation
can be run with code/run_local_sim.R. This runs the simulation for all
parameter combinations assessed in the manuscrip and is set up to be run in parallel using `furrr`.
Since the simulations are computationally intensive, especially for settings with large numbers
of imputations, we recommend reproducing the results for
select parameter combinations by modifying the `params` dataframe at the
beginning of the script. For example, to run 500 simulations associated
with the LIHC setting, a sample size of 100, varying imputation number
(10 and 20), and beta = 0, define params as follows:

``` r
params <- crossing(
  cancer_type = c("LIHC"),
  sample_size = c(100),
  n_imps = c(10, 20),
  n_sims = 500,
  beta = c(0)
)
```

For users wanting to reproduce the results from the entire simulation
study, we strongly recommend using a cluster. `tcga_simulations.slurm` provides
an example batch submission script. The user and account fields, as well
as the cluster directory path should be filled in so the simulations can
run in your cluster environment. It works by defining an array over all
the simulation settings and interfacing with `run_cluster_sim.R` to run
submit jobs associated with each setting’s parameters. As before, a
subset of the simulations can be run by modifying the `params` dataframe
in `run_cluster_sim.R`.

Simulation results will be saved to results/ in subfolders associated
with each cancer type from the TCGA-CDR data on which the simulations
were based. File names are tagged with the structure
`_<cancer_type>_n<sample_size>_M<n_imps>_b<beta>` to identify the
simulation parameters. Files starting with `bias_data` contain the
simulation level estimates, truth, and bias for the state occupation
probability estimates over time. Files starting with `bias_summary`
aggregate these results across simulations and compare thee average bias
for the estimates over time. Files starting with `coverage_data` report
the joint coverage probabilities over time for each of the studies
methods.

For users who do not want to fun the full simulation study, precomputed
simulation results are provided in `results_precomputed`. These results
can be used to reproduce the tables and figures presented in the paper.

# Run Time

Computation time for the simulation settings depends on the number of imputations
used in the proposed algorithms. When using M=10 imputations, generating results for a 
single simulation setting takes approximately 15 minutes. However, settings with 
M=200 imputations takes upwards of 4 hours to run.

# Reproducing the Figures

After running the simulations, or by using the pre-computed output, the
figures from the manuscript can be reproduced using:

The resulting figures will be saved to figures/.


# Citation

If you use the code or simulation results from this repository, please
cite the accompanying paper.

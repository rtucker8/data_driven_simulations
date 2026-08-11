#Run all of the simulation settings in parallel on a local machine-
library(tidyverse)
library(furrr)

source("Code/2_data_driven_simulations.R") #function to run one simulation

params <- crossing(
  cancer_type = c("LIHC"),
  sample_size = c(100),
  n_imps = c(50, 100, 150),
  n_sims = 500,
  beta = c(0, 0.75)
)

plan(multisession, workers = availableCores() - 1)

results <- future_pmap(
  params,
  function(sample_size, n_sims, n_imps, cancer_type, beta) {
    run_sim(
      sample_size = sample_size,
      n_sims = n_sims,
      n_imps = n_imps,
      cancer_type = cancer_type,
      beta = beta
    )
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

plan(sequential)

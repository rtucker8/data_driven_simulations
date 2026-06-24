#Run all of the simulation settings in parallel on a local machine-
library(tidyverse)
library(furrr)

#For each cancer type, vary the number of people, imputations, and logHR for T1
params <- crossing(
  cancer_type = c("LIHC", "MESO", "PRAD", "UCEC", "KIRC", "KICH"),
  sample_size = c(50, 100, 200, 500),
  n_imps = c(10, 20, 30, 50),
  n_sims = 500,
  beta = c(-0.75, 0)
)

plan(multisession, workers = availableCores() - 1)

results <- future_pmap(
  params[which(params$cancer_type %in% c("LIHC", "PRAD", "MESO")),],
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

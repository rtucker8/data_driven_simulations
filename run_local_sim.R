#Run all of the simulation settings in parallel on a local machine-
library(tidyverse)
library(furrr)


params <- crossing(
  cancer_type = c("LIHC", "MESO", "PRAD", "UCEC", "KIRC", "KICH"),
  sample_size = c(50, 100, 200, 500),
  n_imps = c(10, 20, 30, 50),
  n_sims = 500,
  beta = c(-0.75, 0),
  approach = c("trWald", "agresticoull")
)

plan(multisession, workers = availableCores() - 1)

results <- future_pmap(
  params[1:2,],
  function(sample_size, n_sims, n_imps, cancer_type, beta, approach) {
    run_sim(
      sample_size = sample_size,
      n_sims = n_sims,
      n_imps = n_imps,
      cancer_type = cancer_type,
      beta = beta,
      approach = approach
    )
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

plan(sequential)

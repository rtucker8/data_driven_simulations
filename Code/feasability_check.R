#Check how many simulations out of 500 have at least one observed ill to dead
#transition. This is important to ensure the imputation methods are applicable to the
#simulated datasets

  #Load Libraries
  library(tidyverse)
  library(survival)
  library(mstate) #implements Aalen-Johansen estimation
  library(msmi) #proposed methods
  library(ggtern) #plot on the simplex
  library(sp) #check if point is inside of a convex hull
  library(PWEXP) #piecewise exponential distribution

  # Load functions, true parameters, and random seeds ---------------------------------------------------------------
  source("Code/0_data_driven_simulation_helper.R") #helper functions for this simulation study
  load("true_parameters.RData") #true piecewise exponential parameters for each setting
  random.seeds <- read_csv("randomSeeds.csv")$simulationSeeds #random seeds for data generation to ensure reproducibility

  params <- crossing(
    cancer_type = c("LIHC", "MESO", "PRAD", "UCEC", "KIRC", "KICH"),
    sample_size = c(50, 100, 200, 500),
    n_imps = c(2),
    n_sims = 500,
    beta = c(-0.75, 0),
    approach = c("trWald"))

  test <- function(cancer_type, n_imps, n_sims, beta, approach, sample_size) {

    #Get random seeds
    seeds = random.seeds[1:n_sims]

    #Simulate datasets with reproducible random seeds using helper function
    d.sim <- map(seeds, function(s) {
      simulate_illness_death(
        n = sample_size,
        params = truth[[cancer_type]],
        seed = s,
        beta = beta,
        return_latent_data = FALSE,
        return_censored_data = TRUE
      )
    })

    warning_indices <- integer(0) #simulations with no observed ill to dead transitions so imputation is not complete

    imps_marginal <- map2(seq_along(d.sim), d.sim, function(i, d) {

      withCallingHandlers(
        msmi.impute(
          dat = d,
          M = n_imps,
          method = "marginal",
          seed = seeds[i]
        ),
        warning = function(w) {
          if (grepl("There are no uncensored transitions between the first and second event in this dataset.
            The second layer of imputation cannot be performed.", conditionMessage(w))) {
            warning_indices <<- c(warning_indices, i)
          }
          invokeRestart("muffleWarning")
        }
      )

    })

    #Only compute estimates for imputations that could be completed
    good_settings <- setdiff(1:n_sims, warning_indices)

    return(length(good_settings))
  }

  results <- pmap(
    params,
    function(sample_size, n_sims, n_imps, cancer_type, beta, approach) {
      test(
        sample_size = sample_size,
        n_sims = n_sims,
        n_imps = n_imps,
        cancer_type = cancer_type,
        beta = beta,
        approach = approach
      )
    })


params$good_simulations <- unlist(results, use.names=FALSE)

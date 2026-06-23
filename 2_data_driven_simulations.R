#Purpose: Simulation Study to accompany Gonzalez, Dempsey, and Boonstra (2027+)
#Author: Rachel Gonzalez
#Date: Sep 11 2025

#Input: true_occupation_probabilities.RData, a list containing the true state occupation probabilities for each simulation setting
#Input: true_parameters.RData, a list containing the true piecewise exponential parameters for each transition and censoring distribution in each simulation setting
#Input: randomSeeds.csv, a csv file containing 500 random seeds to use for reproducibility

#Load Libraries
library(tidyverse)
library(survival)
library(mstate) #implements Aalen-Johansen estimation
library(msmi) #proposed methods
library(ggtern) #plot on the simplex
library(sp) #check if point is inside of a convex hull
library(furrr) #parallel processing
library(PWEXP)

# Load functions, true parameters, and random seeds ---------------------------------------------------------------
source("0_data_driven_simulation_helper.R") #helper functions for this simulation study
load("true_occupation_probabilities.RData") #true state occupation probabilities for each setting
load("true_parameters.RData") #true piecewise exponential parameters for each setting
random.seeds <- read_csv("randomSeeds.csv")$simulationSeeds #random seeds for data generation to ensure reproducibility

# Simulation Setting ------------------------------------------------------

# Some Sample Parameters for Debugging
sample_size = 100
n_sims = 500
n_imps = 10
cancer_type = "KIRC"
beta = 0

#Main Function: runs a simulation for a given setting
  #sample_size: number of patients in each simulated dataset
  #n_sims: number of simulated datasets to generate
  #n_imps: number of imputations to use for the MSMI methods
  #cancer_type: one of "LIHC", "PRAD", "MESO", "KICH", "KIRC", "UCEC"
  #beta: logHR associated with time of arrival to the illness state for the ill to dead transition hazard

run_sim <- function(sample_size, n_sims, n_imps, cancer_type, beta) {

  ########################
  ##   Preliminaries   ##
  ########################

  #Evaluate methods at 1,2,3,4 and 5 years
  eval_times <- seq(365, 365 * 5, 365)

  #lookup true state occupation probabilities for this setting
  truth_setting <- empirical_truth[[cancer_type]][[paste0("beta",as.character(beta))]] %>%
    pivot_wider(names_from = state, values_from = probability)
  truth_setting <- truth_setting %>% filter(time %in% eval_times)

  #Diagnostic Check (not saved): plot how truth changes over time
  suppressWarnings({
    p <- ggtern::ggtern(
      truth_setting,
      ggplot2::aes(x = Healthy, y = Ill, z = Dead, color = factor(time))
    ) +
      ggplot2::geom_point(size = 2) +
      ggtern::theme_bw() +
      ggtern::theme_showarrows() +
      ggplot2::xlab("Healthy") +
      ggplot2::ylab("Ill") +
      ggtern::zlab("Dead") +
      ggplot2::ggtitle("True State Occupation Probabilities over Time") +
      ggplot2::labs(color = "Time (days)")
  })
  p

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

  #Diagnostic Check (not saved): look at the event rates
  check_events <- function(df) {
    df.long <- msprep(
      data = df,
      trans = transMat(
        x = list(c(2, 3), c(3), c()),
        names = c("Healthy", "Ill", "Death")
      ),
      time = c(NA, "t1", "t2"),
      status = c(NA, "event1", "event2")
    )
    events(df.long)$Proportions
  }

check_events(d.sim[[1]])

#file identifier
tag <- paste0(cancer_type,
              "_",
              as.character(sample_size),
              "_size_",
              as.character(n_imps),
              "_imps_",
              as.character(beta),
              "beta")

  ################################
  ##   Method: Aalen-Johansen   ##
  ################################

  #create transition matrix for illness death model
  tmat <- transMat(
    x = list(c(2, 3), c(3), c()),
    names = c("Healthy", "Ill", "Death")
  )

  #function to compute AJ state occupation probabilities for a given dataset
  aj_estimation <- function(df) {
    #prep data
    df.long <- msprep(
      data = df,
      trans = tmat,
      time = c(NA, "t1", "t2"),
      status = c(NA, "event1", "event2")
    )

    #Markov model without covariates- fully nonparametric
    c0 <- coxph(
      Surv(Tstart, Tstop, status) ~ strata(trans),
      data = df.long,
      method = "breslow"
    )
    msf0 <- msfit(object = c0, vartype = "aalen", trans = tmat)

    #state occupation probabilities and confidence region
    pt <- probtrans(msf0, predt = 0, covariance = TRUE)
    pt0 <- pt[[1]]
    aj.cov <- pt[["varMatrix"]][
      c("from1to1", "from1to2", "from1to3"),
      c("from1to1", "from1to2", "from1to3"),
    ]
    return(list(pt0 = pt0, cov = aj.cov))
  }

  aj <- map(d.sim, aj_estimation)

  #AJ estimates for each state at eval_times
  aj_probabilities <- map(aj, function(df) {
    tibble(
      time = eval_times,
      pstate1 = map_dbl(
        eval_times,
        ~ get_step_value(df[['pt0']]$time, df[["pt0"]]$pstate1, .x)
      ),
      pstate2 = map_dbl(
        eval_times,
        ~ get_step_value(df[['pt0']]$time, df[['pt0']]$pstate2, .x)
      ),
      pstate3 = map_dbl(
        eval_times,
        ~ get_step_value(df[['pt0']]$time, df[['pt0']]$pstate3, .x)
      )
    )
  })

  aj_results <- bind_rows(aj_probabilities, .id = "simulation") %>%
    select(simulation, time, pstate1, pstate2, pstate3)

  #Diagnostic Check (not saved): plot AJ estimates on top of the truth
  suppressWarnings({
    p_aj <- ggtern::ggtern(
      truth_setting,
      ggplot2::aes(x = Healthy, y = Ill, z = Dead)
    ) +
      ggplot2::geom_point(
        data = aj_results,
        aes(x = pstate1, y = pstate2, z = pstate3, color = factor(time)),
        size = 1
      ) +
      ggplot2::geom_point(size = 2, shape = 18) +
      ggtern::theme_bw() +
      ggtern::theme_showarrows() +
      ggplot2::xlab("Healthy") +
      ggplot2::ylab("Ill") +
      ggtern::zlab("Dead") +
      ggplot2::ggtitle(
        "State Occupation Probabilities over Time for AJ Method"
      ) +
      ggplot2::labs(color = "Time (days)")
  })
  p_aj

  #Get AJ confidence regions for state occupation probability vector at eval_times
  RegionAJ <- map2(aj, aj_probabilities, function(df, probs) {
    map(eval_times, function(t) {
      idx <- max(which(df[['pt0']]$time <= t), na.rm = TRUE)

      AJ_Region(
        probs %>%
          filter(time == t) %>%
          dplyr::select(pstate1, pstate2, pstate3) %>%
          as.matrix() %>%
          as.numeric(),

        df[['cov']][,, idx]
      )
    }) %>%
      setNames(as.character(eval_times))
  })


  #################################################################
  ##   Diagnostic Check: Behavior when no imputation is needed   ##
  #################################################################

  # #Calculate Point Estimates
  # p_est <- map(d.sim, function(sim) {
  #   pt_est <- empirical_transition_probs(sim, s = 0, from = 1, times = seq(365, 365*5, 365)) %>% select(-stime)
  # }) %>% bind_rows(.id = "simulation") %>%
  #   group_by(simulation) %>%
  #   group_split()
  #
  # #Compute variance
  # variance <- function(p, n = 100) {
  #   (diag(p) - tcrossprod(p)) / n
  # }
  #
  # #Get regions
  # RegionTest <- map(p_est, function(probs) {
  #   map(eval_times, function(t) {
  #    est = probs %>%
  #      filter(time == t) %>%
  #      dplyr::select(pstate1, pstate2, pstate3) %>%
  #      as.matrix() %>%
  #      as.numeric()
  #
  #     AJ_Region(
  #       est,
  #       variance(est)
  #     )
  #   }) %>%
  #     setNames(as.character(eval_times))
  # })

  ################################
  ##   Method: MSMI- Marginal   ##
  ################################

  #Create multiple imputed datasets
  imps_marginal <- map2(d.sim, seeds, function(d, s) {
    msmi.impute(
      dat = d,
      M = n_imps,
      type = "km",
      method = "marginal",
      concentration = 1,
      seed = s
    )
  })

  #Estimated state occupation probabilities from the imputations
  marginal_tprobs <- map(imps_marginal, function(imp) {
    msmi.tprobs(
      imp_obj = imp,
      times = eval_times,
      #int.type = "agresticoull",
      int.type = "trWald",
      alpha = 0.05
    )
  })

  mi_estimate <- map(marginal_tprobs, function(obj) {
    obj$mi_estimate
  })

  mi_results <- bind_rows(mi_estimate, .id = "simulation") %>%
    select(simulation, time, p1, p2, p3) %>%
    rename(pHealthy = p1, pIll = p2, pDead = p3)

  #Diagnostic Check (Not Saved): Plot MSMI point estimates compared to the truth
  suppressWarnings({
    p_msmi <- ggtern::ggtern(
      truth_setting,
      ggplot2::aes(x = Healthy, y = Ill, z = Dead)
    ) +
      ggplot2::geom_point(
        data = mi_results,
        aes(x = pHealthy, y = pIll, z = pDead, color = factor(time)),
        size = 1
      ) +
      ggplot2::geom_point(size = 2, shape = 18) +
      ggtern::theme_bw() +
      ggtern::theme_showarrows() +
      ggplot2::xlab("Healthy") +
      ggplot2::ylab("Ill") +
      ggtern::zlab("Dead") +
      ggplot2::ggtitle(
        "State Occupation Probabilities over Time for MSMI-Marginal Method"
      ) +
      ggplot2::labs(color = "Time (days)")
  })
  p_msmi

  #Diagnostic Check (Not Saved): Visualize correlation of MSMI point estimates and AJ point estimates
  mi_results |>
    rename(pstate1 = pHealthy, pstate2 = pIll, pstate3 = pDead) |>
    left_join(
      aj_results,
      by = c("simulation", "time"),
      suffix = c("_msmi", "_aj")
    ) |>
    pivot_longer(
      cols = starts_with("p"),
      names_to = c("state", "method"),
      names_pattern = "p(.*)_(.*)",
      values_to = "estimate"
    ) |>
    pivot_wider(
      names_from = method,
      values_from = estimate
    ) |>
    ggplot(aes(x = aj, y = msmi)) +
    geom_point(alpha = 0.6) +
    facet_grid(
      state ~ time,
      labeller = labeller(
        time = c(
          "365" = "365d",
          "730" = "730d",
          "1095" = "1095d",
          "1460" = "1460d",
          "1825" = "1825d"
        ),
        state = c("Healthy" = "Healthy", "Ill" = "Ill", "Dead" = "Dead")
      )
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "red",
      alpha = 0.5
    ) +
    labs(
      title = "MSMI vs AJ State Occupation Probability Estimates",
      x = "AJ Estimate",
      y = "MSMI Estimate"
    ) +
    coord_equal()

  ################################
  ##   Method: MSMI-Cox        ##
  ################################

  imps_cox <- map2(d.sim, seeds, function(d, s) {
    msmi.impute.debug(dat = d, M = n_imps, method = "cox", seed = s)
  })

  #Estimated state occupation probabilities (Rubin's Rules Estimate)
  cox_tprobs <- map(imps_cox, function(imp) {
    msmi.tprobs(
      imp_obj = imp,
      times = eval_times,
      int.type = "trWald",
      alpha = 0.05
    )
  })

  cox_mi_estimate <- map(cox_tprobs, function(obj) {
    obj$mi_estimate
  })

  cox_mi_results <- bind_rows(cox_mi_estimate, .id = "simulation") %>%
    select(simulation, time, p1, p2, p3) %>%
    rename(pHealthy = p1, pIll = p2, pDead = p3)


  ################################
  ##   Compare Methods         ##
  ################################

  #Create and save bias dataset
  bias <- bind_rows(

    #Aalen Johansen
    aj_results %>%
      mutate(method = "Aalen-Johansen"),

    #MSMI-KM
    mi_results %>%
      rename(
        pstate1 = pHealthy,
        pstate2 = pIll,
        pstate3 = pDead
      ) %>%
      mutate(method = "Marginal MI"),

    #MSMI-Cox
    cox_mi_results %>%
      rename(
        pstate1 = pHealthy,
        pstate2 = pIll,
        pstate3 = pDead
      ) %>%
      mutate(method = "Cox MI")
  ) %>%
    left_join(truth_setting, by = "time") %>%
    pivot_longer(
      cols = starts_with("pstate"),
      names_to = "state",
      values_to = "estimate"
    ) %>%
    mutate(
      truth = case_when(
        state == "pstate1" ~ Healthy,
        state == "pstate2" ~ Ill,
        state == "pstate3" ~ Dead
      ),
      bias = estimate - truth,
      state = sub("pstate", "", state),
      time = factor(round(time, 2))
    ) %>%
    select(time, simulation, method, state, bias, estimate, truth) %>%
    filter(time != 0)


  write_csv(
    bias,
    paste0(
      "Output/bias_data_",
      tag,
      ".pdf"
    )
  )

  #Bias Figure
  ggplot(bias) +
    geom_boxplot(
      aes(
        x = factor(as.numeric(as.character(time)) / 365),
        y = bias,
        fill = method
      ),
      outliers = F
    ) +
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
    facet_wrap(~state, nrow = 3, ncol = 1) +
    labs(
      title = "Bias of State Occupation Probability Estimates",
      x = "Time (Years)",
      y = "Bias"
    ) +
    scale_fill_brewer(type = "qual", palette = 6)

  ggsave(
    paste0(
      "Output/Figures/bias_",
      tag,
      "pdf"
    ),
    height = 10,
    width = 8
  )

  #Coverage

  #Build a dataframe that contains the points defining the covex_hulls that represent confidence regions for each estimate
  is_missing <- function(x) is.null(x) || (length(x) > 0 && all(is.na(x)))

  # truth lookup (one row per time)
  truth_lookup <- truth_setting |>
    dplyr::select(time, Healthy, Ill, Dead)

  # helper: does hull contain truth at time t?
  contains_truth_hull <- function(hull_xy, t, truth_tbl) {

    tr <- truth_tbl |>
      dplyr::filter(time == t)

    if (is_missing(hull_xy) || nrow(hull_xy) < 3) {
      return(FALSE)
    }

    sp::point.in.polygon(
      point.x = tr$Healthy,
      point.y = tr$Ill,
      pol.x = hull_xy[, 1],
      pol.y = hull_xy[, 2]
    ) >
      0
  }

  #helper: makes the hull dataframe for a given method
  make_hull_df <- function(hull_obj, simulation, time, method) {

    if (is_missing(hull_obj))
      return(NULL)

    hull <- as.matrix(hull_obj)

    if (is_missing(hull) || nrow(hull) < 3)
      return(NULL)

    time_copy = time

    tibble(
      simulation = as.character(simulation),
      time = time,
      method = method,
      x = hull[, 1],
      y = hull[, 2],
      z = if (ncol(hull) >= 3) {
        hull[, 3]
      } else {
        1 - hull[, 1] - hull[, 2]
      },
      contains_truth = contains_truth_hull(
        hull[, 1:2, drop = FALSE],
        time_copy,
        truth_lookup
      )
    )
  }

  #main dataset that contains the confidence regions for each method at each eval_time
  polygon_df <-
    map_dfr(seq_along(marginal_tprobs), \(i) {

      map_dfr(eval_times, \(t) {

        bind_rows(
          make_hull_df(
            marginal_tprobs[[i]][["cr_list"]][[as.character(t)]][["p.space"]],
            i, t, "MSMIKM"
          ),
          make_hull_df(
            cox_tprobs[[i]][["cr_list"]][[as.character(t)]][["p.space"]],
            i, t, "MSMICox"
          ),
          make_hull_df(
            RegionAJ[[i]][[as.character(t)]],
            i, t, "AJ"
          )
        )

      })

    }) %>%
    mutate(
      poly_id = paste(simulation, time, method, sep = "_")
    )

  #Calculate coverage with same denominator across methods
  #Exclude timepoint/simulation combinations where both confidence regions are not defined
  coverage = polygon_df %>%
    group_by(simulation, time, method) %>%
    slice_head(n=1) %>%
    ungroup() %>%
    select(-poly_id) %>%
    pivot_wider(names_from = "method", values_from = c(x,y,z, contains_truth) ) %>%
    group_by(time) %>%
    mutate(both_defined = contains_truth_AJ + contains_truth_MSMIKM + contains_truth_MSMICox) %>%
    filter(!is.na(both_defined)) %>%
    summarise(coverageAJ = mean(contains_truth_AJ),
              coverageMSMI = mean(contains_truth_MSMIKM),
              coverageMSMICox = mean(contains_truth_MSMICox))

  write_csv(
    coverage,
    paste0(
      "Output/coverage_data_",
      tag,
      "pdf"
    )
  )

  #Diagnostic Check: plot the point and interval estimates for each simulation and under each method

  plot_sim_regions <- function(sim_id, facet_cols = 3) {
    sim_str <- as.character(sim_id)
    print(paste0("Progress: Simulation ", sim_id))
    hulls <- polygon_df %>% filter(simulation == sim_str)
    hulls_fill <- hulls |>
      dplyr::filter(contains_truth %in% TRUE)
    hulls_outline <- hulls |>
      dplyr::filter(!(contains_truth %in% TRUE))
    pts <- truth_setting
    ms_pts <- mi_results %>% filter(simulation == sim_str)
    aj_pts <- aj_results %>% filter(simulation == sim_str)

    suppressWarnings({
      ggtern::ggtern() +
        # filled if contains truth
        ggplot2::geom_polygon(
          data = hulls_fill,
          ggplot2::aes(x = x, y = y, z = z, group = poly_id, fill = method),
          alpha = 0.35,
          color = NA
        ) +
        # outline only if does not contain truth
        ggplot2::geom_polygon(
          data = hulls_outline,
          ggplot2::aes(x = x, y = y, z = z, group = poly_id, color = method),
          fill = NA,
          linewidth = 0.5
        ) +
        ggplot2::geom_point(
          data = pts,
          ggplot2::aes(x = Healthy, y = Ill, z = Dead),
          color = "black",
          shape = 18,
          size = 2
        ) +
        ggplot2::geom_point(
          data = ms_pts,
          ggplot2::aes(x = pHealthy, y = pIll, z = pDead, color = "MSMI"),
          size = 1.4
        ) +
        ggplot2::geom_point(
          data = aj_pts,
          ggplot2::aes(x = pstate1, y = pstate2, z = pstate3, color = "AJ"),
          size = 1.4
        ) +
        ggplot2::facet_wrap(~ factor(time), ncol = facet_cols) +
        ggplot2::scale_fill_manual(
          values = c(MSMI = "#1b9e77", AJ = "#d95f02")
        ) +
        ggplot2::scale_color_manual(
          values = c(MSMI = "#1b9e77", AJ = "#d95f02")
        ) +
        ggtern::theme_bw() +
        ggtern::theme_showarrows() +
        ggplot2::labs(
          title = paste("Confidence regions - simulation", sim_str),
          x = "Healthy",
          y = "Ill",
          z = "Dead",
          fill = "Filled region",
          color = "Method"
        )
    })
  }

  # #Save all figures to one PDF
  # out_file <- paste0("confidence_regions_agresti_coull", cancer_type, ".pdf")
  # sim_ids <- seq_along(marginal_tprobs)
  #
  # grDevices::pdf(out_file, width = 12, height = 8, onefile = TRUE)
  # purrr::walk(sim_ids, \(i) {
  #   p <- plot_sim_regions(i)
  #   plot(p)
  # })
  # grDevices::dev.off()


  return()
}

#Errors caused by msmi_cox method in new settings (KIRC for one):
#basically, zero survival probability past a certain point ($surv) due to large times to illness fed into cox model
# Error in `map2()`:
#   ℹ In index: 335.
# Caused by error in `purrr::map()`:
#   ℹ In index: 1.
# Caused by error in `sample.int()`:
#   ! too few positive probabilities
# Run `rlang::last_trace()` to see where the error occurred.
# Called from: signal_abort(cnd, .file)

#No error with larger sample sizes (200 vs 100 for HR = 0) or protective HR (-0.75)

#Run all of the simulation settings in parallel-

# #For each cancer type, vary the number of people, imputations, and logHR for T1
# params <- crossing(
#   cancer_type = c("LIHC", "MESO", "PRAD", "UCEC", "KIRC", "KICH"),
#   sample_size = c(50, 100, 200, 500),
#   n_imps = c(10, 20, 30, 50),
#   n_sims = 500,
#   beta = c(-0.75, 0)
# )
#
# plan(multisession, workers = availableCores() - 1)
#
# results <- future_pmap(
#   params,
#   function(sample_size, n_sims, n_imps, cancer_type, beta) {
#     run_sim(
#       sample_size = sample_size,
#       n_sims = n_sims,
#       n_imps = n_imps,
#       cancer_type = cancer_type,
#       beta = beta
#     )
#   },
#   .options = furrr_options(seed = TRUE),
#   .progress = TRUE
# )
#
# plan(sequential)

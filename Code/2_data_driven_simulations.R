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
library(PWEXP) #piecewise exponential distribution

# Load functions, true parameters, and random seeds ---------------------------------------------------------------
source("R/simulation_helper.R") #helper functions for this simulation study
load("Results/true_occupation_probabilities.RData") #true state occupation probabilities for each setting
load("Results/true_parameters.RData") #true piecewise exponential parameters for each setting
random.seeds <- read_csv("Data/randomSeeds.csv")$simulationSeeds #random seeds for data generation to ensure reproducibility

# Simulation Setting ------------------------------------------------------

#run_sim: runs a simulation for a given setting
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
  truth_setting <- empirical_truth[[cancer_type]][[paste0(
    "beta",
    as.character(beta)
  )]] %>%
    pivot_wider(names_from = state, values_from = probability)
  truth_setting <- truth_setting %>% filter(time %in% eval_times)

  #Get random seeds
  seeds <- random.seeds[1:n_sims]

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

  #Diagnositc Check (not saved): examine event rates for each transition
  event_rates <- map(d.sim, function(df) {
    df.long <- msprep(
      data = df,
      trans = transMat(
        x = list(c(2, 3), c(3), c()),
        names = c("Healthy", "Ill", "Death")
      ),
      time = c(NA, "t1", "t2"),
      status = c(NA, "event1", "event2")
    )
    events(df.long)$Proportions %>%
      as.data.frame() %>%
      pivot_wider(names_from = to, values_from = Freq) %>%
      filter(from != "Death") %>%
      select(-Healthy)
  }) %>%
    bind_rows(.id = "simulation")

  event_rates %>%
    group_by(from) %>%
    summarise(Ill = mean(Ill), Death = mean(Death), noevent = mean(`no event`))

  sum(is.na(event_rates[event_rates$from == 'Ill', ]$`no event`)) #num of sims where no one became ill
  sum(event_rates[event_rates$from == 'Ill', ]$`no event` == 1, na.rm = T) #num of sims where all ill to dead transitions are censored

  #file identifier
  tag <- paste0(
    cancer_type,
    "_n",
    as.character(sample_size),
    "_M",
    as.character(n_imps),
    "_b",
    as.character(beta)
  )

  ################################
  ##   Method: MSMI- Marginal   ##
  ################################

  run_imputation <- function(method, bootstrap) {
    map2(
      d.sim,
      seeds,
      ~ msmi.impute(
        dat = .x,
        M = n_imps,
        method = method,
        seed = .y,
        bootstrap = bootstrap
      )
    )
  }

  get_tprobs <- function(imps, approach) {
    map(
      imps,
      ~ msmi.tprobs.v2(
        imp_obj = .x,
        times = eval_times,
        method = approach,
        alpha = 0.05
      )
    )
  }

  get_results <- function(tprobs, bootstrap) {
    map(tprobs, "mi_estimate") %>%
      bind_rows(.id = "simulation") %>%
      select(simulation, time, p1, p2, p3) %>%
      rename(
        pHealthy = p1,
        pIll = p2,
        pDead = p3
      ) %>%
      mutate(bootstrap = bootstrap)
  }

  # No bootstrap
  imps_marginal <- run_imputation("marginal", FALSE)

  marginal_tprobs_goodman <- get_tprobs(imps_marginal, "goodman")
  marginal_tprobs_trWald <- get_tprobs(imps_marginal, "trWald")

  marginal_mi_results <- get_results(marginal_tprobs_goodman, FALSE)

  # Bootstrap
  imps_marginal_boot <- run_imputation("marginal", TRUE)

  marginal_tprobs_goodman_boot <- get_tprobs(imps_marginal_boot, "goodman")
  marginal_tprobs_trWald_boot <- get_tprobs(imps_marginal_boot, "trWald")

  marginal_mi_results_boot <- get_results(marginal_tprobs_goodman_boot, TRUE)

  marginal_mi_results <- bind_rows(
    marginal_mi_results,
    marginal_mi_results_boot
  )

  ################################
  ##   Method: MSMI-Cox        ##
  ################################

  # No bootstrap
  imps_cox <- run_imputation("cox", FALSE)

  cox_tprobs_goodman <- get_tprobs(imps_cox, "goodman")
  cox_tprobs_trWald <- get_tprobs(imps_cox, "trWald")

  cox_mi_results <- get_results(cox_tprobs_goodman, FALSE)

  # Bootstrap
  imps_cox_boot <- run_imputation("cox", TRUE)

  cox_tprobs_goodman_boot <- get_tprobs(imps_cox_boot, "goodman")
  cox_tprobs_trWald_boot <- get_tprobs(imps_cox_boot, "trWald")

  cox_mi_results_boot <- get_results(cox_tprobs_goodman_boot, TRUE)

  cox_mi_results <- bind_rows(cox_mi_results, cox_mi_results_boot)

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

  ################################
  ##   Compare Methods         ##
  ################################

  #Create and save bias dataset
  bias <- bind_rows(
    #Aalen Johansen
    aj_results %>%
      mutate(method = "AJ"),

    #MSMI-KM
    marginal_mi_results %>%
      rename(
        pstate1 = pHealthy,
        pstate2 = pIll,
        pstate3 = pDead
      ) %>%
      mutate(method = if_else(bootstrap == FALSE, "MSMI-KM", "MSMI-KM-Boot")),

    #MSMI-Cox
    cox_mi_results %>%
      rename(
        pstate1 = pHealthy,
        pstate2 = pIll,
        pstate3 = pDead
      ) %>%
      mutate(method = if_else(bootstrap == FALSE, "MSMI-Cox", "MSMI-Cox-Boot"))
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
      "Results/",
      tolower(cancer_type),
      "/bias_data_",
      tag,
      ".csv"
    )
  )

  bias_summary <- bias %>%
    group_by(time, method, state) %>%
    summarise(
      time = first(time),
      method = first(method),
      state = first(state),
      avg_bias = mean(bias)
    ) %>%
    pivot_wider(names_from = state, values_from = avg_bias)

  write_csv(
    bias_summary,
    paste0(
      "Results/",
      tolower(cancer_type),
      "/bias_summary_",
      tag,
      ".csv"
    )
  )

  #Bias Figure
  ggplot(
    bias,
    aes(
      x = factor(as.numeric(as.character(time)) / 365),
      y = bias,
      fill = method
    )
  ) +
    geom_boxplot(outliers = FALSE) +
    #geom_violin(scale = "width" ) +
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
    facet_wrap(~state, nrow = 3, ncol = 1) +
    labs(
      title = tag,
      x = "Time (Years)",
      y = "Bias (estimate - truth)"
    ) +
    scale_fill_brewer(type = "qual", palette = 6)

  ggsave(
    paste0(
      "Results/",
      tolower(cancer_type),
      "/Figures/bias_",
      tag,
      ".pdf"
    ),
    height = 10,
    width = 8
  )

  #Coverage

  #Build a dataframe that contains the points defining the covex_hulls that represent confidence regions for each estimate
  is_missing <- function(x) is.null(x) || (length(x) > 0 && all(is.na(x)))

  # truth lookup (one row per time)
  truth_lookup <- truth_setting %>%
    dplyr::select(time, Healthy, Ill, Dead)

  # helper: does hull contain truth at time t?
  contains_truth_hull <- function(hull_xy, t, truth_tbl) {
    tr <- truth_tbl %>%
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
  make_hull_df <- function(hull_obj, simulation, time, method, approach) {
    if (is_missing(hull_obj)) {
      return(NULL)
    }

    hull <- as.matrix(hull_obj)

    if (is_missing(hull) || nrow(hull) < 3) {
      return(NULL)
    }

    time_copy <- time

    tibble(
      simulation = simulation,
      time = time,
      method = method,
      approach = approach,
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
    map_dfr(1:n_sims, \(i) {
      sim <- i

      map_dfr(eval_times, \(t) {
        bind_rows(
          make_hull_df(
            marginal_tprobs_goodman[[i]][["cr_list"]][[as.character(t)]],
            sim,
            t,
            "MSMI-KM",
            "goodman"
          ),
          make_hull_df(
            marginal_tprobs_trWald[[i]][["cr_list"]][[as.character(t)]],
            sim,
            t,
            "MSMI-KM",
            "trWald"
          ),
          make_hull_df(
            marginal_tprobs_goodman_boot[[i]][["cr_list"]][[as.character(t)]],
            sim,
            t,
            "MSMI-KM-Boot",
            "goodman"
          ),
          make_hull_df(
            marginal_tprobs_trWald_boot[[i]][["cr_list"]][[as.character(t)]],
            sim,
            t,
            "MSMI-KM-Boot",
            "trWald"
          ),
          make_hull_df(
            cox_tprobs_goodman[[i]][["cr_list"]][[as.character(t)]],
            sim,
            t,
            "MSMI-Cox",
            "goodman"
          ),
          make_hull_df(
            cox_tprobs_trWald[[i]][["cr_list"]][[as.character(t)]],
            sim,
            t,
            "MSMI-Cox",
            "trWald"
          ),
          make_hull_df(
            cox_tprobs_goodman_boot[[i]][["cr_list"]][[as.character(t)]],
            sim,
            t,
            "MSMI-Cox-Boot",
            "goodman"
          ),
          make_hull_df(
            cox_tprobs_trWald_boot[[i]][["cr_list"]][[as.character(t)]],
            sim,
            t,
            "MSMI-Cox-Boot",
            "trWald"
          ),
          make_hull_df(
            RegionAJ[[i]][[as.character(t)]],
            sim,
            t,
            "AJ",
            "NA"
          )
        )
      })
    }) %>%
    mutate(
      poly_id = paste(simulation, time, method, approach, sep = "_")
    )

  #Calculate coverage with same denominator across methods
  #Exclude simulations where there were no observed ill to dead transitions
  #Exclude timepoint/simulation combinations where both confidence regions are not defined

  # One row per simulation/time/method/approach

  coverage_long <- polygon_df %>%
    group_by(simulation, time, method, approach) %>%
    summarise(contains_truth = first(contains_truth), .groups = "drop") %>%
    pivot_wider(
      names_from = c(method, approach),
      values_from = contains_truth
    ) %>%
    mutate(across(matches("MSMI|AJ"), ~ replace_na(.x, FALSE)))

  coverage <- coverage_long %>%
    group_by(time) %>%
    summarise(
      AJ = mean(AJ_NA),
      `MSMI-KM-Goodman` = mean(`MSMI-KM_goodman`),
      `MSMI-KM-trWald` = mean(`MSMI-KM_trWald`),
      `MSMI-Cox-Goodman` = mean(`MSMI-Cox_goodman`),
      `MSMI-Cox-trWald` = mean(`MSMI-Cox_trWald`),
      `MSMI-KM-Boot-Goodman` = mean(`MSMI-KM-Boot_goodman`),
      `MSMI-KM-Boot-trWald` = mean(`MSMI-KM-Boot_trWald`),
      `MSMI-Cox-Boot-Goodman` = mean(`MSMI-Cox-Boot_goodman`),
      `MSMI-Cox-Boot-trWald` = mean(`MSMI-Cox-Boot_trWald`)
    )

  write_csv(
    coverage,
    paste0(
      "Results/",
      tolower(cancer_type),
      "/coverage_data_",
      tag,
      ".csv"
    )
  )

  #Misc. Checks and Visualizations

  # #Diagnostic Check (not saved): plot estimates on top of the truth
  # suppressWarnings({
  #   p_estimates <- ggtern::ggtern(
  #     truth_setting,
  #     ggplot2::aes(x = Healthy, y = Ill, z = Dead)
  #   ) +
  #     ggplot2::geom_point(
  #       data = bias %>%
  #         select(simulation, time, method, state, estimate) %>%
  #         pivot_wider(
  #           names_from = state,
  #           values_from = estimate,
  #           names_prefix = "pstate"
  #         ),
  #       aes(x = pstate1, y = pstate2, z = pstate3, color = factor(time)),
  #       size = 1,
  #       alpha = 0.5
  #     ) +
  #     facet_wrap(~ as.factor(method)) +
  #     ggplot2::geom_point(size = 2, shape = 18) +
  #     ggtern::theme_bw() +
  #     ggtern::theme_showarrows() +
  #     ggplot2::xlab("Healthy") +
  #     ggplot2::ylab("Ill") +
  #     ggtern::zlab("Dead") +
  #     ggplot2::ggtitle(
  #       "State Occupation Probability Estimates over Time"
  #     ) +
  #     ggplot2::labs(color = "Time (days)")
  # })
  # p_estimates

  # #Diagnostic Check (Not Saved): Visualize correlation of MSMI point estimates and AJ point estimates
  # bias %>%
  #   select(simulation, time, state, method, estimate) %>%
  #   filter(method %in% c("AJ", "MSMI-KM")) %>%
  #   mutate(
  #     state = case_when(
  #       state == 1 ~ "Healthy",
  #       state == 2 ~ "Ill",
  #       state == 3 ~ "Dead"
  #     ) %>%
  #       factor(levels = c("Healthy", "Ill", "Dead"), ordered = T)
  #   ) %>%
  #   pivot_wider(names_from = method, values_from = estimate) %>%
  #   ggplot(aes(x = AJ, y = `MSMI-KM`)) +
  #   geom_point(alpha = 0.2) +
  #   facet_grid(
  #     state ~ time,
  #     labeller = labeller(
  #       time = c(
  #         "365" = "365d",
  #         "730" = "730d",
  #         "1095" = "1095d",
  #         "1460" = "1460d",
  #         "1825" = "1825d"
  #       ),
  #       state = c("Healthy" = "Healthy", "Ill" = "Ill", "Dead" = "Dead")
  #     )
  #   ) +
  #   geom_abline(
  #     intercept = 0,
  #     slope = 1,
  #     linetype = "dashed",
  #     color = "red",
  #     alpha = 0.5
  #   ) +
  #   labs(
  #     title = "MSMI-KM vs AJ State Occupation Probability Estimates",
  #     x = "AJ Estimate",
  #     y = "MSMI-KM Estimate"
  #   ) +
  #   coord_equal()

  # #Diagnostic Check: plot the point and interval estimates for each simulation and under each method
  # plot_sim_regions <- function(sim_id, time_num) {
  #   sim_str <- as.character(sim_id)
  #   print(paste0("Progress: Simulation ", sim_id))
  #   hulls <- polygon_df %>%
  #     filter(
  #       method %in% c("AJ", "MSMIKM"),
  #       simulation == sim_str,
  #       time == time_num
  #     )

  #   hulls_fill <- hulls %>%
  #     dplyr::filter(contains_truth %in% TRUE)
  #   hulls_outline <- hulls %>%
  #     dplyr::filter(!(contains_truth %in% TRUE))
  #   pts <- truth_setting %>% filter(time == time_num)

  #   est_pts <- bind_rows(
  #     mi_results %>%
  #       filter(simulation == sim_str, time == time_num) %>%
  #       mutate(method = "MSMIKM"),
  #     aj_results %>%
  #       rename(pHealthy = pstate1, pIll = pstate2, pDead = pstate3) %>%
  #       filter(simulation == sim_str, time == time_num) %>%
  #       mutate(method = "AJ")
  #   )

  #   suppressWarnings({
  #     ggtern::ggtern() +
  #       # filled if contains truth
  #       ggplot2::geom_polygon(
  #         data = hulls_fill,
  #         ggplot2::aes(x = x, y = y, z = z, group = poly_id, fill = method),
  #         alpha = 0.35,
  #         color = NA
  #       ) +
  #       # outline only if does not contain truth
  #       ggplot2::geom_polygon(
  #         data = hulls_outline,
  #         ggplot2::aes(x = x, y = y, z = z, group = poly_id, color = method),
  #         fill = NA,
  #         linewidth = 0.5
  #       ) +
  #       ggplot2::geom_point(
  #         data = pts,
  #         ggplot2::aes(x = Healthy, y = Ill, z = Dead),
  #         color = "black",
  #         shape = 18,
  #         size = 2
  #       ) +
  #       ggplot2::geom_point(
  #         data = est_pts,
  #         ggplot2::aes(x = pHealthy, y = pIll, z = pDead, color = method),
  #         size = 1.4
  #       ) +
  #       ggtern::theme_bw() +
  #       ggtern::theme_showarrows() +
  #       ggplot2::labs(
  #         title = paste("Goodman Confidence regions - simulation", sim_str),
  #         x = "Healthy",
  #         y = "Ill",
  #         z = "Dead",
  #       )
  #   })
  # }

  # plot_order <- polygon_df %>%
  #   filter(method %in% c("AJ", "MSMIKM")) %>%
  #   group_by(simulation, time, method) %>%
  #   slice_head(n = 1) %>% # one row per simulation/time/method
  #   ungroup() %>%
  #   select(simulation, time, method, contains_truth) %>%
  #   pivot_wider(
  #     names_from = method,
  #     values_from = contains_truth
  #   ) %>%
  #   mutate(
  #     plot_group = case_when(
  #       AJ & !MSMIKM ~ 1, # AJ only
  #       !AJ & MSMIKM ~ 2, # MSMIKM only
  #       AJ & MSMIKM ~ 3, # both
  #       !AJ & !MSMIKM ~ 4, #neigher
  #       TRUE ~ 5 # (missing)
  #     )
  #   ) %>%
  #   arrange(time, plot_group, simulation)

  # #Save all figures to one PDF
  # out_file <- paste0("confidence_regions_goodman_t1825", cancer_type, ".pdf")
  # sim_ids <- plot_order %>% filter(time == 1825) %>% pull(simulation)

  # grDevices::pdf(out_file, width = 12, height = 8, onefile = TRUE)
  # purrr::walk(sim_ids, \(i) {
  #   p <- plot_sim_regions(i, time_num = 1825)
  #   plot(p)
  # })
  # grDevices::dev.off()

  return()
}

# #Plot true state occupation probabilities for various settings
# truth_df <-
#   imap_dfr(empirical_truth, function(beta_list, cancer_type) {
#     imap_dfr(beta_list, function(df, beta_name) {
#       df %>%
#         filter(time %in% eval_times) %>%
#         mutate(
#           cancer_type = cancer_type,
#           beta = str_remove(beta_name, "^beta")
#         )
#     })
#   }) %>%
#   pivot_wider(
#     names_from = state,
#     values_from = probability
#   ) %>%
#   filter(
#     #cancer_type %in% c("LIHC", "PRAD", "MESO")
#     beta <= 0
#   )

# ggtern(
#   truth_df,
#   aes(
#     x = Healthy,
#     y = Ill,
#     z = Dead,
#     colour = cancer_type,
#     shape = cancer_type,
#     linetype = beta,
#     group = interaction(cancer_type, beta)
#   )
# ) +
#   geom_path(
#     linewidth = 0.8,
#     #arrow = arrow(
#     #  type = "open",
#     #  length = unit(0.08, "in")
#     #)
#   ) +
#   geom_point(size = 2) +
#   theme_bw() +
#   theme_showarrows() +
#   labs(
#     colour = "Cancer Type",
#     shape = "Cancer Type",
#     linetype = "beta",
#     x = "Healthy",
#     y = "Ill",
#     z = "Dead",
#     title = "True State Occupation Probabilities Over Time"
#   )

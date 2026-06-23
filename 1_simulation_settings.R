################################################################################
# Author: Rachel Gonzalez
# Date: 2026-04-16
#
# Purpose:
# Generate data-driven simulation parameters for multiple TCGA cancer types
# using illness-death multi-state models with piecewise exponential transition
# hazards.
#
# The script:
#   1. Loads and preprocesses TCGA survival data
#   2. Fits piece-wise exponential models to each transition hazard
#   3. Fits a censoring distribution for each cancer type
#   4. Stores the resulting "true" simulation parameters
#   5. Generates large simulated datasets to approximate true state
#      occupation probabilities over time
#   6. Creates plots showing model dynamics
#
# Inputs:
#   - tcga_clean.csv
#   - data_driven_simulation_helper.R
#
# Main Outputs:
#   - true_parameters.RData
#   - true_occupation_probabilities.RData
#   - model_dynamics.pdf
################################################################################

# ==============================================================================
# Load required packages
# ==============================================================================

# tidyverse  : data manipulation and visualization
# mstate     : illness-death model utilities
# PWEXP      : piecewise exponential model fitting
# patchwork  : combine ggplots
pacman::p_load(tidyverse, mstate, PWEXP, patchwork)

# Load custom helper functions used for fitting, simulation,
# and state occupation probability estimation.
source("data_driven_simulations/0_data_driven_simulation_helper.R")

# ==============================================================================
# Load and clean TCGA data
# ==============================================================================

# Cancer types included in the analysis
# These will each become separate simulation settings
cancer_types <- c("LIHC", "MESO", "PRAD", "KICH", "KIRC", "UCEC")

tcga <- read_csv("/Users/rachelgonzalez/Documents/Dissertation/Chapter 1/TCGA/tcga_clean.csv") %>%
  filter(type %in% cancer_types) %>%
  # Remove logically inconsistent observations
  filter(ifelse(OS == 0 & PFS == 1 & PFS.time > OS.time, FALSE, TRUE))

tcga <- split(tcga, tcga$type)

# ==============================================================================
# Convert data to long format for multi-state modeling
# ==============================================================================

# Define illness-death transition matrix:
#
# Transition 1: Diagnosis -> Progression
# Transition 2: Diagnosis -> Death
# Transition 3: Progression -> Death
tmat <- trans.illdeath(names = c("Dx", "Progression", "Death"))

# Convert patient-level data into long-format multi-state data
make_long <- function(d, t = tmat) {

  long <- msprep(
    data = d,
    trans = t,

    # Times to progression and death
    time = c(NA, "time_to_progression", "time_to_death"),

    # Event indicators
    status = c(NA, "progression", "death")
  )

  return(long)
}

# Generate long-format data for every cancer type
tcga.long <- map(tcga, make_long)

# Quick check of transition counts and censoring rates
map(tcga.long, events)

# Notes from exploratory analysis:
#
#PRAD: very high censoring (81% do not have either event; very few of those who do progress die)
#KICH: very high censoring (84% do not have either event; over half of those who do progress die)
#UCEC: high censoring (75% do not have either event; 46% of those who progress are observed to die)
#KIRC: medium censoring (58% do not have either event)
#LIHC: lower censoring (38% do not have either event)
#MESO: very low censoring (7% do not have either event)

# ==============================================================================
# Fit generative models for all transitions
# ==============================================================================

# Helper: Fit a piecewise exponential survival model.
#
# Arguments:
#   data   : dataframe containing time and status variables
#   breaks : number of change points
#   plot   : optionally visualize fit
#
# Returns:
#   PWEXP model object
myFit <- function(data, breaks = 3, plot = FALSE) {

  fit <- pwexp.fit(
    data$time,
    data$status,
    nbreak = breaks,
    min_pt_tail = 10
  )

  if (plot) {

    plot_survival(
      data$time,
      data$status,
      main = "Piecewise Exponential Fit"
    )

    plot_survival(
      fit,
      col = "blue",
      lwd = 3,
      show_breakpoint = TRUE
    )

    legend(
      "topright",
      c("TCGA-CDR Data", "piecewise-exponential"),
      lwd = 3,
      col = c("black", "blue")
    )
  }

  return(fit)
}

# Fit piecewise exponential models to each cancer type
  # Default number of breakpoints is 3.
  # Certain cancer types/transitions require fewer breakpoints because
  # the data are too sparse for stable estimation.

# Number of breakpoints used to fit censoring distribution when default=3 is not used
censor_breaks <- c(MESO = 1)

# Structure storing custom breakpoint choices
# for each transition and cancer type when default=3 is not used
breaks <- replicate(
  3,
  setNames(rep(NA_real_, length(cancer_types)), cancer_types),
  simplify = FALSE
)

breaks[[1]][["KICH"]] <- 2

# Container for all fitted model parameters.
truth <- list()

for (cancer in cancer_types) {

  temp_list <- vector("list", 4)

  # Fit transition hazards
  for (t in 1:3) {

    temp_list[[t]] <-
      myFit(
        tcga.long[[cancer]] %>% filter(trans == t),

        breaks = ifelse(
          !is.na(breaks[[t]][cancer]),
          breaks[[t]][cancer],
          3
        ),

        plot = FALSE
      ) %>%
      as.data.frame() %>%
      select(-AIC, -BIC, -likelihood) %>%

      mutate(transition = t)
  }

  # Fit censoring distribution
  # Reverse Kaplan-Meier approach: censoring becomes the event of interest
  censor_fit <- myFit(
    data.frame(
      time = tcga[[cancer]]$OS.time,
      status = 1 - tcga[[cancer]]$OS
    ),
    breaks = ifelse(
      cancer %in% names(censor_breaks),
      censor_breaks[cancer],
      3
    ),
    plot = FALSE
  )

  temp_list[[4]] <-
    censor_fit %>%
    as.data.frame() %>%
    select(-AIC, -BIC, -likelihood) %>%
    mutate(transition = 4)

  # Store all transition and censoring models
  truth[[cancer]] <- bind_rows(temp_list)
}

# Save fitted simulation parameters
save(
  truth,
  file = "data_driven_simulations/true_parameters.RData"
)

# ==============================================================================
# Simulation Diagnostic Checks
# ==============================================================================

# Compare:
#   1. Real TCGA data
#   2. Simulated data generated from fitted model
#   3. Piece-wise exponential truth

# Visual check of the three components above, can change sample size and setting
make_plots <- function(setting, n = 100, seed = NULL) {

  data <- tcga.long[[setting]]

  d.sim <- msprep(
    data = simulate_illness_death(
      n,
      params = truth[[setting]],
      seed = seed
    ),
    trans = tmat,
    time = c(NA, "t1", "t2"),
    status = c(NA, "event1", "event2")
  ) %>%
    group_split(trans)

  for (t in 1:3) {

    sim <- d.sim[[t]]
    real <- data[[t]]

    fit <- pwexp.fit(
      data[[t]]$time,
      data[[t]]$status,
      nbreak = 3,
      min_pt_tail = 10
    )

    title <- paste0(setting, " Transition ", t)

    plot_survival(
      real$time,
      real$status,
      main = title,
      mark.time = FALSE
    )

    plot_survival(
      sim$time,
      sim$status,
      col = "red",
      add = TRUE,
      mark.time = FALSE
    )

    plot_survival(
      fit,
      col = "blue",
      lwd = 3,
      show_breakpoint = TRUE
    )

    legend(
      "bottomleft",
      c(
        "TCGA-CDR data",
        paste0("simulated data, n = ", n),
        "generative model"
      ),
      lwd = 3,
      col = c("black", "red", "blue")
    )
  }
}

# Load random seeds used to generate reproducible simulation examples.
random.seeds <- read_csv(
  "data_driven_simulations/randomSeeds.csv"
)$simulationSeeds

# The following block creates PDF diagnostics for every cancer type and simulations.
# Currently disabled because generating all figures can be time-consuming.

# for (cancer in cancer_types) {
#
# pdf(
# file = paste0(
# "data_driven_simulations/",
# cancer,
# "_sample_simulations.pdf"
# ),
# width = 8,
# height = 3.5
# )
#
# par(mfrow = c(1, 3))
#
# for (i in seq_along(random.seeds)) {
# make_plots(cancer, n = 100, seed = random.seeds[i])
# }
#
# par(mfrow = c(1, 1))
# dev.off()
# }

# ==============================================================================
# Estimate true occupation probabilities
# ==============================================================================

# Maximum follow-up horizon used for each cancer type
max_time <- function(setting) {
  case_when(
    setting %in% c("LIHC", "KIRC", "MESO", "UCEC") ~ 365.25 * 5,
    setting %in% c("PRAD", "KICH") ~ 365.25 * 10,
    TRUE ~ 365.25 * 5
  )
}

# Container for empirical state occupation probabilities.
empirical_truth <- list()

for (setting in names(truth)) {

  t <- max_time(setting)

  # Evaluate two scenarios:
    # beta = 0   -> semi-Markov scenario
    # beta %in% c(-0.75, 0.75) -> extended state-arrival semi-markov scenario
  betas <- c(-0.75, 0, 0.75)

  # Simulate a very large dataset so Monte Carlo error is negligible.
  d.large <- map(betas, function(b) {
    simulate_illness_death(
      n = 1000000,
      params = truth[[setting]],
      return_latent_data = TRUE,
      return_censored_data = FALSE,
      beta = b,
      seed = 123
    )
  })

  names(d.large) <- paste0("beta", betas)

  empirical_truth[[setting]] <- map(
    d.large,
    function(L) {
      get_empirical_probs(
        L,
        times = seq(0, t, 1)
      ) %>%
        select(time, pHealthy, pIll, pDead) %>%
        pivot_longer(
          cols = c("pHealthy", "pIll", "pDead"),
          names_to = "state",
          names_prefix = "p",
          values_to = "probability"
        )
    }
  )

  names(empirical_truth[[setting]]) <- paste0("beta", betas)
}

# Save estimated occupation probabilities.
save(
  empirical_truth,
  file = "data_driven_simulations/true_occupation_probabilities.RData"
)


# ==============================================================================
# Create state occupation probability plots
# ==============================================================================

load("data_driven_simulations/true_occupation_probabilities.RData")
load("data_driven_simulations/true_parameters.RData")

# Convert nested list structure into one plotting dataframe
empirical_truth_long <- purrr::imap_dfr(
  empirical_truth,
  ~ bind_rows(.x, .id = "setting") %>%
    mutate(cancer_type = .y)
) %>%
  mutate(
    state = factor(
      state,
      levels = c("Healthy", "Ill", "Dead"),
      ordered = TRUE
    )
  )

# Generate state occupation probability plot for a specific cancer type and beta value
plot_truth <- function(cancer, beta) {

  t <- max_time(cancer)

  df <- filter(
    empirical_truth_long,
    cancer_type == cancer,
    setting == paste0("beta", as.character(beta))
  ) %>%
    select(time, probability, state)

  ggplot(df) +
    geom_line(
      aes(
        x = time,
        y = probability,
        color = state,
        linetype = state
      ),
      linewidth = 3
    ) +
    labs(
      title = paste0(cancer, " HR(T1) = ", round(exp(beta), 2)),
      x = "Time (years)",
      y = "Occupation Probability"
    ) +
    scale_color_brewer(type = "qual", palette = 1) +
    theme_minimal(base_size = 25) +
    scale_x_continuous(
      breaks = seq(0, t, 365.25),
      labels = seq(0, t / 365.25),
      limits = c(0, t)
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90")
    )
}

# Create every cancer:beta combination.
values <- crossing(
  beta = c(-0.75, 0, 0.75),
  cancer = names(truth)
)

truth_plots <- pmap(values, plot_truth) %>%
  setNames(
    values %>%
      transmute(id = paste0(cancer, "_beta", beta)) %>%
      pull(id)
  )

# Add a table showing the proportion still at risk
# for progression and death at selected time points.
add_risk_table <- function(
    p,
    cancer_type,
    beta,
    year_interval = 1) {

  # Generate large censored dataset.
  d <- simulate_illness_death(
    n = 1000000,
    params = truth[[cancer_type]],
    beta = beta,
    return_latent_data = FALSE,
    return_censored_data = TRUE,
    seed = 123
  )

  # Extract x-axis information directly from plot.
  x_breaks <- ggplot_build(p)$layout$panel_params[[1]]$x$breaks
  x_limits <- ggplot_build(p)$layout$panel_params[[1]]$x.range

  # Optionally reduce number of displayed time points.
  x_breaks_filtered <-
    x_breaks[seq(1, length(x_breaks), by = year_interval)]

  # Calculate risk quantities at each displayed time.
  t.risk <- map_dfr(x_breaks_filtered, function(x) {

    denom <- nrow(d)

    tibble(
      time = x,
      risk_prog =
        nrow(filter(d, t1 > x)) / denom,

      risk_death =
        nrow(filter(d, t2 > x)) / denom,

      cum_prog =
        nrow(filter(d, t1 <= x & event1 == 1)) / denom,

      cum_death =
        nrow(filter(d, t2 <= x & event2 == 1)) / denom
    )
  })

  # Reshape for plotting.
  t.risk.long <- t.risk %>%
    pivot_longer(
      cols = -time,
      names_to = c("metric", "event_type"),
      names_sep = "_"
    ) %>%
    pivot_wider(
      names_from = metric,
      values_from = value
    ) %>%
    mutate(
      label = sprintf("%.2f", risk),

      event_type = factor(
        event_type,
        levels = c("prog", "death"),
        labels = c("Progression", "Death")
      )
    )

  # Plot risk table.
  p.risk <- ggplot(
    t.risk.long,
    aes(x = time, y = event_type, label = label)
  ) +
    geom_text(size = 8) +
    scale_x_continuous(
      breaks = x_breaks_filtered,
      labels = x_breaks_filtered / 365.25,
      limits = x_limits
    ) +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 25) +
    theme(
      plot.margin = margin(t = 0)
    ) +
    ggtitle("Proportion at Risk")

  # Combine main figure and risk table.
  p / p.risk +
    plot_layout(heights = c(3, 1))
}

# Custom spacing between displayed years.
year_intervals <- c(
  LIHC = 1,
  PRAD = 2,
  MESO = 1,
  KICH = 2,
  KIRC = 1,
  UCEC = 1
)

# Generate final plots.
plots <- vector("list", length = nrow(values))

for (j in 1:nrow(values)) {

  cancer <- as.character(values[j, "cancer"])
  beta <- as.numeric(values[j, "beta"])

  plots[[j]] <- add_risk_table(
    p = truth_plots[[paste0(cancer, "_beta", beta)]],
    cancer_type = cancer,
    beta = beta,
    year_interval = year_intervals[cancer]
  )
}

# Exportfigure collection


pdf(
  "data_driven_simulations/model_dynamics.pdf",
  height = 8,
  width = 10
)

map(plots, function(p) {
  p
})

dev.off()

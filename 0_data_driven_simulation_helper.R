#Purpose: Build helper functions for TGCA data-driven MSMI simulations
#Author: Rachel Gonzalez
#Date: April 20 2026

# General Purpose Functions -----------------------------------------------

#Wrapper for sample that gives the desired behavior when n=1
resample <- function(x, ...) x[sample.int(length(x), ...)]


#' Return value of a step function at a given time
#'
#' @param step_times a numeric vector of times where the step function jumps
#' @param step_values a numeric vector of values that correspond to a step function evaluated at the step times
#' @param t a time point of interest, usually not a value of step_times
#'
#' @returns the value of the step function at time t
get_step_value <- function(step_times, step_values, t) {

  idx <- max(which(step_times <= t), na.rm = TRUE)

  if (is.na(idx)) {
    return(NA)
  }

  step_values[idx]
}


#' Calculate distance between two discrete probability distributions
#'
#' @param df_reference data frame with columns corresponding to the values of the reference distribution over evaluation times
#' @param df_estimate data frame with columns corresponding to the values of the distribution to compare over evaluation times
#' @param cols_reference character vector containing the column names defining the probability distribution for `df_reference`
#' @param cols_estimate character vector containing the column names defining the probability distribution for `df_estimate`
#' @param times numeric vector containing time points when distance calculations are requested
#' @param epsilon small positive numeric value for computational stability (default = 1e-12)
#'
#' @returns a dataframe with the KL and Helinger distance at each time point
compute_distance <- function(
  df_reference,
  df_estimate,
  cols_reference = c("pHealthy", "pIll", "pDead"),
  cols_estimate = c("pstate1", "pstate2", "pstate3"),
  times,
  epsilon = 1e-12
) {
  distance_df <- lapply(times, function(t) {
    p_ref <- as.numeric(df_reference[
      which(df_reference$time == t),
      cols_reference
    ]) +
      epsilon
    p_est <- as.numeric(df_estimate[
      which(df_estimate$time == t),
      cols_estimate
    ]) +
      epsilon

    distance_df = tibble(
      time = t,
      KL = philentropy::KL(rbind(p_est, p_ref), unit = "log2"),
      hellinger = (1 / sqrt(2)) * sqrt(sum((sqrt(p_est) - sqrt(p_ref))^2))
    )
  })

  return(bind_rows(distance_df))
}
# Data Generation Process -------------------------------------------------


#' Simulate piece-wise exponential data
#'
#' @param params dataframe with piecewise exponential parameters for each transition
#' @param trans transition of interest (numeric)
#' @param n sample size (numeric)
#'
#' @returns a numeric vector with random draws from the specified PW exponential distribution
simPW <- function(params, trans, n) {
  temp <- params %>% filter(transition == trans)

  myRates <- temp %>% select(starts_with("lam")) %>% as.numeric() %>% na.omit()
  myBreaks <- temp %>% select(starts_with("brk")) %>% as.numeric() %>% na.omit()
  s <- PWEXP::rpwexp(n = n, rate = myRates, breakpoint = myBreaks)

  return(s)
}


#' Simulate Cox-Pieceiwise exponential survival times
#'
#' @param params dataframe with the piece-wise exponential parameters for each transition
#' @param beta true parameters for the cox model, logHR associated with columns of `X` (px1)
#' @param X Covariate matrix (nxp)
#'
#' @returns a length n vector with random survival times according to the coviariates X, associated logHR beta, and parameters for the baseline hazard
simCoxPW <- function(params, beta, X) {

  lambda <- params %>% filter(transition == 3) %>%
    select(starts_with("lam")) %>% as.numeric()
  cuts <- params %>% filter(transition == 3) %>%
    select(starts_with("brk")) %>% as.numeric()

  #Sample random survival times based on covariate(s) X
  n = nrow(X)
  times <- numeric(n)

  #Transform log(HR) for yearly time scale to days time scale since X is specified in days
  beta = beta/365

  for (i in seq_len(n)) {

    eta <- sum(X[i, ] * beta)
    r <- exp(eta)

    U <- runif(1)
    Z <- -log(U) / r

    # inverse cumulative hazard function at Z, using H^-1(x) = S^-t(exp(-x))
    times[i] <- qpwexp(1-exp(-Z), lambda, cuts)

  }
  times
}

# #Check that SimCoxPW functions are the same as rpwexp when beta = 0
# params <- truth$LIHC
# trans = 3
# beta = 0
# X = as.matrix(rnorm(100000), nrow = 1)
# n = 100000
#
# test = data.frame(pkg = simPW(params, trans, n),
#                   myMethod = simCoxPW(params, beta, X, seed = 123)) %>%
#   pivot_longer(cols = c("pkg", "myMethod"), values_to = "value", names_to = "method")
#
# ggplot(test) + geom_density(aes(x = value, color = method), alpha = 0.5)


#' Simulate data from a multi-state model with three states: healthy (0), ill (1), dead (2)
#'
#' @param n desired sample size (numeric)
#' @param params dataframe with piecewise exponential parameters for each transition
#' @param beta logHR for time from illness to death associated with time of arrival to the illness state
#' @param return_latent_data indicates if the function should return the underlying survival times (\tilde{X}) (logical)
#' @param return_censored_data indicates if the function should return the censored data that is observed in practice (X). (logical)
#' @param seed random seed for reproducibility. Default is a random seed. (numeric)
#'
#' @returns dataframe with n rows and columns id, t1, event1, t2, event2 corresponding to event times and indicators for time to illness and time to death,
#' or a list of dataframes if `return_latent_data` and `return_censored_data` are both TRUE.
simulate_illness_death <- function(
  n,
  params,
  beta = NA,
  return_latent_data = FALSE,
  return_censored_data = TRUE,
  seed = NULL
) {
  #Set up random seed for data generation
  if (is.null(seed)) {
    warning(
      "Please provide a seed for reproducibility. A random seed has been generated."
    )
    seed = sample(1:.Machine$integer.max, 1)
  }

  set.seed(seed)

  #Determine first event (illness or death)
  sojourn12 <- simPW(params, t = 1, n)
  sojourn13 <- simPW(params, t = 2, n)

  #Generate time from illness to death only for those with illness
  sojourn23 <- rep(NA, n)
  ill_ids <- which(sojourn12 < sojourn13)
  if (is.na(beta) | beta == 0) {
    #sojourn23[ill_ids] <- simPW(params, t = 3, n = length(ill_ids))
    X = as.matrix(sojourn12[ill_ids],ncol = 1)
    sojourn23[ill_ids] <- simCoxPW(params, beta = 0, X)
  } else {
    X = as.matrix(sojourn12[ill_ids],ncol = 1)
    sojourn23[ill_ids] <- simCoxPW(params, beta, X)
  }

  temp = data.frame(
    sojourn12 = sojourn12,
    sojourn13 = sojourn13,
    sojourn23 = sojourn23
  )

  #Time of entry into each state and event indicators
  temp = temp %>%
    mutate(
      id = 1:n,
      t1 = if_else(sojourn12 < sojourn13, sojourn12, sojourn13),
      event1 = if_else(sojourn12 < sojourn13, 1, 0),
      t2 = if_else(sojourn12 < sojourn13, sojourn12 + sojourn23, sojourn13),
      event2 = 1
    ) %>%
    select(id, t1, event1, t2, event2)

  #Add in censoring times
  temp_censoring <- temp %>%
    mutate(
      C = simPW(params, t = 4, n),
      #C = runif(n, min = 0, max = 5000),
      event1 = if_else(C < t1, 0, event1),
      t1 = if_else(C < t1, C, t1),
      event2 = if_else(C < t2, 0, event2),
      t2 = if_else(C < t2, C, t2)
    )

  temp_censoring <- temp_censoring %>%
    select(id, t1, event1, t2, event2)

  if (return_latent_data == FALSE & return_censored_data == TRUE) {
    return(d.observed = temp_censoring)
  } else if (return_censored_data == FALSE & return_latent_data == TRUE) {
    return(d.latent = temp)
  } else {
    return(list(
      d.observed = data.frame(temp_censoring),
      d.latent = data.frame(temp)
    ))
  }
}

# Empirical Probability in State ------------------------------------------

#' Calculate state occupation probabilities from a dataset without censoring
#'
#' @param df dataframe with columns, t1, event1, t2, and event2 that describe data from an illness-death model
#' @param times times to calculate state occupation probabilities
#'
#' @returns a dataframe with with columns corresponding to empirical state occupation probabilities at each time
get_empirical_probs <- function(df, times) {
  # Initialize an empty list to store results for each time
  results <- lapply(times, function(t) {
    # logical conditions for each state
    is_healthy <- df$t1 > t & df$t2 > t
    is_ill <- df$t1 <= t & df$t2 > t
    is_dead_no_illness <- df$t2 <= t & df$event1 == 0
    is_dead_with_illness <- df$t2 <= t & df$event1 == 1

    # Calculate proportions in each state
    n <- nrow(df)
    pHealthy <- sum(is_healthy) / n
    pIll <- sum(is_ill) / n
    pDeadMinusIll <- sum(is_dead_no_illness) / n
    pDeadWithIll <- sum(is_dead_with_illness) / n
    pDead <- pDeadMinusIll + pDeadWithIll

    tibble(
      time = t,
      pHealthy = pHealthy,
      pIll = pIll,
      pDeadMinusIll = pDeadMinusIll,
      pDeadWithIll = pDeadWithIll,
      pDead = pDead
    )
  })

  bind_rows(results)
}


# Empirical Transition Probability ----------------------------------------


#' Calculate empirical transition probabilities
#'
#' @param df dataframe with columns, t1, event1, t2, and event2 that describe data from an illness-death model
#' @param s  starting time >=0, fixed
#' @param from starting state, fixed
#' @param times times to evaluate transition probabilities, must be greater than `s`
#'
#' @returns a dataframe with rows corresponding to different time points and columns corresponding to arrival states k = {1,2,3}
empirical_transition_probs <- function(df, s, from, times) {
  results <- map_dfr(times, function(t) {
    # subsets representing people in state=from at time=s
    if (from == 1) {
      is_healthy <- df$t1 > s & df$t2 > s
      sub <- df[is_healthy, ]
    } else if (from == 2) {
      is_ill <- df$t1 <= s & df$t2 > s
      sub <- df[is_ill, ]
    } else if (from == 3) {
      is_dead <- df$t2 <= s
      sub <- df[is_dead, ]
    }

    # Calculate proportions in each state at time t from subgroup
    is_state1 <- sub$t1 > t & sub$t2 > t
    is_state2 <- sub$t1 <= t & sub$t2 > t
    is_state3 <- sub$t2 <= t

    n <- nrow(sub)
    pHealthy <- sum(is_state1) / n
    pIll <- sum(is_state2) / n
    pDead <- sum(is_state3) / n

    tibble(
      stime = s,
      time = t,
      pstate1 = pHealthy,
      pstate2 = pIll,
      pstate3 = pDead
    )
  })
}

#' Calculate state transition probability matrix from all possible starting states
#'
#' @param df dataframe with columns, t1, event1, t2, and event2 that describe data from an illness-death model
#' @param s  starting time >=0, fixed
#' @param times times to evaluate transition probabilities, must be greater than `s`
#'
#' @returns a list of dataframes rows corresponding to different time points and columns corresponding to arrival states k = {1,2,3}
#'  over all possible starting states
get_empirical_transition_probs <- function(df, s, times) {
  from = c(1, 2, 3) #states specific to illness death model
  names(from) = c("Healthy", "Ill", "Dead")

  L_results <- map(from, function(f) {
    empirical_transition_probs(from = f, df = df, s = s, times = times)
  })

  return(L_results)
}

# Functions for Checking Coverage -----------------------------------------


#' Multinomial logit transformation
#'
#' @param p a length k vector of probabilities
#'
#' @returns a length k-1 vector of real numbers
multinomial_logit <- function(p) {
  eps <- 1e-8
  p[p == 0] <- eps
  p = p / sum(p)

  k <- length(p)
  x <- numeric(k - 1)

  for (i in 1:(k - 1)) {
    x[i] <- log(p[i] / p[k])
  }

  return(x)
}

#' Inverse Multinomial Logit Transformation
#'
#' @param x a vector of real numbers of length k-1
#'
#' @returns a length k vector of probabilities
multinomial_logit_inverse <- function(x) {
  k <- length(x) + 1
  p <- numeric(k)

  denom <- 1 + sum(exp(x))

  for (i in 1:(k - 1)) {
    p[i] <- exp(x[i]) / denom
  }
  p[k] <- 1 / denom

  return(p)
}

#' Compute Jacobian Matrix for the Multinomial Logit Transformation
#'
#' @param p vector of probabilities of length k
#'
#' @returns the (k-1, k) dimensional Jacobian matrix
jacobian <- function(p) {
  J <- matrix(0, nrow = 2, ncol = 3)
  J[1, 1] <- 1 / p[1]
  J[1, 3] <- -1 / p[3]
  J[2, 2] <- 1 / p[2]
  J[2, 3] <- -1 / p[3]
  return(J)
}


#' Generate Confidence Region for Aalen-Johansen estimates of state occupation probabilities
#'
#' @param p Aalen-Johansen estimates of the state occupation probabilities
#' @param Sigma covariance matrix associated with the AJ estimates
#' @param alpha type 1 error, default (0.05)
#'
#' @returns a dataframe with points in the 2-simplex representing the boundary of the 1-alpha Wald confidence region
AJ_Region <- function(p, Sigma, alpha = 0.05) {

  z <- multinomial_logit(p)
  J <- jacobian(p)
  Sigma_z <- J %*% Sigma %*% t(J)

  if (any(!is.finite(Sigma_z))) {
    return(NULL)
  }
  if (is.null(p) || is.null(Sigma_z)) {
    return(NULL)
  }

  ellipse_pts_z <- mixtools::ellipse(
    mu = z,
    sigma = Sigma_z,
    alpha = alpha,
    npoints = 500,
    draw = FALSE
  )

  ellipse_pts <- t(apply(ellipse_pts_z, 1, multinomial_logit_inverse))

  return(ellipse_pts)
}

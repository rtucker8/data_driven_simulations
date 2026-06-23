#Author: Rachel Gonzalez
#Date: 23 January 2026
#Purpose: Create R script that will run a host of simulation settings on the cluster

library(tidyverse)

source("data_driven_simulations/2_data_driven_simulations.R")

params <- crossing(
  cancer_type = c("LIHC", "MESO", "PRAD", "UCEC", "KIRC", "KICH"),
  sample_size = c(50, 100, 200, 500),
  n_imps = c(10, 20, 30, 50),
  n_sims = 500,
  beta = c(-0.75, 0)
)

args <- commandArgs(trailingOnly = TRUE)

task_id <- as.integer(args[1])

if (is.na(task_id)) {
  stop("No task ID supplied")
}

setting <- params[task_id, ]

run_sim(
  sample_size = setting$sample_size,
  n_sims      = setting$n_sims,
  n_imps      = setting$n_imps,
  cancer_type = setting$cancer_type,
  beta        = setting$beta
)

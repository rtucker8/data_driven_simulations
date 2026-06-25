#Purpose: Compare methods to compute state occupation probabilities from an illness-death model
#Author: Rachel Gonzalez
#Date: May 26, 2026

#Input: data_driven_simulations/./coverage_data_XXX.csv, csv files containing the coverage probabilities for each setting

#Libraries
library(tidyverse)

# Load simulation results ---------------------------------------------------------------

coverage_data <- list.files(
  path = "Results/trWald/kich",
  pattern = "^coverage_data_.*\\.csv$",
  full.names = TRUE,
  recursive = TRUE
) %>%
  set_names() %>%
  map(read_csv) %>%
  list_rbind(names_to = "filename") %>%
  mutate(
    file = basename(filename),
    cancer_type = factor(str_extract(file, "(?<=coverage_data_)[^_]+(?=_n)")),
    n = factor(str_extract(file, "(?<=_n)\\d+(?=_M)")) %>% fct_inseq(),
    M = factor(str_extract(file, "(?<=_M)\\d+(?=_b)")) %>% fct_inseq(),
    beta = factor(str_extract(file, "(?<=_b)-?\\d*\\.?\\d+(?=_)")) %>% fct_inseq(),
    approach = factor(str_extract(file, "[^_]+(?=\\.csv$)"))
  ) %>%
  select(-filename, -file)


#Make coverage plot for each combination of cancer type and beta
coverage_groups <- coverage_data %>%
  filter(cancer_type %in% c("MESO", "LIHC", "PRAD"),
         approach = "trWald") %>%
  group_by(cancer_type, beta) %>%
  pivot_longer(cols = c("coverageAJ", "coverageMSMI", "coverageMSMICox"), names_to = "method", names_prefix = "coverage", values_to = "coverage")

group_data <- coverage_groups %>% group_split()
group_keys <- coverage_groups %>% group_keys()

for (i in seq_along(group_data)) {

  dat <- group_data[[i]]
  cancer <- group_keys$cancer_type[i]
  beta <- group_keys$beta[i]

  y_min <- min(min(dat$coverage), 0.8)

  p <- ggplot(
    dat,
    aes(
      x = time / 365,
      y = coverage,
      color = method,
      shape = method
    )
  ) +
    geom_point() +
    geom_line() +
    geom_hline(
      yintercept = 0.95,
      color = "black",
      linetype = "dashed"
    ) +
    scale_x_continuous(breaks = seq(0, 10, 1)) +
    scale_y_continuous(limits = c(y_min, 1)) +
    facet_grid(n ~ M, labeller = label_both) +
    theme_bw(base_size = 15) +
    labs(
      title = paste0(
        "Coverage for ", cancer,
        " (beta = ", beta, ")"
      ),
      x = "Time (years)",
      y = "Empirical Coverage Probability",
      color = "Method",
      shape = "Method"
    )

  ggsave(
    filename = paste0(
      "Output/Figures/",
      tolower(cancer),
      "coverage_",
      cancer,
      "_beta_",
      beta,
      ".pdf"
    ),
    plot = p,
    height = 8,
    width = 10
  )
}

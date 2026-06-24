#Purpose: Compare methods to compute state occupation probabilities from an illness-death model
#Author: Rachel Gonzalez
#Date: May 26, 2026

#Input: data_driven_simulations/coverage_data_XXX.csv, csv files containing the coverage probabilities for each setting

#Libraries
library(tidyverse)

#TO DO: update to accommodate different betas

# Load simulation results ---------------------------------------------------------------
coverage_data <- list.files(
  path = "Output",
  pattern = "^coverage_data_.*\\.csv$",
  full.names = TRUE
) |>
  set_names() |>
  map(read_csv) |>
  list_rbind(names_to = "filename") |>
  mutate(
    n = factor(str_extract(filename, "(?<=_)(\\d+)(?=_size_)")) %>% fct_inseq(),
    M = factor(str_extract(filename, "(?<=_size_)(\\d+)(?=_imps)")) %>%
      fct_inseq(),
  ) |>
  select(-filename)

coverage_by_cancer <- coverage_data |>
  filter(cancer_type %in% c("MESO", "LIHC", "PRAD")) |>
  group_by(cancer_type) |>
  group_split() |>
  set_names(sort(c("MESO", "LIHC", "PRAD")))
  #set_names(unique(sort(coverage_data$cancer_type)))


#Make plots that compare the coverage across different settings
for (cancer in names(coverage_by_cancer)) {
  y_min = pmin(min(coverage_by_cancer[[cancer]]$coverage), 0.8)
  ggplot(
    data = coverage_by_cancer[[cancer]],
    aes(x = time / 365, y = coverage, color = method, shape = method)
  ) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.950, color = "black", linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 10, 1)) +
    scale_y_continuous(limits = c(y_min, 1)) +
    facet_grid(n ~ M, labeller = label_both) +
    theme_bw(base_size = 15) +
    labs(
      title = paste0("Coverage for ", cancer, " Setting"),
      x = 'Time (years)',
      y = "Empirical Coverage Probability",
      color = "Method",
      shape = "Method"
    )

  ggsave(
    paste0("Output/Figures/coverage_", cancer, ".pdf"),
    height = 8,
    width = 10
  )
}

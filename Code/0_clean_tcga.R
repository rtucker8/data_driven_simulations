#Author: Rachel Gonzalez
#Date:   March 2026
#Purpose: Perform data cleaning for TCGA-CDR Raw Data

# Load necessary libraries
library(tidyverse, readxl)

# Read in the data
d.all <- read_excel("Data/tcga_raw.xlsx", sheet = "TCGA-CDR")[,
  -1
]
more.edpoints <- read_excel(
  "Data/tcga_raw.xlsx",
  sheet = "ExtraEndpoints"
)[, -1]

d.all <- d.all %>%
  left_join(more.edpoints, by = c("bcr_patient_barcode", "type"))
rm(more.edpoints)

#Exclude those with no followup
d.all <- d.all %>%
  filter(
    !(OS.time == 0 | PFS.time == 0) &
      !is.na(OS.time) &
      !is.na(PFS.time)
  ) #exclude those with no followup

#Tease out semi-competing risks data using time to progression and time to death)
d.all <- d.all %>%
  mutate(
    time_to_progression = PFS.time,
    progression = if_else(PFS.time != OS.time & PFS == 1, 1, 0),
    time_to_death = OS.time,
    death = OS
  )

#Save cleaned dataset
write_csv(d.all, "Data/tcga_clean.csv")


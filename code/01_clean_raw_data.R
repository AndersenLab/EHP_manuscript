#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load effect data taken from 'master' tab in Master_File.xlxs
raw_dat <- readxl::read_excel("data/Master_File.xlsx", sheet = 'master', na = c("NA", ""))

# clean it
proc_dat <- raw_dat %>%
  dplyr::select(-Standard.Error) %>% # not needed b/c all NA
  janitor::clean_names(case = "snake") %>% # fix var names to all be snake case
  dplyr::select(cas, chem_name, group, latin_name, strain, test_type, test_statistic, duration_h, effect_value) # reorder for clarity

# save it
readr::write_csv(proc_dat, file = "data/01_clean.csv")

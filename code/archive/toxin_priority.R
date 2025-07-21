#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load cleaned data
dat <- data.table::fread("data/processed/00_data.csv")

# get the Boyd toxins with high effect 
boyd1 <- dat %>%
  dplyr::filter(source == "Boyd et al. 2016") %>%
  dplyr::mutate(log_effect = log10(effect_value))

ggplot(boyd1) +
  aes(x = log_effect, color = QC) +
  geom_histogram()

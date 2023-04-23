#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/03_clean.csv")

#==================================================#
# Part 1: Setup the priority Orth Regs
#==================================================#
# Widmayer vs Oncorhynchus mykiss 96 hr LC50
d1 <- dat %>%
  dplyr::filter((latin_name == "Oncorhynchus mykiss" &
                   duration_d == 4 &
                   test_statistic == "LC50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2))

# run the orth regressions across all pairs
debug(pwOrthReg)
d1l <- pwOrthReg(data = d1, group = "latin_name", message = T)
d1df <- data.table::rbindlist(d1l)

                  

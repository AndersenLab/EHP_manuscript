#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data
#dat <- data.table::fread("data/processed/03_clean.csv")
dat <- data.table::fread("data/processed/00_data.csv")
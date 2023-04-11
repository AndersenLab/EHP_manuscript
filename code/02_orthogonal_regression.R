#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data
dat <- data.table::fread("data/processed/03_clean.csv") 

#==================================================#
# Part 1: loop over group pairs with the pairwise 
# orthogonal regression function. This function will
# loop over all pairwise comparisions found in the data
#==================================================#
# Do pariwise orthogonal regression across all groups - 330 unique groups = 54285 unique pairs!
# TAKES ~ 20 min
pw.orth.reg.1 <- pwOrthReg(data = dat, group = "group", message = T)
pw.orth.reg.1.df <- data.table::rbindlist(pw.orth.reg.1)

# look at the failed pairs
failed.1.df <- pw.orth.reg.1.df %>%
  dplyr::filter(is.na(orth.reg.slope))

# look at the good pairs
good.1.df <- pw.orth.reg.1.df %>%
  dplyr::filter(!(is.na(orth.reg.slope)))

best.test <- orthReg(data = dat,
                     x = c("NEMATODE", "EC10", "1"),
                     y = c("NEMATODE", "EC90", "1"),
                     plot = T)
best.test$plot

best.test.dat <- dat %>%
  dplyr::filter((group == "NEMATODE" & test_statistic == "EC10" & duration_d == "1") | (group == "NEMATODE" & test_statistic == "EC90" & duration_d == "1"))

best.test.2 <- orthReg(data = best.test.dat,
                     x = c("NEMATODE", "EC10", "1"),
                     y = c("NEMATODE", "EC90", "1"),
                     plot = T)
best.test.2$plot

# look at failed pair example in raw data - looks like a single value for FISH causes error
failed <- dat %>%
  dplyr::filter((group == "NEMATODE" &  test_statistic == "EC10" & duration_d == "1") | (group == "FISH" &  test_statistic == "LD50" & duration_d == "7"))
test2 <- pwOrthReg(data = failed, group = "group", message = FALSE)
test3 <- pwOrthReg(data = failed, group = "group", message = FALSE)
orth.reg.data.out.test2 <- data.table::rbindlist(test2)

#### OLF CODE ####
# OLD function testing code
# #====================================================#
# # Part 1: test the function for latin_name
# #====================================================#
# test1 <- orthReg(data = dat,
#         x = c("Oncorhynchus mykiss", "LC50", "4"),
#         y = c("Daphnia magna", "EC50", "2"), plot = T)
# 
# # # see the filtered focal data
# # a <- test1$all_focal_dat
# # 
# # # see the reshaped plot data
# # b <- test1$plot_dat
# 
# # see the orthogonal regression model
# #test1$orthog_reg_model
# 
# # see the plot
# test1$plot
# 
# # save a test plot
# cowplot::ggsave2(test1$plot, filename = "plots/rainbow_vs_daphnia.png", width = 6, height = 5)
# 
# #====================================================#
# # Part 2: test the function for group
# #====================================================#
# test2 <- orthReg(data = dat,
#                  x = c("FISH", "LC50", "4"),
#                  y = c("INVERT", "EC50", "2"), plot = T)
# 
# # see the plot
# test2$plot
# 
# # save a test plot
# cowplot::ggsave2(test2$plot, filename = "plots/FISH_vs_INVERT.png", width = 6, height = 5)
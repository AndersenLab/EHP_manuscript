#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data
dat <- data.table::fread("data/processed/03_clean.csv") 

#====================================================#
# Part 1: test the function for latin_name
#====================================================#
test1 <- orthReg(data = dat,
        x = c("Oncorhynchus mykiss", "LC50", "4"),
        y = c("Daphnia magna", "EC50", "2"))

# # see the filtered focal data
# a <- test1$all_focal_dat
# 
# # see the reshaped plot data
# b <- test1$plot_dat

# see the orthogonal regression model
#test1$orthog_reg_model

# see the plot
test1$plot

# save a test plot
cowplot::ggsave2(test1$plot, filename = "plots/rainbow_vs_daphnia.png", width = 6, height = 5)

#====================================================#
# Part 2: test the function for group
#====================================================#
test2 <- orthReg(data = dat,
                 x = c("FISH", "LC50", "4"),
                 y = c("INVERT", "EC50", "2"))

# see the plot
test2$plot

# save a test plot
cowplot::ggsave2(test2$plot, filename = "plots/FISH_vs_INVERT.png", width = 6, height = 5)

#==================================================#
# Test pulling the mse, rmese
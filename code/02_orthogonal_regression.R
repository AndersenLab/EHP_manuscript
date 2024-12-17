#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data
#dat <- data.table::fread("data/processed/03_clean.csv")
dat <- data.table::fread("data/processed/00_data.csv")
#===================================================================================================#
# Part 1: loop over group pairs with the pairwise orthogonal regression function. This function will
# loop over all pairwise comparisons found in the data. There are over 2 million unique pairs if we
# consider species. If we just consider groups there are ~160,000 pairs. Let's look at groups first,
# then I suggest we look at C. elegans and the most important species only.
#===================================================================================================#
pw.orth.reg.list <- pwOrthReg(data = dat, group = "group", message = T)
pw.orth.reg.df <- data.table::rbindlist(pw.orth.reg.list, fill = T)
# Save orth reg data b/c it takes a long time to generate - 1min per 1K pairs
save(pw.orth.reg.df, file = glue::glue("data/processed/{today}.pw.orth.reg.df.rda"))
#load orth reg data (UPDATE TO MOST RECENT IF NEEDED)
load("data/processed/20241029.pw.orth.reg.df.rda")

#==================================================#
# Part 2: Summarize the results from pwOrthReg
#==================================================#
# filter the pairwise regressions to those with >= 5 observations
pw.orth.reg.df.proc <- pw.orth.reg.df %>%
  dplyr::filter(orth.reg.n.observations >= 10) 

# export those regressions
#save(pw.orth.reg.df.proc, file = "data/processed/20230412.pw.orth.reg.df.proc.rda")

# reshape for plotting
pw.orth.reg.df.proc.long <- pw.orth.reg.df.proc %>%
  tidyr::pivot_longer(cols = c(orth.reg.n.observations:orth.reg.r.squared))

# now plot it
p1 <- ggplot(pw.orth.reg.df.proc.long) +
  aes(x = value) +
  geom_histogram(bins = 50) +
  facet_wrap(~name, scales = "free") +
  labs(x = "", subtitle = glue::glue("{nrow(pw.orth.reg.df.proc)} of {nrow(pw.orth.reg.df)} possible orthogonal regressions\nwith >= 10 observations"))
p1  

# save the plot
cowplot::ggsave2(p1, filename = glue::glue("plots/{today}_orth.reg.model.stat.dist.png", width = 7.5, height = 7.5)

#==================================================#
# Part 3: Prioritize elegans comps
#==================================================#
n_orthreg <- pw.orth.reg.df %>%
  dplyr::filter(grepl(x, pattern = "NEMATODE") | grepl(y, pattern = "NEMATODE")) %>%
  dplyr::arrange(orth.reg.slope, orth.reg.n.observations)



#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data
dat <- data.table::fread("data/processed/03_clean.csv")

#===================================================================================================#
# Part 1: loop over group pairs with the pairwise orthogonal regression function. This function will
# loop over all pairwise comparisons found in the data. There are over 2 million unique pairs if we
# consider species. If we just consider groups there are ~160,000 pairs. Let's look at groups first,
# then I suggest we look at C. elegans and the most important species only.
#===================================================================================================#
#pw.orth.reg.list <- pwOrthReg(data = dat, group = "group", message = T)
#pw.orth.reg.df <- data.table::rbindlist(pw.orth.reg.list)
# # Save orth reg data b/c it takes a long time to generate - 1min per 1K pairs
#save(pw.orth.reg.df, file = "data/processed/20230412.pw.orth.reg.df.rda")
# load orth reg data
load("data/processed/20230412.pw.orth.reg.df.rda")

#==================================================#
# Part 2: Summarize the results from pwOrthReg
#==================================================#
# filter the pairwise regressions to those with >= 5 observations
pw.orth.reg.df.proc <- pw.orth.reg.df %>%
  dplyr::filter(orth.reg.n.observations >= 5) 

# export those regressions
save(pw.orth.reg.df.proc, file = "data/processed/20230412.pw.orth.reg.df.proc.rda")

# reshape for plotting
pw.orth.reg.df.proc.long <- pw.orth.reg.df.proc %>%
  tidyr::pivot_longer(cols = c(orth.reg.n.observations:orth.reg.r.squared))

# now plot it
p1 <- ggplot(pw.orth.reg.df.proc.long) +
  aes(x = value) +
  geom_histogram(bins = 50) +
  facet_wrap(~name, scales = "free") +
  labs(x = "", subtitle = glue::glue("{nrow(pw.orth.reg.df.proc)} of {nrow(pw.orth.reg.df)} possible orthogonal regressions\nwith >= 5 five observations"))
p1  

# save the plot
cowplot::ggsave2(p1, filename = "plots/orth.reg.model.stat.dist.png", width = 7.5, height = 7.5)

#==================================================#
# check 1: explore NAs for orth.reg stats
#==================================================#
# look at fraction of pairs without reg data
table(is.na(pw.orth.reg.df$orth.reg.slope))

# setup df to help understand why there are so many NAs. choose three random to explore
set.seed(99)
na.df <- pw.orth.reg.df %>%
  dplyr::filter(is.na(orth.reg.slope)) %>%
  dplyr::slice_sample(n = 1000)

# setup with id in raw data
id.dat <- dat %>%
  dplyr::mutate(id = paste(group, test_statistic, duration_d, endpoint, sep = ":"))

# setup to capture output
out.list <- NULL
for(i in 1:nrow(na.df)){
  # filter to the row we want
  na.df.row <- na.df[i,]
  
  # lets look at the raw data to see what is going on with these guys
  na.raw.df <- id.dat %>%
    dplyr::filter(id == na.df.row$x | id == na.df.row$y) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(x.n = ifelse(id == na.df.row$x, n(), NA_real_),
                  y.n = ifelse(id == na.df.row$y, n(), NA_real_))
  
  out.df <- tibble::tibble(x = na.df.row$x,
                           y = na.df.row$y,
                           x.n = max(na.raw.df$x.n, na.rm = T),
                           y.n = max(na.raw.df$y.n, na.rm = T))
    
  out.list[[i]] <- out.df
}

# bind the list
na.proc.df <- data.table::rbindlist(out.list) %>%
  dplyr::mutate(pair.n.min = ifelse(x.n >= y.n, y.n, x.n)) %>%
  dplyr::filter(pair.n.min > 1)

# explore the x,y pair with the largest min
na.proc.df.2 <- na.proc.df %>%
  dplyr::arrange(desc(pair.n.min)) %>%
  dplyr::slice(1)

# look at all values in raw dat, ok! This looks right, there are no shared drugs, so orthogonal regression is impossible!
na.proc.df.3 <- id.dat %>%
  dplyr::filter(id == na.proc.df.2$x | id == na.proc.df.2$y)

#==================================================#
# check 2: explore all the orth.reg r squared of 1
#==================================================#
# setup df to help understand why there are so many r squared of 1.
# these are all instances where there are only two datapoints - these should be filtered out.
set.seed(99)
rsq.df <- pw.orth.reg.df %>%
  dplyr::filter(orth.reg.r.squared == 1)


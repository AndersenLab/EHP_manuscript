#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/03_clean.csv")

# Widmayer vs Oncorhynchus mykiss 96 hr LC50
d1 <- dat %>%
  dplyr::filter((latin_name == "Oncorhynchus mykiss" &
                   duration_d == 4 &
                   test_statistic == "LC50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2)) %>%
  dplyr::arrange(desc(latin_name))

# run the orth regressions across all pairs
debug(pwOrthReg)
d1l <- pwOrthReg(data = d1, group = "latin_name", min.n = 10, limit.comp = "Caenorhabditis", message = T, plot = T)
d1df <- data.table::rbindlist(d1l$orthregs)
d1plots <- cowplot::plot_grid(plotlist = d1l$plots)
cowplot::ggsave2(d1plots, filename = "plots/Widmayer-Oncorhynchus_mykiss_96_hr_LC50.png", width = 12, height = 6)


# testing pair list
#assign("pairs.list", pairs.list, envir = .GlobalEnv)
pairs.list.test <- pairs.list
pairs.list.test[[7]] <- c(pairs.list.test[[3]][2], pairs.list.test[[3]][1])
limit.comp = "Caenorhabditis"
j <- 7
for(j in 1:length(pairs.list)){
  if((grepl(pairs.list[[j]][1], pattern = limit.comp) == F & grepl(pairs.list[[j]][2], pattern = limit.comp) == T))
  pairs.list[[j]] <- rev(pairs.list[[j]])
}
   
  
#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
#dat <- data.table::fread("data/processed/03_clean.csv")
dat <- data.table::fread("data/processed/00_data.csv")

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

# run the orth regressions across all pairs - NEED TO DEBUG
# to debug - debug(pwOrthReg) debug(orthReg)
d1l <- pwOrthReg(data = d1, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis elegans", message = T, plot = T)
d1df <- data.table::rbindlist(d1l$orthregs)
d1plots <- cowplot::plot_grid(plotlist = d1l$plots)
cowplot::ggsave2(d1plots, filename = "plots/test.png", width = 12, height = 6)

# Widmayer vs Pimephales promelas 96 hr LC50
d2 <- dat %>%
  dplyr::filter((latin_name == "Pimephales promelas" &
                   duration_d == 4 &
                   test_statistic == "LC50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2))

# run the orth regressions across all pairs
d2l <- pwOrthReg(data = d2, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d2df <- data.table::rbindlist(d2l$orthregs)
d2plots <- cowplot::plot_grid(plotlist = d2l$plots)
cowplot::ggsave2(d2plots, filename = "plots/Widmayer-Pimephales_promelas_96_hr_LC50.png", width = 12, height = 6)

# Widmayer vs Lepomis macrochirus 96 hr LC50
d3 <- dat %>%
  dplyr::filter((latin_name == "Lepomis macrochirus" &
                   duration_d == 4 &
                   test_statistic == "LC50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2))

# run the orth regressions across all pairs
d3l <- pwOrthReg(data = d3, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d3df <- data.table::rbindlist(d3l$orthregs)
d3plots <- cowplot::plot_grid(plotlist = d3l$plots)
cowplot::ggsave2(d3plots, filename = "plots/Widmayer-Lepomis_macrochirus_96_hr_LC50.png", width = 12, height = 6)

# Widmayer vs ALL FISH 96 hr LC50
d4 <- dat %>%
  dplyr::filter((group == "FISH" &
                   duration_d == 4 &
                   test_statistic == "LC50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2)) %>%
  dplyr::mutate(latin_name = case_when(group == "FISH" ~ "FISH",
                                       TRUE ~ latin_name)) # handle ALL FISH

# run the orth regressions across all pairs
d4l <- pwOrthReg(data = d4, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d4df <- data.table::rbindlist(d4l$orthregs)
d4plots <- cowplot::plot_grid(plotlist = d4l$plots)
cowplot::ggsave2(d4plots, filename = "plots/Widmayer-FISH_96_hr_LC50.png", width = 12, height = 6)

# Widmayer vs Daphnia magna 48 hr LC50 (it might be called EC50 because it is technically immobility not death)
d5 <- dat %>%
  dplyr::filter((latin_name == "Daphnia magna" &
                   duration_d == 4 &
                   test_statistic %in% c("LC50", "EC50")) |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2))

# run the orth regressions across all pairs
d5l <- pwOrthReg(data = d5, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d5df <- data.table::rbindlist(d5l$orthregs)
d5plots <- cowplot::plot_grid(plotlist = d5l$plots)
cowplot::ggsave2(d5plots, filename = "plots/Widmayer-Daphnia_magna_96_hr_LC50_EC50.png", width = 6, height = 6)

# Widmayer vs Americamysis bahia 96 hr LC50
d6 <- dat %>%
  dplyr::filter((latin_name == "Americamysis bahia" &
                   duration_d == 4 &
                   test_statistic == "LC50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2))

# run the orth regressions across all pairs
d6l <- pwOrthReg(data = d6, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
#d6df <- data.table::rbindlist(d6l$orthregs)
#d6plots <- cowplot::plot_grid(plotlist = d6l$plots)
#cowplot::ggsave2(d6plots, filename = "plots/Widmayer-Americamysis_bahia_96_hr_LC50.png", width = 12, height = 12)

# Widmayer vs ALL INVERT 48 hr LC50 (we may need to play with the duration for this one because each invert might have a different duration and some may be EC50 vs LC50 -- they should be a mortality effect)
d7 <- dat %>%
  dplyr::filter((group == "INVERT" &
                   duration_d >= 1 &
                   duration_d <= 4 &
                   endpoint == "Mortality") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2)) %>%
  dplyr::mutate(latin_name = case_when(group == "INVERT" ~ "INVERT",
                                       TRUE ~ latin_name)) # handle ALL INVERT


# run the orth regressions across all pairs
d7l <- pwOrthReg(data = d7, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d7df <- data.table::rbindlist(d7l$orthregs)
d7plots <- cowplot::plot_grid(plotlist = d7l$plots)
cowplot::ggsave2(d7plots, filename = "plots/Widmayer-INVERT_24-96_hr_Mortality.png", width = 24, height = 24)

# # Widmayer vs Lemna gibba (7 days or 168 hours) EC50 (growth) NOT IN DAT!
# d8 <- dat %>%
#   dplyr::filter((latin_name == "Lemna gibba" &
#                    duration_d == 7 &
#                    test_statistic == "EC50") |
#                   (latin_name == "Caenorhabditis elegans" &
#                      duration_d == 2))
# 
# # run the orth regressions across all pairs
# d8l <- pwOrthReg(data = d8, group = "latin_name", min.cases = 3, limit.comp = "Caenorhabditis", message = T, plot = T)
# d8df <- data.table::rbindlist(d8l$orthregs)
# d8plots <- cowplot::plot_grid(plotlist = d8l$plots)
# cowplot::ggsave2(d8plots, filename = "plots/Widmayer-Lemna_gibba_7_days_EC50.png", width = 12, height = 12)

# Widmayer vs ALGAE 96 hr LC50 (could also be a 72hr EC50 or IC50 for growth)
d9 <- dat %>%
  dplyr::filter((group == "ALGAE" &
                   duration_d >= 3 &
                   duration_d <= 4 &
                   test_statistic %in% c("LC50", "EC50", "IC50")) |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2)) %>%
  dplyr::mutate(latin_name = case_when(group == "ALGAE" ~ "ALGAE",
                                       TRUE ~ latin_name)) # handle ALL ALGAE


# run the orth regressions across all pairs
d9l <- pwOrthReg(data = d9, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d9df <- data.table::rbindlist(d9l$orthregs)
d9plots <- cowplot::plot_grid(plotlist = d9l$plots)
cowplot::ggsave2(d9plots, filename = "plots/Widmayer-ALGAE_3-4_days_LC50_EC50_IC50.png", width = 18, height = 18)

# Widmayer vs zebrafish embryo data
d10 <- dat %>%
  dplyr::filter((group %in% c("ZF_EMBRYO1", "ZF_EMBRYO2")) |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2)) %>%
  dplyr::mutate(latin_name = case_when(group == "ZF_EMBRYO1" ~ "ZF_EMBRYO1",
                                       group == "ZF_EMBRYO2" ~ "ZF_EMBRYO2",
                                       TRUE ~ latin_name)) # handle ZF_EMBRYO
# run the orth regressions across all pairs
d10l <- pwOrthReg(data = d10, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d10df <- data.table::rbindlist(d10l$orthregs)
d10plots <- cowplot::plot_grid(plotlist = d10l$plots)
cowplot::ggsave2(d10plots, filename = "plots/Widmayer-ZF_EMBRYO.png", width = 12, height = 12)

# Widmayer vs RAT
d11 <- dat %>%
  dplyr::filter((group %in% c("RAT_TEST", "RAT_EMPIRICAL")) |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2)) %>%
  dplyr::mutate(latin_name = case_when(group == "RAT_TEST" ~ "RAT_TEST",
                                       group == "RAT_EMPIRICAL" ~ "RAT_EMPIRICAL",
                                       TRUE ~ latin_name)) # handle RAT
# run the orth regressions across all pairs
d11l <- pwOrthReg(data = d11, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d11df <- data.table::rbindlist(d11l$orthregs)
d11plots <- cowplot::plot_grid(plotlist = d11l$plots)
cowplot::ggsave2(d11plots, filename = "plots/Widmayer-RAT.png", width = 12, height = 12)

#==================================================#
# Part 2: Consolidate priority orth_regs and sort
#==================================================#
# combine the data sets
proc_orthregs <- dplyr::bind_rows(d1df, d2df, d3df, d4df, d5df, d7df, d9df, d11df) %>% #d8df, #d6df,#d10df
  dplyr::filter(!is.na(orth.reg.slope)) 

# run the distance function - probably not the best function to find good comparisons, but cool attempt
proc_orthregs2 <- distFromIdeal(proc_orthregs)

# extract function outputs
dist.df <- proc_orthregs2$data
dist.plot <- proc_orthregs2$plot

# save these as an R object for now, so I can play with them
save(proc_orthregs, file = "data/processed/proc_orthregs.rda")
save(proc_orthregs2, file = "data/processed/proc_orthregs2.rda")

#================================================#
# Look at Boyd data vs imageXpress data
#================================================#
# Widmayer vs boyd vs Oncorhynchus mykiss 96 hr LC50
x1 <- dat %>%
  dplyr::filter((latin_name == "Oncorhynchus mykiss" &
                   duration_d == 4 &
                   test_statistic == "LC50") |
                  (latin_name == "Caenorhabditis elegans"))

# run the orth regressions across all pairs
x1l <- pwOrthReg(data = x1, group = "group", min.n = 5, limit.comp = "NEMATODE_COPAS1", message = T, plot = T)
x1df <- data.table::rbindlist(x1l$orthregs)
x1plots <- cowplot::plot_grid(plotlist = x1l$plots)
cowplot::ggsave2(x1plots, filename = glue::glue("plots/{today}_Widmayer-Boyd-Oncorhynchus_mykiss_96_hr_LC50_2.png"), width = 12, height = 12)

# Test the orthReg function again with single comp
x2 <- d1 %>%
  dplyr::mutate(keep = case_when((group == "NEMATODE" & test_statistic == "EC10") ~ "keep",
                                 (group == "FISH" &  test_statistic == "LC50" & duration_d == 4 & endpoint == "Mortality") ~ "keep",
                                 TRUE ~ "remove")) %>%
  dplyr::filter(keep == "keep")

test <- orthReg(x2, x = c("Caenorhabditis elegans", "EC10", "2", "Growth"), y = c("Oncorhynchus mykiss", "LC50", "4", "Mortality"), plot = T)
test$plot                  
  

#================================================#
# Look at all Algae endpoints
#================================================#
# Widmayer vs ALGAE 96 hr LC50 (could also be a 72hr EC50 or IC50 for growth)
alage1 <- dat %>%
  dplyr::filter((group == "ALGAE") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2)) %>%
  dplyr::mutate(latin_name = case_when(group == "ALGAE" ~ "ALGAE",
                                       TRUE ~ latin_name)) # handle ALL ALGAE


# run the orth regressions across all pairs
alage1l <- pwOrthReg(data = alage1, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
alage1df <- data.table::rbindlist(alage1l$orthregs) %>%
  dplyr::arrange(desc(orth.reg.n.observations), desc(orth.reg.r.squared))
alage1plots <- cowplot::plot_grid(plotlist = alage1l$plots)
cowplot::ggsave2(alage1plots, filename = "plots/Widmayer-ALGAE_all.png", width = 18, height = 18)

# save to send to scott
rio::export(alage1df, file = "data/processed/Widmayer_vs_ALGAE_orthReg.csv")

#================================================#
# Look at all group endpoints vs Widmayer
#================================================#
all.groups <- dat %>%
  dplyr::filter(group !="NEMATODE_COPAS1") %>%
  dplyr::mutate(latin_name = case_when(group != "NEMATODE" ~ group,
                                       group == "NEMATODE" ~ latin_name,
                                       TRUE ~ "wtf?")) # handle ALL

# run the orth regressions across all pairs
al <- pwOrthReg(data = all.groups, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
adf <- data.table::rbindlist(al$orthregs) %>%
  dplyr::arrange(desc(orth.reg.n.observations), desc(orth.reg.r.squared))
alplots <- cowplot::plot_grid(plotlist = al$plots)
cowplot::ggsave2(alplots, filename = "plots/Widmayer-allGROUPS.png", width = 49, height = 49)

# save to send to scott
rio::export(adf, file = "data/processed/Widmayer_vs_allGROUPS_orthReg.csv")

#================================================#
# Look at all species endpoints vs Widmayer
#================================================#
all.sp <- dat %>%
  dplyr::filter(group !="NEMATODE_COPAS1" & !(group == "INVERT" & latin_name == "Caenorhabditis elegans")) %>%
  dplyr::mutate(latin_name = paste0(group, "_", latin_name))

# run the orth regressions across all pairs
a.spl <- pwOrthReg(data = all.sp, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
a.spdf <- data.table::rbindlist(a.spl$orthregs) %>%
  dplyr::arrange(desc(orth.reg.n.observations), desc(orth.reg.r.squared))
a.spplots <- cowplot::plot_grid(plotlist = a.spl$plots)
cowplot::ggsave2(a.spplots, filename = "plots/Widmayer-allLATIN_NAME.png", width = 49, height = 49)

# save to send to scott
rio::export(a.spdf, file = "data/processed/Widmayer_vs_allLATIN_NAME_orthReg.csv")

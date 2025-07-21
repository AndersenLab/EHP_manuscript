#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/00_data.csv")

#===========================================================#
# Figure 4 - C. elegans toxicity data compared with data from
#  the "EnviroTox DB": fish 4-day mortality,
# invertebrate 1 and/or 2 day mortality, and algae mortality.
#===========================================================#
# 4a
# Filter to Boyd, Widmayer, and day 4 fish mortality from the EnviroTox DB
# Look at regression between widmayer and fish then Boyd and fish, plot both and update scott
# To merge fish endpoints - need LC50, which is mortality by definition
# Idea is to compare our data to Table 2, which was showing that Boyd is poor. Scott thinks should
# have our own approach to this.
df4a <- dat %>%
  dplyr::filter(source %in% c("Widmayer et al. 2022", "Boyd et al. 2016", "EnviroTox DB")) %>%
  dplyr::filter(group == "NEMATODE" | group == "FISH" | group == "NEMATODE_COPAS1") %>%
  dplyr::filter(case_when(group == "NEMATODE" ~ T,
                          group == "NEMATODE_COPAS1" ~ T,
                          source == "EnviroTox DB" & duration_d == 4 & test_statistic %in% c("LC50") ~ T, # removing EC50s has a big filtering effect, even though some are coded as mortality
                          TRUE ~ F)) %>%
  dplyr::mutate(endpoint2 = dplyr::case_when(endpoint %in% c("Mortality/Growth",
                                                           "Mortality, Mortality",
                                                           "Mortality, Survival") ~ "Mortality",
                                           endpoint == "Immobilization: Change in the failure to respond or lack of movement after mechanical stimulation." ~ "Immobility",
                                           endpoint == "Intoxication, Immobile" ~ "Immobility",
                                           TRUE ~ endpoint)) %>%
  dplyr::mutate(endpoint = endpoint2) %>%
  dplyr::select(-endpoint2)

# run the orth regressions across all pairs to NEMATODE with QC filter
or4a_QC <- pwOrthReg(data = df4a, group = "group",  limit.comp = "NEMATODE", min.n = 5, message = T, QC = "filter", plot = T)
or4adf_QC <- data.table::rbindlist(or4a_QC$orthregs)
print(glue::glue("slope = {round(or4adf_QC[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4adf_QC[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4adf_QC[1]$orth.reg.r.squared, digits = 2)}"))
or4ap_QC <- or4a_QC$plots[[1]]

# run the orth regressions across all pairs to NEMATODE without QC filter
or4a <- pwOrthReg(data = df4a, group = "group",  limit.comp = "NEMATODE", min.n = 5, QC = "ignore", message = T, plot = T)
or4adf <- data.table::rbindlist(or4a$orthregs)
print(glue::glue("slope = {round(or4adf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4adf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4adf[1]$orth.reg.r.squared, digits = 2)}"))
or4ap <- or4a$plots[[1]]
or4a$plots[[1]]

# Look at Boyd chems shared with Widmayer and compare to FISH
wb_overlap_fish <- df4a %>%
  dplyr::filter(chem_name %in% or4a$plots[[1]]$data$chem_name)
or4a2 <- pwOrthReg(data = wb_overlap_fish, group = "group",  limit.comp = "NEMATODE", min.n = 5, QC = "ignore", message = T, plot = T)
or4adf2 <- data.table::rbindlist(or4a2$orthregs)
or4a2$plots[[3]]
print(glue::glue("slope = {round(or4adf2[3]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4adf2[3]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4adf2[3]$orth.reg.r.squared, digits = 2)}"))

# need to add option to color QC fail data or filter QC fail data to orthReg function. - DONE, QC filter makes fit worse
# this will help visualize issues with the Boyd data
# consider adding a QC filter to enviroTox DB if endpoint changed? and flag that too?
# add notes to EHP manuscript DB in notion

#===========================================================#
# 4b
# Widmayer vs INVERTEBRATES 
# The Invert mortality data can sometimes be labelled as immobility or immobilization or survival select
# LC50 and EC50 could both be used b/c of immobilization endppoint which is essentially mortality
# plot comps like in 4b above and update scott
df4b <- dat %>%
  dplyr::filter(source %in% c("Widmayer et al. 2022", "Boyd et al. 2016", "EnviroTox DB")) %>%
    dplyr::filter(group == "NEMATODE" | group == "INVERT" | group == "NEMATODE_COPAS1") %>%
  dplyr::filter(case_when(group == "NEMATODE" ~ T,
                          group == "NEMATODE_COPAS1" ~ T,
                          source == "EnviroTox DB" & test_statistic %in% c("LC50", "EC50") ~ T, # removing EC50s has a big filtering effect, even though some are coded as mortality
                          TRUE ~ F)) %>%
  dplyr::mutate(endpoint2 = dplyr::case_when(endpoint %in% c("Mortality/Growth",
                                                             "Mortality, Mortality",
                                                             "Mortality, Survival") ~ "Mortality",
                                             endpoint == "Immobilization: Change in the failure to respond or lack of movement after mechanical stimulation." ~ "Immobility",
                                             endpoint == "Intoxication, Immobile" ~ "Immobility",
                                             TRUE ~ endpoint)) %>%
  dplyr::mutate(endpoint = endpoint2) %>%
  dplyr::select(-endpoint2) %>%
  dplyr::filter(case_when(source == "EnviroTox DB" & !(endpoint %in% c("Mortality", "Immobility"))  ~ F, # take either Mort or immobilization
                           TRUE ~ T)) %>%
  dplyr::mutate(endpoint = case_when(source == "EnviroTox DB" & endpoint %in% c("Mortality", "Immobility")  ~ "Mortality",
                                      TRUE ~ endpoint), # set both to mort
                 test_statistic = case_when(source == "EnviroTox DB"  ~ "LC50",
                                      TRUE ~ test_statistic)) # set test_stat to LC50 for both
  # Some improvement to merge EC50 Immobilization to LC50 Mort - 13 to 17 compounds max at LC50 Mort 2day, with imporved R^2

# run the orth regressions across all pairs without QC filter
or4b <- pwOrthReg(data = df4b, group = "group", limit.comp = "NEMATODE", min.n = 5, message = F, QC = "ignore", plot = T)
or4bdf <- data.table::rbindlist(or4b$orthregs)
or4bp <- or4b$plots[[3]]
# Widmayer vs. LC50 2d INVERT
or4b$plots[[1]]
print(glue::glue("slope = {round(or4bdf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4bdf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4bdf[1]$orth.reg.r.squared, digits = 2)}"))
# Widmayer vs. LC50 1d INVERT 
or4b$plots[[3]]
print(glue::glue("slope = {round(or4bdf[3]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4bdf[3]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4bdf[3]$orth.reg.r.squared, digits = 2)}"))
# Boyd vs. LC50 2d INVERT
or4b$plots[[9]]
print(glue::glue("slope = {round(or4bdf[9]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4bdf[9]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4bdf[9]$orth.reg.r.squared, digits = 2)}"))
# Boyd vs. LC50 1d INVERT 
or4b$plots[[11]]
print(glue::glue("slope = {round(or4bdf[11]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4bdf[11]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4bdf[11]$orth.reg.r.squared, digits = 2)}"))
# run the orth regressions across all pairs with QC filter, makes it worse for Boyd again.



#===========================================================#
# put them together
fig4 <- cowplot::plot_grid(or4ap + labs(subtitle = ""),
                           or4bp + labs(subtitle = ""),
                           or4cp + labs(subtitle = ""),
                           labels = c("A", "B", "C"), align = "vh", ncol = 1)
cowplot::ggsave2(fig4, filename = glue::glue("figures/{today}_figure4.png"), width = 3, height = 9)

#===========================================================#
# Figure 4 - C. elegans toxicity data compared with common
#            surrogate species, D. magna, C. dubia, O. mykiss,
#            D. rerio, P. promelas, and L. macrochirus
#===========================================================#
# 4a
# Widmayer vs daphnia manga
df4a <- dat %>%
  dplyr::filter(group == "NEMATODE" | grepl(latin_name, pattern = " magna"))
# 4b
# Widmayer vs daphnia manga
df4b <- dat %>%
  dplyr::filter(group == "NEMATODE" | grepl(latin_name, pattern = " dubia"))
# 4c
# Widmayer vs daphnia manga
df4c <- dat %>%
  dplyr::filter(group == "NEMATODE" | grepl(latin_name, pattern = " mykiss"))
# 4d
# Widmayer vs daphnia manga
df4d <- dat %>%
  dplyr::filter(group == "NEMATODE" | grepl(latin_name, pattern = " rerio"))
# 4e
# Widmayer vs daphnia manga
df4e <- dat %>%
  dplyr::filter(group == "NEMATODE" | grepl(latin_name, pattern = " promelas"))
# 4f
# Widmayer vs daphnia manga
df4f <- dat %>%
  dplyr::filter(group == "NEMATODE" | grepl(latin_name, pattern = " macrochirus"))

#===========================================================#
# run the orth regressions across all pairs
or4a <- pwOrthReg(data = df4a, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.n = 3, message = F, plot = T)
or4adf <- data.table::rbindlist(or4a$orthregs)
or4ap <- or4a$plots[[1]] # this is std.
print(glue::glue("slope = {round(or4adf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4adf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4adf[1]$orth.reg.r.squared, digits = 2)}"))

# run the orth regressions across all pairs
or4b <- pwOrthReg(data = df4b, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.n = 3, message = F, plot = T)
or4bdf <- data.table::rbindlist(or4b$orthregs)
or4bp <- or4b$plots[[2]] # std.
print(glue::glue("slope = {round(or4bdf[2]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4bdf[2]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4bdf[2]$orth.reg.r.squared, digits = 2)}"))

# run the orth regressions across all pairs
or4c <- pwOrthReg(data = df4c, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.n = 3, message = F, plot = T)
or4cdf <- data.table::rbindlist(or4c$orthregs)
or4cp <- or4c$plots[[1]]
print(glue::glue("slope = {round(or4cdf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4cdf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4cdf[1]$orth.reg.r.squared, digits = 2)}"))

# run the orth regressions across all pairs
or4d <- pwOrthReg(data = df4d, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.n = 3, message = F, plot = T)
or4ddf <- data.table::rbindlist(or4d$orthregs)
or4dp <- or4d$plots[[3]] # std
print(glue::glue("slope = {round(or4ddf[3]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4ddf[3]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4ddf[3]$orth.reg.r.squared, digits = 2)}"))

# run the orth regressions across all pairs
or4e <- pwOrthReg(data = df4e, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.n = 3, message = F, plot = T)
or4edf <- data.table::rbindlist(or4e$orthregs)
or4ep <- or4e$plots[[1]]
print(glue::glue("slope = {round(or4edf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4edf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4edf[1]$orth.reg.r.squared, digits = 2)}"))

# run the orth regressions across all pairs
or4f <- pwOrthReg(data = df4f, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.n = 3, message = F, plot = T)
or4fdf <- data.table::rbindlist(or4f$orthregs)
or4fp <- or4f$plots[[1]]
print(glue::glue("slope = {round(or4fdf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4fdf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4fdf[1]$orth.reg.r.squared, digits = 2)}"))

#===========================================================#
# put them together
fig4 <- cowplot::plot_grid(or4ap, or4bp, or4cp, or4dp, or4ep, or4fp, labels = c("A", "B", "C", "D", "E", "F"), align = "vh", nrow = 3)

cowplot::ggsave2(fig4, filename = "figures/figure4.png", width = 10, height = 15)


# old algae
# temp export for scott
#===========================================================#
# 4d POSSIBLY CUT 
# Widmayer vs ALGAE 1 d, 2d, 3d, 4d possible and there are not clear guidelines can test all?
# Try to merge obvious mortality categories (population?)
# Same issue with EC50 and LC50, and IC50
# update scott
df4c_test <- dat %>%
  dplyr::filter(source %in% c("Widmayer et al. 2022", "Boyd et al. 2016", "EnviroTox DB")) %>%
  dplyr::filter(group == "NEMATODE" | group == "ALGAE" | group == "NEMATODE_COPAS1") %>%
  dplyr::filter(endpoint == "Population")
table(df4c_test$test_statistic) # just pull the EC50 population, because it is by far the biggest group

# filter it
df4c <- dat %>%
  dplyr::filter(source %in% c("Widmayer et al. 2022", "Boyd et al. 2016", "EnviroTox DB")) %>%
  dplyr::filter(group == "NEMATODE" | group == "ALGAE" | group == "NEMATODE_COPAS1") %>%
  dplyr::filter(case_when(group == "NEMATODE" ~ T,
                          group == "NEMATODE_COPAS1" ~ T,
                          source == "EnviroTox DB" & test_statistic == "EC50" & endpoint == "Population"  ~ T, 
                          TRUE ~ F))

# run the orth regressions across all pairs
or4c <- pwOrthReg(data = df4c, group = "group", limit.comp = "NEMATODE", min.n = 3, message = F, QC = "ignore", plot = T)
or4cdf <- data.table::rbindlist(or4c$orthregs)
or4cp <- or4c$plots[[1]]
#Widmayer vs. 1 day EC50 population
or4c$plots[[1]]
print(glue::glue("slope = {round(or4cdf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4cdf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4cdf[1]$orth.reg.r.squared, digits = 2)}"))

# Boyd vs. 1 day EC50 population
or4c$plots[[11]]
print(glue::glue("slope = {round(or4cdf[11]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4cdf[11]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4cdf[11]$orth.reg.r.squared, digits = 2)}"))

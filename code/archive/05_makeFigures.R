#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/00_data.csv")
#===========================================================#
# Figure 4 - C. elegans toxicity data compared with rats,
#            fish, invertebrates, and algae.
#===========================================================#
# 3a
# Widmayer vs RAT
df3a <- dat %>%
  dplyr::filter(group == "NEMATODE" | group == "RAT_EMPIRICAL")

# run the orth regressions across all pairs
or3a <- pwOrthReg(data = df3a, group = "group", min.n = 5, message = T, plot = T)
or3adf <- data.table::rbindlist(or3a$orthregs)
print(glue::glue("slope = {round(or3adf[3]$orth.reg.slope, digits = 2)}, y-intercept = {round(or3adf[3]$orth.reg.intercept, digits = 2)}, r^2 = {round(or3adf[3]$orth.reg.r.squared, digits = 2)}"))
or3ap <- or3a$plots[[3]]

#===========================================================#
# 3b
# Widmayer vs FISH
df3b <- dat %>%
  dplyr::filter(group == "NEMATODE" | group == "FISH")

# run the orth regressions across all pairs
or3b <- pwOrthReg(data = df3b, group = "group", limit.comp = "NEMATODE", min.n = 5, message = F, plot = T)
or3bdf <- data.table::rbindlist(or3b$orthregs)
print(glue::glue("slope = {round(or3bdf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or3bdf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or3bdf[1]$orth.reg.r.squared, digits = 2)}"))
or3bp <- or3b$plots[[1]]
or3bp2 <- or3b$plots[[8]]

#===========================================================#
# 3c
# Widmayer vs INVERTEBRATES
df3c <- dat %>%
  dplyr::filter(group == "NEMATODE" | group == "INVERT")

# run the orth regressions across all pairs
or3c <- pwOrthReg(data = df3c, group = "group", limit.comp = "NEMATODE", min.n = 5, message = F, plot = T)
or3cdf <- data.table::rbindlist(or3c$orthregs)
or3cp <- or3c$plots[[2]]
print(glue::glue("slope = {round(or3cdf[2]$orth.reg.slope, digits = 2)}, y-intercept = {round(or3cdf[2]$orth.reg.intercept, digits = 2)}, r^2 = {round(or3cdf[2]$orth.reg.r.squared, digits = 2)}"))

# temp export for scott
#cowplot::ggsave2(or3cp, filename = "plots/Nematode-EC10-2d-Growth_Invert-LC50-1d-Mortality.png", width = 6, height = 6)
#cowplot::ggsave2(or3bp, filename = "plots/Nematode-EC10-2d-Growth_Fish-LC50-4d-Mortality.png", width = 6, height = 6)
#===========================================================#
# 3d
# Widmayer vs ALGAE
df3d <- dat %>%
  dplyr::filter(group == "NEMATODE" | group == "ALGAE")

# run the orth regressions across all pairs
or3d <- pwOrthReg(data = df3d, group = "group", limit.comp = "NEMATODE", min.n = 3, message = F, plot = T)
or3ddf <- data.table::rbindlist(or3d$orthregs)
or3dp <- or3d$plots[[3]]
print(glue::glue("slope = {round(or3ddf[3]$orth.reg.slope, digits = 2)}, y-intercept = {round(or3ddf[3]$orth.reg.intercept, digits = 2)}, r^2 = {round(or3ddf[3]$orth.reg.r.squared, digits = 2)}"))
or3dp2 <- or3d$plots[[1]]
or3dp3 <- or3d$plots[[4]]

#===========================================================#
# put them together
figure3 <- cowplot::plot_grid(or3ap, or3bp, or3cp, or3dp, labels = c("A", "B", "C", "D"), align = "vh")
cowplot::ggsave2(figure3, filename = "figures/figure3.png", width = 10, height = 10)

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

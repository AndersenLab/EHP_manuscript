#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/03_clean.csv")

#===========================================================#
# Figure 2 - Boyd COPAS data vs. ImageXpress
#===========================================================#
# # Widmayer vs boyd
# df2 <- dat %>%
#   dplyr::filter(latin_name == "Caenorhabditis elegans")
# 
# # run the orth regressions across all pairs
# or2 <- pwOrthReg(data = df2, group = "latin_name", min.n = 3, message = T, plot = T)
# or2df <- data.table::rbindlist(or2$orthregs)
# or2p <- or2$plots[[4]]
# cowplot::ggsave2(or2p, filename = "figures/figure2.png", width = 6, height = 6)

#==============================================================================#
# step 1: pull Widmayer vs boyd EC10 vs AC50 w/ Carbryl and shape for plotting
#==============================================================================#
# Widmayer vs boyd EC10 vs AC50, leave in Carbaryl for plotting
df2 <- dat %>%
  dplyr::filter(latin_name == "Caenorhabditis elegans" &
                  test_statistic %in% c("EC10", "AC50"))

# set args
# x is a vector of specific latin_name/group, test_statistic, duration_d, endpoint, so is y
x = c("Caenorhabditis elegans", "EC10", "2", "Growth")
y = c("Caenorhabditis elegans", "AC50", "NA", "Growth")

or2_proc <- df2 %>%
  dplyr::mutate(duration_d = ifelse(is.na(duration_d), "NA", duration_d)) %>% # THIS STEP HANDELS NAs in duration data
  dplyr::mutate(endpoint = ifelse(is.na(endpoint), "NA", endpoint)) %>% # THIS STEP HANDELS NAs in endpoint data???????
  dplyr::mutate(pair = dplyr::case_when((latin_name == x[1] | group == x[1]) & test_statistic == x[2] & duration_d == x[3] & endpoint == x[4] ~ x[1],
                                        (latin_name == y[1] | group == y[1]) & test_statistic == y[2] & duration_d == y[3] & endpoint == y[4] ~ y[1],
                                        TRUE ~ NA_character_),
                pair_gen = dplyr::case_when((latin_name == x[1] | group == x[1]) & test_statistic == x[2] & duration_d == x[3] & endpoint == x[4] ~ "x",
                                            (latin_name == y[1] | group == y[1]) & test_statistic == y[2] & duration_d == y[3] & endpoint == y[4] ~ "y",
                                            TRUE ~ NA_character_)) %>% # label pairs, should handle giving a group or a latin_name since they are unique
  dplyr::filter(!is.na(pair_gen)) %>% # filter to pairs with labels NEW pair_gen OLD pair
  dplyr::group_by(cas, pair_gen) %>% # NEW pair_gen OLD pair
  dplyr::mutate(gm_mean = gm_mean(effect_value),
                min = min(effect_value),
                max = max(effect_value)) %>% # get geometric mean for chemical and pair
  dplyr::ungroup()

# reshap for plotting
plot_dat <- or2_proc %>%
  dplyr::distinct(chem_name, pair_gen, gm_mean, min, max) %>% # just get geom mean and ranges
  tidyr::pivot_wider(names_from = pair_gen, values_from = c(gm_mean, min, max)) %>% # give us x and a y vars
  dplyr::filter(complete.cases(.)) %>% # keep only complete cases
  dplyr::mutate(fill = ifelse(chem_name == "Carbaryl", "grey", "red"))
#==============================================================================#
# Step 2: take out Carbaryl to run the analysis
#==============================================================================#
df2_noCarb <- dat %>%
  dplyr::filter(latin_name == "Caenorhabditis elegans" &
                  test_statistic %in% c("EC10", "AC50") &
                  chem_name != "Carbaryl") %>%
  dplyr::mutate(duration_d = ifelse(is.na(test_statistic), ""))

# shape data and calc geometric mean
or2_noCarb_proc <- df2_noCarb %>%
  dplyr::mutate(duration_d = ifelse(is.na(duration_d), "NA", duration_d)) %>% # THIS STEP HANDELS NAs in duration data
  dplyr::mutate(endpoint = ifelse(is.na(endpoint), "NA", endpoint)) %>% # THIS STEP HANDELS NAs in endpoint data???????
  dplyr::mutate(pair = dplyr::case_when((latin_name == x[1] | group == x[1]) & test_statistic == x[2] & duration_d == x[3] & endpoint == x[4] ~ x[1],
                                        (latin_name == y[1] | group == y[1]) & test_statistic == y[2] & duration_d == y[3] & endpoint == y[4] ~ y[1],
                                        TRUE ~ NA_character_),
                pair_gen = dplyr::case_when((latin_name == x[1] | group == x[1]) & test_statistic == x[2] & duration_d == x[3] & endpoint == x[4] ~ "x",
                                            (latin_name == y[1] | group == y[1]) & test_statistic == y[2] & duration_d == y[3] & endpoint == y[4] ~ "y",
                                            TRUE ~ NA_character_)) %>% # label pairs, should handle giving a group or a latin_name since they are unique
  dplyr::filter(!is.na(pair_gen)) %>% # filter to pairs with labels NEW pair_gen OLD pair
  dplyr::group_by(cas, pair_gen) %>% # NEW pair_gen OLD pair
  dplyr::mutate(gm_mean = gm_mean(effect_value),
                min = min(effect_value),
                max = max(effect_value)) %>% # get geometric mean for chemical and pair
  dplyr::ungroup()

# reshape for plotting and orthogonal regression
noCarb_plot_dat <- or2_noCarb_proc %>%
  dplyr::distinct(chem_name, pair_gen, gm_mean, min, max) %>% # just get geom mean and ranges
  tidyr::pivot_wider(names_from = pair_gen, values_from = c(gm_mean, min, max)) %>% # give us x and a y vars
  dplyr::filter(complete.cases(.)) # keep only complete cases

# ORIGINAL REGRESSION BASED ON gm_mean
orthog_reg_model_log10 <- pracma::odregress(x = log10(noCarb_plot_dat$gm_mean_x), y = log10(noCarb_plot_dat$gm_mean_y))

# # old orthogonal regression
oldOR_df <- tibble::tibble(x = log10(noCarb_plot_dat$gm_mean_x), y = log10(noCarb_plot_dat$gm_mean_y))
pcObject <- princomp(oldOR_df)
myLoadings <- unclass(loadings(pcObject))[,1]
OR.slope <- myLoadings["y"]/myLoadings["x"]
OR.int <- mean(log10(noCarb_plot_dat$gm_mean_y))-OR.slope*mean(log10(noCarb_plot_dat$gm_mean_x))
###Rsq is empirical in this setting.  The first principal component will always explain at least half
###of the total variation, since the first two components will explain 100% of it in this simple setting,
###and the first must explain at least as much as the second.
orthog_reg_model_r_squared <- as.double(((summary(pcObject)$sdev[1]^2)/sum(summary(pcObject)$sdev^2)-.5)/.5)
# corrCoef <- cor(log10(plot_dat$gm_mean_x),log10(plot_dat$gm_mean_y))
# corr.rsq <- corrCoef^2

# run the extraction functions
orthog_reg_model_mse <- mse_odreg(orthog_reg_model_log10)
# NEW Rsq!!!!! DIFFERNET THAN OLD - using old for now
#orthog_reg_model_r_squared <- r_squared_odreg(orthog_reg_model_log10, log10(plot_dat$gm_mean_y))

# build output for orthogonal regression
noCarb_orthreg_df <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                             y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                             orth.reg.n.observations = nrow(noCarb_plot_dat),
                             orth.reg.slope = orthog_reg_model_log10$coeff[1],
                             orth.reg.intercept = orthog_reg_model_log10$coeff[2],
                             orth.reg.ssq = orthog_reg_model_log10$ssq[1],
                             orth.reg.mse = orthog_reg_model_mse,
                             orth.reg.r.squared = orthog_reg_model_r_squared)

# plot all data with orthogonal regression from Carbaryl removed
fig2 <- ggplot2::ggplot(plot_dat) +
  ggplot2::aes(x = gm_mean_x, y = gm_mean_y) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
  ggplot2::geom_abline(slope = orthog_reg_model_log10$coeff[1], intercept = orthog_reg_model_log10$coeff[2], size = 0.5, color = "red") +
  ggplot2::geom_errorbar(aes(ymin = min_y, ymax = max_y), width = 0, size = 0.25, color = "grey70") +
  ggplot2::geom_errorbarh(aes(xmin = min_x, xmax = max_x), height = 0, size = 0.25, color = "grey70") +
  ggplot2::geom_point(aes(fill = fill), show.legend = F, shape = 21, color = "black") +
  ggplot2::scale_fill_manual(values = c("grey" = "grey", "red" = "red")) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = bquote(~italic(.(x[1]))~"toxicity"~.(x[2])~"(μg/L)"),
                y = bquote(~italic(.(y[1]))~"toxicity"~.(y[2])~"(μg/L)")) +
                #subtitle = glue::glue("x={x[1]}_{x[2]}_{x[3]}d_{x[4]}\ny={y[1]}_{y[2]}_{y[3]}d_{y[4]}")) +
  ggplot2::scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(0.01, 300)
  ) +
  ggplot2::scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(0.01, 300)
  ) +
  ggplot2::annotation_logticks() 
fig2

# save the figure
cowplot::ggsave2(fig2, filename = "figures/figure2.png", width = 6, height = 6)

# save the data for ploting and the regression coefficients
rio::export(plot_dat, file = "data/processed/fig2.plot.data.csv")
rio::export(noCarb_orthreg_df, file = "data/processed/fig2.orth.reg.csv")
#===========================================================#
# Figure 3 - C. elegans toxicity data compared with rats,
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
or3d <- pwOrthReg(data = df3d, group = "group", limit.comp = "NEMATODE", min.cases = 3, message = F, plot = T)
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
or4a <- pwOrthReg(data = df4a, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.cases = 3, message = F, plot = T)
or4adf <- data.table::rbindlist(or4a$orthregs)
or4ap <- or4a$plots[[1]] # this is std.
print(glue::glue("slope = {round(or4adf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4adf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4adf[1]$orth.reg.r.squared, digits = 2)}"))

# run the orth regressions across all pairs
or4b <- pwOrthReg(data = df4b, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.cases = 3, message = F, plot = T)
or4bdf <- data.table::rbindlist(or4b$orthregs)
or4bp <- or4b$plots[[2]] # std.
print(glue::glue("slope = {round(or4bdf[2]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4bdf[2]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4bdf[2]$orth.reg.r.squared, digits = 2)}"))

# run the orth regressions across all pairs
or4c <- pwOrthReg(data = df4c, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.cases = 3, message = F, plot = T)
or4cdf <- data.table::rbindlist(or4c$orthregs)
or4cp <- or4c$plots[[1]]
print(glue::glue("slope = {round(or4cdf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4cdf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4cdf[1]$orth.reg.r.squared, digits = 2)}"))

# run the orth regressions across all pairs
or4d <- pwOrthReg(data = df4d, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.cases = 3, message = F, plot = T)
or4ddf <- data.table::rbindlist(or4d$orthregs)
or4dp <- or4d$plots[[3]] # std
print(glue::glue("slope = {round(or4ddf[3]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4ddf[3]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4ddf[3]$orth.reg.r.squared, digits = 2)}"))

# run the orth regressions across all pairs
or4e <- pwOrthReg(data = df4e, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.cases = 3, message = F, plot = T)
or4edf <- data.table::rbindlist(or4e$orthregs)
or4ep <- or4e$plots[[1]]
print(glue::glue("slope = {round(or4edf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4edf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4edf[1]$orth.reg.r.squared, digits = 2)}"))

# run the orth regressions across all pairs
or4f <- pwOrthReg(data = df4f, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.cases = 3, message = F, plot = T)
or4fdf <- data.table::rbindlist(or4f$orthregs)
or4fp <- or4f$plots[[1]]
print(glue::glue("slope = {round(or4fdf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4fdf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4fdf[1]$orth.reg.r.squared, digits = 2)}"))

#===========================================================#
# put them together
fig4 <- cowplot::plot_grid(or4ap, or4bp, or4cp, or4dp, or4ep, or4fp, labels = c("A", "B", "C", "D", "E", "F"), align = "vh", nrow = 3)

cowplot::ggsave2(fig4, filename = "figures/figure4.png", width = 10, height = 15)

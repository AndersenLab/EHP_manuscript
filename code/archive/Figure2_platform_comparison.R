#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/00_data.csv")
dat2 <- data.table::fread("data/processed/00_data2.csv")
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
df2 <- dat2 %>%
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
                  chem_name != "Carbaryl")

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

# just make a lm to test
lmfig2 <- lm(data = noCarb_plot_dat, gm_mean_y ~ gm_mean_x)
summary(lmfig2)

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
  ggplot2::labs(x = bquote(~"Nematode Imager"~.(x[2])~"(mg/L)"),
                y = bquote(~"Nematode COPAS"~.(y[2])~"(mg/L)")) +
  #ggplot2::labs(x = bquote(~italic(.(x[1]))~"toxicity"~.(x[2])~"(mg/L)"),
  #              y = bquote(~italic(.(y[1]))~"toxicity"~.(y[2])~"(mg/L)")) +
  #subtitle = glue::glue("x={x[1]}_{x[2]}_{x[3]}d_{x[4]}\ny={y[1]}_{y[2]}_{y[3]}d_{y[4]}")) +
  ggplot2::scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(0.00075, 500)
  ) +
  ggplot2::scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(0.00075, 500)
  ) +
  annotate(geom = "text", x = 30, y = 0.02, label = glue::glue("carbaryl"), size = 4) +
  #annotate(geom = "text", x = 3, y = 0.02, label = glue::glue("r = {corr.all}"), size = 3) +
  ggplot2::annotation_logticks() 
fig2

# save the figure
cowplot::ggsave2(fig2, filename = "figures/figure2_v2.png", width = 4.5, height = 4.5)

# save the data for ploting and the regression coefficients
rio::export(plot_dat, file = "data/processed/fig2.plot.data.csv")
rio::export(noCarb_orthreg_df, file = "data/processed/fig2.orth.reg.csv")

#==============================================================================#
# Step 3: Take out Carbaryl and the other QC FAILED compounds to run the analysis
#==============================================================================#
# look at QC status for Boyd
df2_boyd_QC_fail <- or2_proc %>%
  dplyr::filter(QC == "FAIL") %>%
  dplyr::filter(chem_name %in% plot_dat$chem_name) %>%
  dplyr::pull(chem_name)

plot_dat_QC_fail <- or2_proc %>%
  dplyr::distinct(chem_name, pair_gen, gm_mean, min, max) %>% # just get geom mean and ranges
  tidyr::pivot_wider(names_from = pair_gen, values_from = c(gm_mean, min, max)) %>% # give us x and a y vars
  dplyr::filter(complete.cases(.)) %>% # keep only complete cases
  dplyr::mutate(fill = ifelse(chem_name == "Carbaryl", "grey", 
                              ifelse(chem_name %in% df2_boyd_QC_fail, "black", "red")))

# shape the data
df2_QC_clean <- dat %>%
  dplyr::filter(latin_name == "Caenorhabditis elegans" &
                  test_statistic %in% c("EC10", "AC50") &
                  !(chem_name %in% c("Carbaryl", df2_boyd_QC_fail)))

# shape data and calc geometric mean
or2_QC_clean <- df2_QC_clean %>%
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
QC_clean_plot_dat <- or2_QC_clean %>%
  dplyr::distinct(chem_name, pair_gen, gm_mean, min, max) %>% # just get geom mean and ranges
  tidyr::pivot_wider(names_from = pair_gen, values_from = c(gm_mean, min, max)) %>% # give us x and a y vars
  dplyr::filter(complete.cases(.)) # keep only complete cases

# ORIGINAL REGRESSION BASED ON gm_mean
orthog_reg_model_log10_QC_clean <- pracma::odregress(x = log10(QC_clean_plot_dat$gm_mean_x), y = log10(QC_clean_plot_dat$gm_mean_y))

# # old orthogonal regression
oldOR_df_QC_clean <- tibble::tibble(x = log10(QC_clean_plot_dat$gm_mean_x), y = log10(QC_clean_plot_dat$gm_mean_y))
pcObject_QC_clean <- princomp(oldOR_df_QC_clean)
myLoadings_QC_clean <- unclass(loadings(pcObject_QC_clean))[,1]
OR.slope_QC_clean <- myLoadings_QC_clean["y"]/myLoadings_QC_clean["x"]
OR.int_QC_clean <- mean(log10(QC_clean_plot_dat$gm_mean_y))-OR.slope*mean(log10(QC_clean_plot_dat$gm_mean_x))
###Rsq is empirical in this setting.  The first principal component will always explain at least half
###of the total variation, since the first two components will explain 100% of it in this simple setting,
###and the first must explain at least as much as the second.
orthog_reg_model_r_squared_QC_clean <- as.double(((summary(pcObject_QC_clean)$sdev[1]^2)/sum(summary(pcObject_QC_clean)$sdev^2)-.5)/.5)
# corrCoef <- cor(log10(plot_dat$gm_mean_x),log10(plot_dat$gm_mean_y))
# corr.rsq <- corrCoef^2

# run the extraction functions
orthog_reg_model_mse_QC_clean <- mse_odreg(orthog_reg_model_log10_QC_clean)
# NEW Rsq!!!!! DIFFERNET THAN OLD - using old for now
#orthog_reg_model_r_squared <- r_squared_odreg(orthog_reg_model_log10, log10(plot_dat$gm_mean_y))


# build output for orthogonal regression
QC_clean_orthreg_dfthreg_df <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                                              y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                                              orth.reg.n.observations = nrow(QC_clean_plot_dat),
                                              orth.reg.slope = orthog_reg_model_log10_QC_clean$coeff[1],
                                              orth.reg.intercept = orthog_reg_model_log10_QC_clean$coeff[2],
                                              orth.reg.ssq = orthog_reg_model_log10_QC_clean$ssq[1],
                                              orth.reg.mse = orthog_reg_model_mse_QC_clean,
                                              orth.reg.r.squared = orthog_reg_model_r_squared_QC_clean)

# plot all data with orthogonal regression from Carbaryl removed
fig2_QC_clean <- ggplot2::ggplot(plot_dat_QC_fail) +
  ggplot2::aes(x = gm_mean_x, y = gm_mean_y) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
  ggplot2::geom_abline(slope = orthog_reg_model_log10_QC_clean$coeff[1], intercept = orthog_reg_model_log10_QC_clean$coeff[2], size = 0.5, color = "red") +
  ggplot2::geom_errorbar(aes(ymin = min_y, ymax = max_y), width = 0, size = 0.25, color = "grey70") +
  ggplot2::geom_errorbarh(aes(xmin = min_x, xmax = max_x), height = 0, size = 0.25, color = "grey70") +
  ggplot2::geom_point(aes(fill = fill), show.legend = F, shape = 21, color = "black") +
  ggplot2::scale_fill_manual(values = c("grey" = "grey", "black" = "black",  "red" = "red")) +
  ggplot2::theme_bw() +
  #ggplot2::labs(x = bquote(~italic(.(x[1]))~"toxicity"~.(x[2])~"(mg/L)"),
  #              y = bquote(~italic(.(y[1]))~"toxicity"~.(y[2])~"(mg/L)")) +
  ggplot2::labs(x = bquote(~"Nematode Imager"~"toxicity"~.(x[2])~"(mg/L)"),
                y = bquote(~"Nematode COPAS"~"toxicity"~.(y[2])~"(mg/L)")) +
  #subtitle = glue::glue("x={x[1]}_{x[2]}_{x[3]}d_{x[4]}\ny={y[1]}_{y[2]}_{y[3]}d_{y[4]}")) +
  ggplot2::scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(0.00075, 500)
  ) +
  ggplot2::scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(0.00075, 500)
  ) +
  annotate(geom = "text", x = 30, y = 0.02, label = glue::glue("carbaryl"), size = 4) +
  ggplot2::annotation_logticks() 
fig2_QC_clean

cowplot::ggsave2(fig2_QC_clean, file = glue::glue("plots/{today}_fig2_QC_clean.png"), width = 4.5, height = 4.5)
# The QC failed compounds do not appear to be outliers - Could remove this section if we stick with QC failed compounds
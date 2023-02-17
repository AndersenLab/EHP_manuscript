#!/usr/bin/env Rscript
library(tidyverse)

#==================================================#
# Geometric mean function from P. McMurdie 
#==================================================#
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

#======================================================================#
# function for data selection, shaping, orthogonal regression, plotting 
#======================================================================#
orthReg <- function(data, x, y){
  # make a list to hold it all
  out <- NULL
  
  # filter data to focal taxa and and geometric mean for group
  all_focal_dat <- data %>%
    dplyr::mutate(pair = dplyr::case_when((latin_name == x[1] | group == x[1]) & test_statistic == x[2] & duration_h == x[3] ~ x[1],
                                          (latin_name == y[1] | group == y[1]) & test_statistic == y[2] & duration_h == y[3] ~ y[1],
                                          TRUE ~ NA_character_),
                  pair_gen = dplyr::case_when((latin_name == x[1] | group == x[1]) & test_statistic == x[2] & duration_h == x[3] ~ "x",
                                              (latin_name == y[1] | group == y[1]) & test_statistic == y[2] & duration_h == y[3] ~ "y",
                                              TRUE ~ NA_character_)) %>% # label pairs, should handle giving a group or a latin_name since they are unique
    dplyr::filter(!is.na(pair)) %>% # filter to pairs with labels
    dplyr::group_by(cas, pair) %>% 
    dplyr::mutate(gm_mean = gm_mean(effect_value),
                  min = min(effect_value),
                  max = max(effect_value)) %>% # get geometric mean for chemical and pair
    dplyr::ungroup()
  
  # reshape for plotting
  plot_dat <- all_focal_dat %>%
    dplyr::distinct(chem_name, pair_gen, gm_mean, min, max) %>% # just get geom mean and ranges
    tidyr::pivot_wider(names_from = pair_gen, values_from = c(gm_mean, min, max)) %>% # give us x and a y vars
    dplyr::filter(complete.cases(.)) # keep only complete cases
  
  # get the orthogonal regression values with odregress function - use log transform of x and y
  # coeff[1] = slope and coeff[2] = intercept
  orthog_reg_model_log10 <- pracma::odregress(x = log10(plot_dat$gm_mean_x), y = log10(plot_dat$gm_mean_y))
  
  # get orthogonal regression coeffs with old code
  logMeanData <- data.frame(x = log10(plot_dat$gm_mean_x), y = log10(plot_dat$gm_mean_y), chem = plot_dat$chem_name)
  indata <- logMeanData[,c("x","y")]
  #print(summary(pcObject <- princomp(indata)))
  print(summary(pcObject <- princomp(indata)))
  #myLoadings <- unclass(loadings(pcObject))[,1]
  myLoadings <- unclass(loadings(pcObject))[,1]
  #OR.slope <- myLoadings["y"]/myLoadings["x"]
  OR.slope <- myLoadings["y"]/myLoadings["x"]
  ### oreg line goes through the center of the data, so...
  #OR.int <- mean(indata$y)-OR.slope*mean(indata$x)
  OR.int <- mean(indata$y)-OR.slope*mean(indata$x)
  ###Rsq is empirical in this setting.  The first principal component will always explain at least half
  ###of the total variation, since the first two components will explain 100% of it in this simple setting,
  ###and the first must explain at least as much as the second.
  princomp.rsq <- ((summary(pcObject)$sdev[1]^2)/sum(summary(pcObject)$sdev^2)-.5)/.5
  old_OR_dat <- tibble::tibble(OR.slope = OR.slope,
                               OR.int = OR.int,
                               princomp.rsq = princomp.rsq)
  
  # plot it with log ticks
  plot <- ggplot2::ggplot(plot_dat) +
    ggplot2::aes(x = gm_mean_x, y = gm_mean_y) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
    ggplot2::geom_abline(slope = orthog_reg_model_log10$coeff[1], intercept = orthog_reg_model_log10$coeff[2], size = 0.5, color = "red") +
    ggplot2::geom_errorbar(aes(ymin = min_y, ymax = max_y), size = 0.25, color = "grey70") +
    ggplot2::geom_errorbarh(aes(xmin = min_x, xmax = max_x), size = 0.25, color = "grey70") +
    ggplot2::geom_point(shape = 21, fill = "red", color = "black") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = bquote(~italic(.(x[1]))~"toxicity"~.(x[2])~"(μg/L)"),
                  y = bquote(~italic(.(y[1]))~"toxicity"~.(y[2])~"(μg/L)")) +
    ggplot2::scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::annotation_logticks()
  plot
  
  # setup output list
  out$all_focal_dat <- all_focal_dat
  out$plot_dat <- plot_dat
  out$orthog_reg_model_log10 <- orthog_reg_model_log10
  out$old_OR_dat <- old_OR_dat
  out$plot <- plot
  
  # return it
  return(out)
  
} 

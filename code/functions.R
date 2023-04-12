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
# Pull out mean squared error and r-squared from orthReg model
#======================================================================#
# r squared - https://stackoverflow.com/questions/75067630/orthogonal-linear-regression-total-least-squares-fit-get-rmse-and-r-squared-i
r_squared_odreg <- function(object, y) {
  denom <- sum((y - mean(y))^2)
  1 - object$ssq/denom
}
# meas squeared error - https://stackoverflow.com/questions/75067630/orthogonal-linear-regression-total-least-squares-fit-get-rmse-and-r-squared-i
mse_odreg <- function(object){
  mean(object$resid^2)
}

#======================================================================#
# function for data selection, shaping, orthogonal regression, plotting 
#======================================================================#
orthReg <- function(data, x, y, plot = F){
  # make a list to hold it all
  out <- NULL
  
  # filter data to focal taxa and and geometric mean for group
  all_focal_dat <- data %>%
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
  
  # reshape for plotting
  plot_dat <- all_focal_dat %>%
    dplyr::distinct(chem_name, pair_gen, gm_mean, min, max) %>% # just get geom mean and ranges
    tidyr::pivot_wider(names_from = pair_gen, values_from = c(gm_mean, min, max)) %>% # give us x and a y vars
    dplyr::filter(complete.cases(.)) # keep only complete cases
  
  # ORIGINAL REGRESSION BASED ON gm_mean
  orthog_reg_model_log10 <- pracma::odregress(x = log10(plot_dat$gm_mean_x), y = log10(plot_dat$gm_mean_y))
  
  # run the extraction functions
  orthog_reg_model_mse <- mse_odreg(orthog_reg_model_log10)
  orthog_reg_model_r_squared <- r_squared_odreg(orthog_reg_model_log10, log10(plot_dat$gm_mean_y))
  
  # build output
  orthreg_df <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                                   y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                                   orth.reg.n.observations = nrow(plot_dat),
                                   orth.reg.slope = orthog_reg_model_log10$coeff[1],
                                   orth.reg.intercept = orthog_reg_model_log10$coeff[2],
                                   orth.reg.ssq = orthog_reg_model_log10$ssq[1],
                                   orth.reg.mse = orthog_reg_model_mse,
                                   orth.reg.r.squared = orthog_reg_model_r_squared)
  
  if(plot == T){
  # plot it with log ticks
  plot <- ggplot2::ggplot(plot_dat) +
    ggplot2::aes(x = gm_mean_x, y = gm_mean_y) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
    ggplot2::geom_abline(slope = orthog_reg_model_log10$coeff[1], intercept = orthog_reg_model_log10$coeff[2], size = 0.5, color = "red") +
    ggplot2::geom_errorbar(aes(ymin = min_y, ymax = max_y), width = 0, size = 0.25, color = "grey70") +
    ggplot2::geom_errorbarh(aes(xmin = min_x, xmax = max_x), height = 0, size = 0.25, color = "grey70") +
    ggplot2::geom_point(shape = 21, fill = "red", color = "black") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = bquote(~italic(.(x[1]))~"toxicity"~.(x[2])~"(μg/L)"),
                  y = bquote(~italic(.(y[1]))~"toxicity"~.(y[2])~"(μg/L)"),
                  subtitle = glue::glue("{x[1]}_{x[2]}_{x[3]}_{x[4]} : {y[1]}_{y[2]}_{y[3]}_{y[4]}")) +
    ggplot2::scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::annotation_logticks()
  
  # setup output list
  out$orthreg_df <- orthreg_df
  out$plot <- plot
  
  # return list with data and plot
  return(out)
  }
  else{
    # return orthReg df
    return(orthreg_df)
  }
} 

#============================================================#
# Pairwise orthogonal Regression function
#============================================================#
pwOrthReg <- function(data, group, message = F){
  # lets get an id for the pair
  dat.id <- data %>%
    dplyr::rename_at(vars(matches(group)), ~ "group") %>%
    dplyr::mutate(id = paste(group, test_statistic, duration_d, endpoint, sep = ":"))
  
  # setup unique paris of group, test_stattistic, duration, endpoint
  unique.pairs <- dat.id %>%
    dplyr::distinct(id)
  
  pairs.list <- combn(x = unique.pairs$id, m = 2, simplify = F)  
  
  # Create the progress bar
  pb <- progress::progress_bar$new(total = length(pairs.list))
  
  # make a safe funciton for looping
  orthReg.safe <- purrr::safely(orthReg)
  
  # setup a list to hold output
  orthReg.list <- NULL
  
  # loop over pairs
  for(i in 1:length(pairs.list)){
    # get the pair data we care about
    pair.dat <- dat.id %>%
      dplyr::filter(id %in% pairs.list[[i]])
    
    # setup x and y, set NA's to 0 - NEED TO BE CERTAIN THIS IS KOSHER 
    x <- stringr::str_split_fixed(pairs.list[[i]][1], n = 4, pattern = ":")
    y <- stringr::str_split_fixed(pairs.list[[i]][2], n = 4, pattern = ":")
    
    # make a message for running
    if(message == T){
      message(glue::glue("running orthReg for {i} of {length(pairs.list)} pairs - {paste0(pairs.list[[i]], collapse =' ')}"))
    }
    if(message == F){
      # Update the progress bar
      pb$tick()
    }
    # perform orthogonal regression without plotting
    orthReg.out <- orthReg.safe(data = pair.dat, x = x, y = y, plot = F)
    
    # handle orthReg.safe output with errors
    if(is.null(orthReg.out$result)){
      orthReg.out$result <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                                           y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                                           orth.reg.n.observations = NA_real_,
                                           orth.reg.slope = NA_real_,
                                           orth.reg.intercept = NA_real_,
                                           orth.reg.ssq = NA_real_,
                                           orth.reg.mse = NA_real_,
                                           orth.reg.r.squared = NA_real_)
    }
    # add to the output list
    orthReg.list[[i]] <- orthReg.out$result
  }
  
  # return
  return(orthReg.list)
}

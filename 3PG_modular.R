# Modular 3PG TinyDancer Analysis
# Main function approach for multi-wave history matching of 3PG forest model

####################################
###Load Required Libraries###
####################################
source("R/InitialiseRunFunctions.R")
library(SoNWaL)
library(lubridate)
library(coda)
library(BayesianTools)
library(miscTools)
library(ggpubr)
library(matrixStats)
library(future)
library(furrr)
library(parallel)
library(lhs)
library(dgpsi)
library(tinydancer)
library(tictoc)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(tidyr)
library(GGally)

####################################
###Core 3PG Functions###
####################################

# Initialize 3PG parameters and climate data
initialize_3pg_data <- function(timeStep = "weekly") {
  
  # Get climate data
  clm_df_full <- read.csv(paste0("data/clm_df_", timeStep, "_Harwood_CHESS-MET_1973-2024.csv"))
  
  # Get calibration data from Harwood
  flxdata_daily <- read.csv("data/harwood_daily_2015_2020.csv")
  
  # Load default parameters for Sitka 
  sitka <- getParmsSitka(weather = clm_df_full,
                        waterBalanceSubMods = T,
                        timeStp = timeStep)
  
  # Parameter names to calibrate
  nm <- c("sigma_zR", "E_S1", "E_S2", "shared_area", "maxRootDepth", "K_drain", 
          "pFS2", "pFS20", "aS", "nS", "pRx", "pRn", "gammaFx", "gammaF0", "tgammaF", 
          "Rttover", "mF", "mR", "mS", "SLA0", "SLA1", "tSLA", "k", "alpha", "Y", "totNx",
          "kYmax", "kOmax", "hc", "qi", "qh", "qb", "eY", "fYC", "fYN")
  
  # Lower bounds
  f.decrease <- c(0.1, 0.001, 0.01, 1, 0.5, 0.01, 0.3, 0.1, 0.01, 2, 0.3, 0.2, 
                  0.002, 0.00005, 20, 0.003, 0.001, 0.1, 0.1, 4, 2, 2, 0.4, 0.03, 
                  0.43, 2, 0.001, 0.00001, 0.01, 90, 3, 0.6, 0.01, 0.1, 0.1)
  
  # Upper bounds
  f.increase <- c(3, 1, 1, 6, 2, 1, 1.6, 1.5, 0.3, 3, 1, 0.4, 0.1, 0.05, 110, 
                  0.16, 0.8, 0.8, 0.5, 9, 6, 10, 0.7, 0.06, 0.49, 20, 0.1, 0.01, 
                  0.8, 600, 50, 20, 0.9, 0.9, 0.9)
  
  # Create parameter bounds dataframe
  para <- data.frame(rbind(f.decrease, f.increase))
  colnames(para) <- nm
  
  # Subset of parameters for calibration
  para_sub <- para[, c("shared_area", "K_drain", "pFS2", "pFS20", "gammaFx", "SLA0", "SLA1", 
                      "tSLA", "k", "alpha", "Y", "fYC", "fYN", "totNx", "kYmax", "kOmax", 
                      "hc", "qi", "qh", "qb", "eY")]
  
  # NPP modifier
  modif <- 7.142857
  
  return(list(
    sitka = sitka,
    para_sub = para_sub,
    modif = modif,
    clm_df_full = clm_df_full,
    flxdata_daily = flxdata_daily
  ))
}

# Generate initial training data
#Modify 2000
generate_initial_data <- function(para_sub, sitka, modif, n_samples = 100) {
  
  cat("Generating initial training data...\n")
  
  # Generate LHS design
  lhd <- lhs::maximinLHS(n_samples, ncol(para_sub))
  for (i in 1:ncol(lhd)) {
    min_val <- para_sub[1, i]
    max_val <- para_sub[2, i]
    lhd[, i] <- lhd[, i] * (max_val - min_val) + min_val
  }
  
  lhd <- data.frame(lhd)
  names(lhd) <- names(para_sub)
  
  # Run 3PG for initial design
  all_output <- data.frame()
  n <- nrow(lhd)
  results <- vector("list", n)
  
  for (i in seq_len(n)) {
    sitka[colnames(para_sub)] <- lhd[i, ]
    output <- do.call(SoNWaL, sitka)
    output$id <- i
    results[[i]] <- output[, c(1, 2, 3, 46, 103)]
    if (i %% 100 == 0) cat("Iteration:", i, "\n")
  }
  
  all_output <- rbindlist(results)
  
  # Process output
  all_output <- all_output %>%
    mutate(NPP = NPP * modif) %>%
    group_by(Year, Month, id) %>%
    summarise(NPP = mean(NPP))
  
  # Filter for calibration years and months
  all_output_sum <- all_output[all_output$Year %in% c(2016, 2018), ]
  all_output_sum <- na.omit(all_output_sum)
  all_output_sum <- all_output_sum %>% mutate(YM = paste(Year, Month, sep = ""))
  
  # Normalize input data for emulators
  lhd_norm <- lhd
  for (i in 1:ncol(lhd_norm)) {
    min_val <- para_sub[1, i]
    max_val <- para_sub[2, i]
    lhd_norm[, i] <- (lhd_norm[, i] - min_val) / (max_val - min_val)
  }
  
  return(list(
    design_original = lhd,
    design_normalized = lhd_norm,
    output_data = all_output_sum,
    X = unname(as.matrix(lhd_norm))
  ))
}

# Build initial emulators for all months
build_initial_emulators <- function(initial_data) {
  
  cat("Building initial emulators...\n")
  
  # Month combinations for calibration
  months <- c("20181", "20183", "20185", "20187", "20189", 
              "20161", "20163", "20165", "20167", "20169")
  
  emulators <- list()
  Y_data <- list()
  
  for (month in months) {
    cat(sprintf("Building emulator for %s...\n", month))
    
    # Extract NPP data for this month
    month_data <- initial_data$output_data[initial_data$output_data$YM == month, ]$NPP
    Y_month <- unname(as.matrix(month_data))
    
    # Build emulator
    emulator <- dgp(initial_data$X, Y_month, name = c('matern2.5', 'sexp'), 
                   vecchia = T, verb = F, node = 1, share = F)
    
    emulators[[month]] <- emulator
    Y_data[[month]] <- Y_month
  }
  
  return(list(emulators = emulators, Y_data = Y_data))
}

# General implausibility function
general_imp <- function(NewData, Emulator, Obs, Disc, ObsErr) {
  tem <- predict(Emulator, NewData)$results
  abs(Obs - tem$mean) / sqrt(tem$var + Disc + ObsErr)
}

# Multi-month implausibility function
multimonth_imp <- function(x, targetLevel, levels, waves, 
                          EmulatorList = NULL,
                          ObsList = NULL, 
                          DiscList = NULL, 
                          ObsErrList = NULL) {
  
  # Set defaults if not provided
  if (is.null(ObsList)) ObsList <- rep(1.0, waves)
  if (is.null(DiscList)) DiscList <- rep(0.01, waves)
  if (is.null(ObsErrList)) ObsErrList <- rep(0, waves)
  
  ans <- rep(Inf, levels)
  waveFail <- FALSE
  this.level <- 1
  wave.num <- 1
  Timp <- NA
  
  while ((this.level <= levels) & !waveFail) {
    Timp <- general_imp(
      NewData = x, 
      Emulator = EmulatorList[[wave.num]], 
      Obs = ObsList[wave.num],
      Disc = DiscList[wave.num],
      ObsErr = ObsErrList[wave.num]
    )
    wave.num <- wave.num + 1
    if ((Timp > targetLevel[this.level]) & (wave.num <= waves)) {
      waveFail <- TRUE
    }
    if ((!waveFail) & (wave.num > waves)) {
      ans[this.level:levels] <- Timp
      this.level <- levels + 1
    } else {
      ans[this.level] <- Timp
      this.level <- this.level + 1
    }
  }
  return(ans)
}

# 3PG evaluation function for new design points
evaluate_3pg <- function(design_points, para_sub, sitka, modif) {
  
  cat(sprintf("Running 3PG for %d design points...\n", nrow(design_points)))
  
  results <- data.frame()
  
  for (j in 1:nrow(design_points)) {
    x <- as.data.frame(design_points[j, ])
    
    # Scale parameters back to original range
    for (i in 1:ncol(x)) {
      min_val <- para_sub[1, i]
      max_val <- para_sub[2, i]
      x[i] <- x[i] * (max_val - min_val) + min_val
    }
    
    # Update sitka parameters and run model
    sitka[c(colnames(para_sub))] <- x
    output <- do.call(SoNWaL, sitka)
    
    # Process output
    output <- output %>%
      mutate(NPP = NPP * modif) %>%
      group_by(Year, Month) %>%
      summarise(NPP = mean(NPP))
    
    output <- output %>% mutate(YM = paste(Year, Month, sep = ""))
    
    results <- rbind(results, as.data.frame(output[, c("NPP", "YM")]))
    
    if (j %% 100 == 0) cat("Completed:", j, "/", nrow(design_points), "\n")
  }
  
  return(results)
}

# Process a single wave
#same, modify 2000
process_wave <- function(wave_num, EmulatorList, target_level, control_list, 
                        my_obs_values, my_disc_values, my_obserr_values,
                        para_sub, sitka, modif, n_iter = 100) {
  
  cat(sprintf("\n=== Processing Wave %d ===\n", wave_num))
  
  # Create implausibility function for this wave
  implausibility_func <- function(x, target_level) {
    multimonth_imp(x, target_level, 1, 1, 
                  EmulatorList = EmulatorList,
                  ObsList = my_obs_values,
                  DiscList = my_disc_values,
                  ObsErrList = my_obserr_values)
  }
  
  # Construct temperature ladder
  cat("Constructing temperature ladder...\n")
  new_ladder <- construct_temperature_ladder(
    implausibility = implausibility_func,
    dims = 21,
    target_levels = target_level,
    control_list = control_list
  )
  
  # Sample using parallel tempering MCMC
  cat("Sampling using PTMCMC...\n")
  result <- sample_slice_ptmcmc(new_ladder, 
                               n_iter = n_iter,
                               implausibility = implausibility_func,
                               control_list = control_list)
  
  # Extract design points for next wave
  n_chains <- length(result$uniform_sample_list)
  wave_des <- result$uniform_sample_list[[n_chains]][, 1:21]
  wave_des <- as.data.frame(wave_des)
  colnames(wave_des) <- names(para_sub)
  
  # Evaluate 3PG function on new design points
  cat("Evaluating 3PG model...\n")
  wave_output <- evaluate_3pg(wave_des, para_sub, sitka, modif)
  
  # Build new emulators for each month
  cat("Training new emulators...\n")
  months <- c("20181", "20183", "20185", "20187", "20189", 
              "20161", "20163", "20165", "20167", "20169")
  
  new_emulators <- list()
  for (month in months) {
    month_data <- wave_output[wave_output$YM == month, ]$NPP
    if (length(month_data) > 0) {
      emulator <- dgp(as.matrix(wave_des), as.matrix(month_data), 
                     name = c('matern2.5', 'sexp'), vecchia = T, verb = F, node = 1, share = F)
      new_emulators[[month]] <- emulator
    }
  }
  
  # Add ID column for plotting
  wave_output$ID <- rep(1:nrow(wave_des), each = length(months))
  
  return(list(
    ladder = new_ladder,
    result = result,
    design_points = wave_des,
    output_data = wave_output,
    emulators = new_emulators
  ))
}

# Create plots for a wave
create_wave_plots <- function(wave_data, wave_num, observed_NPP = NULL) {
  
  # Prepare data for plotting
  filtered_output <- wave_data$output_data %>%
    filter(grepl("^2015|^2016|^2017|^2018", YM)) %>%
    mutate(
      YM = as.character(YM),
      year = substr(YM, 1, 4),
      month = substr(YM, 5, nchar(YM)),
      month = sprintf("%02d", as.numeric(month)),
      date = ymd(paste0(year, "-", month, "-01"))
    ) %>%
    select(-year, -month)
  
  # Time series plot
  ts_plot <- ggplot(filtered_output, aes(x = date, y = NPP, group = ID)) + 
    geom_line(alpha = 0.3, color = "black") +
    theme_minimal() +
    labs(title = sprintf("Wave %d: 3PG Output Time Series", wave_num),
         x = "Date", y = "NPP") +
    theme(legend.position = "none")
  
  if (!is.null(observed_NPP)) {
    ts_plot <- ts_plot + 
      geom_line(data = observed_NPP, aes(x = date, y = NPP), 
                color = "red", linewidth = 2, inherit.aes = FALSE)
  }
  
  # Parameter space pairplot (first 10 parameters for visibility)
  param_plot <- ggpairs(wave_data$design_points[, 1:10],
                       lower = list(continuous = wrap("points", alpha = 0.3, size = 0.1)),
                       title = sprintf("Wave %d: Parameter Space", wave_num))
  
  # Sample one emulator plot
  first_month <- names(wave_data$emulators)[1]
  emulator_plot <- plot(wave_data$emulators[[first_month]], style = 1) + 
    ggtitle(sprintf("Wave %d Emulator (%s)", wave_num, first_month))
  
  return(list(
    timeseries_plot = ts_plot,
    parameter_plot = param_plot,
    emulator_plot = emulator_plot
  ))
}

# Get observed data for comparison
get_observed_data <- function(flxdata_daily, timeStep = "monthly", startYear = 2015, endYear = 2018) {
  
  observedVals_npp <- function(timeStep, data, sY = 2015, eY = 2018) {
    data <- filter(data, year >= sY & year <= eY)
    data <- data[-1, ]
    
    if (timeStep == "weekly") {
      data <- data %>%
        mutate(grp = week(as.Date(data$yday, origin = paste0(data$year, "-01-01"))))
    }
    
    if (timeStep == "monthly") {
      data <- data %>%
        mutate(grp = month(as.Date(data$yday, origin = paste0(data$year, "-01-01"))))
    }
    
    observed <- c(pull(data %>% 
                      dplyr::group_by(year, grp) %>%
                      dplyr::summarise(npp = mean(npp)) %>%
                      dplyr::select(npp)))
    
    return(observed)
  }
  
  observed <- observedVals_npp(timeStep = timeStep, data = flxdata_daily, 
                              sY = startYear, eY = endYear)
  
  # Create date sequence
  dates <- seq(ymd(paste0(startYear, "-01-01")), 
               ymd(paste0(endYear, "-12-01")), by = "month")
  
  observed_NPP <- data.frame(
    NPP = observed[1:length(dates)],
    date = dates,
    ID = "observed",
    YM = 1
  )
  
  return(observed_NPP)
}

# Main function to run multi-wave 3PG analysis
run_3pg_waves <- function(n_waves, 
                         target_levels = NULL,
                         n_iter = 2000,
                         control_params = NULL,
                         timeStep = "weekly") {
  
  cat("=== Starting 3PG Multi-Wave Analysis ===\n")
  cat(sprintf("Number of waves: %d\n", n_waves))
  
  # Default target levels (progressively tighter)
  if (is.null(target_levels)) {
    if (n_waves == 1) {
      target_levels <- 0.25
    } else if (n_waves == 2) {
      target_levels <- c(0.25, 0.1)
    } else if (n_waves == 3) {
      target_levels <- c(0.25, 0.1, 0.05)
    } else {
      # For more waves, create a sequence
      target_levels <- c(0.25, seq(0.1, 0.05, length.out = n_waves - 1))
    }
  }
  
  # Default control parameters
  if (is.null(control_params)) {
    control_params <- list(
      num_mutations = 10,
      num_iterations = 100,
      box_limits = matrix(rep(c(0, 1), 21), nrow = 21, byrow = TRUE),
      fork_chains = FALSE,
      estimate_volume = TRUE
    )
  }
  
  # Initialize 3PG data and Python
  init_py()
  data_setup <- initialize_3pg_data(timeStep)
  
  # Generate initial training data
  initial_data <- generate_initial_data(data_setup$para_sub, data_setup$sitka, 
                                       data_setup$modif)
  
  # Build initial emulators
  initial_emulators <- build_initial_emulators(initial_data)
  
  # Get observed data for plotting
  observed_NPP <- get_observed_data(data_setup$flxdata_daily)
  
  # Observed values for implausibility (you may want to adjust these)
  my_obs_values <- c(0.7394801260826509, 2.6577939192156195, 4.222251009108864, 
                     6.578773935533921, 4.185232760947738, 1.0677007185458087,
                     2.6150432429649597, 4.998384222085997, 5.442871183494439, 
                     3.6012050730222933)
  my_disc_values <- rep(0.01, 10)
  my_obserr_values <- rep(0, 10)
  
  # Initialize storage for results
  waves_data <- list()
  current_emulators <- initial_emulators$emulators
  all_plots <- list()
  
  cat(sprintf("Target levels: %s\n", paste(target_levels, collapse = ", ")))
  
  # Process each wave
  for (wave in 1:n_waves) {
    
    # For each subsequent wave, update emulator list following TinyDancer pattern
    if (wave > 1) {
      # Add previous wave emulators and repeat the last set
      months <- names(waves_data[[wave-1]]$emulators)
      for (month in months) {
        current_emulators[[paste0(month, "_wave", wave-1)]] <- waves_data[[wave-1]]$emulators[[month]]
      }
    }
    
    # Convert to list format expected by implausibility function
    emulator_list <- as.list(current_emulators)
    
    # Adjust observation values to match number of emulators
    n_emulators <- length(emulator_list)
    obs_values <- rep(my_obs_values, ceiling(n_emulators / length(my_obs_values)))[1:n_emulators]
    disc_values <- rep(my_disc_values, ceiling(n_emulators / length(my_disc_values)))[1:n_emulators]
    obserr_values <- rep(my_obserr_values, ceiling(n_emulators / length(my_obserr_values)))[1:n_emulators]
    
    # Process the wave
    wave_data <- process_wave(
      wave_num = wave,
      EmulatorList = emulator_list,
      target_level = target_levels[wave],
      control_list = control_params,
      my_obs_values = obs_values,
      my_disc_values = disc_values,
      my_obserr_values = obserr_values,
      para_sub = data_setup$para_sub,
      sitka = data_setup$sitka,
      modif = data_setup$modif,
      n_iter = n_iter
    )
    
    # Store results
    waves_data[[wave]] <- wave_data
    
    # Create plots
    wave_plots <- create_wave_plots(wave_data, wave, observed_NPP)
    all_plots[[wave]] <- wave_plots
    
    # Print summary
    cat(sprintf("Wave %d completed. Generated %d design points.\n", 
                wave, nrow(wave_data$design_points)))
  }
  
  # Create summary plots
  summary_plots <- create_summary_plots(waves_data, all_plots, n_waves, observed_NPP)
  
  # Return comprehensive results
  results <- list(
    waves_data = waves_data,
    initial_emulators = initial_emulators,
    plots = all_plots,
    summary_plots = summary_plots,
    target_levels = target_levels,
    initial_data = initial_data,
    data_setup = data_setup,
    observed_NPP = observed_NPP
  )
  
  cat("\n=== Analysis Complete ===\n")
  cat(sprintf("Completed %d waves of 3PG history matching.\n", n_waves))
  
  return(results)
}

# Create summary plots across waves
create_summary_plots <- function(waves_data, all_plots, n_waves, observed_NPP) {
  
  # Convergence plot - number of design points per wave
  n_points <- sapply(waves_data, function(w) nrow(w$design_points))
  convergence_data <- data.frame(
    Wave = 1:n_waves,
    DesignPoints = n_points
  )
  
  convergence_plot <- ggplot(convergence_data, aes(x = Wave, y = DesignPoints)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "red", size = 2) +
    theme_minimal() +
    labs(title = "Design Points Generated Across Waves",
         x = "Wave Number",
         y = "Number of Design Points")
  
  # Combined time series plots
  combined_ts_plots <- wrap_plots(lapply(all_plots, function(p) p$timeseries_plot), 
                                 ncol = min(n_waves, 2))
  
  # Final wave comparison with observed
  if (n_waves > 0) {
    final_wave_data <- waves_data[[n_waves]]$output_data
    final_filtered <- final_wave_data %>%
      filter(grepl("^2015|^2016|^2017|^2018", YM)) %>%
      mutate(
        YM = as.character(YM),
        year = substr(YM, 1, 4),
        month = substr(YM, 5, nchar(YM)),
        month = sprintf("%02d", as.numeric(month)),
        date = ymd(paste0(year, "-", month, "-01"))
      ) %>%
      select(-year, -month)
    
    final_comparison <- ggplot(final_filtered, aes(x = date, y = NPP, group = ID)) + 
      geom_line(alpha = 0.1, color = "black") +
      geom_line(data = observed_NPP, aes(x = date, y = NPP), 
                color = "red", linewidth = 2, inherit.aes = FALSE) +
      theme_minimal() +
      labs(title = "Final Wave: Model vs Observed",
           x = "Date", y = "NPP") +
      theme(legend.position = "none")
  } else {
    final_comparison <- ggplot() + ggtitle("No waves completed")
  }
  
  return(list(
    convergence = convergence_plot,
    timeseries_evolution = combined_ts_plots,
    final_comparison = final_comparison
  ))
}

# Convenience function to display results
display_results <- function(results) {
  cat("\n=== 3PG RESULTS SUMMARY ===\n")
  cat(sprintf("Number of waves: %d\n", length(results$waves_data)))
  cat(sprintf("Total emulators built: %d\n", 
              length(results$initial_emulators$emulators) + 
              sum(sapply(results$waves_data, function(w) length(w$emulators)))))
  
  for (i in 1:length(results$waves_data)) {
    cat(sprintf("Wave %d: %d design points, target level %.3f\n", 
                i, nrow(results$waves_data[[i]]$design_points), 
                results$target_levels[i]))
  }
  
  # Display final comparison plot
  print(results$summary_plots$final_comparison)
  
  return(invisible(results))
}

# Example usage:
results <- run_3pg_waves(n_waves = 3, target_levels = c(0.25, 0.1, 0.05))
display_results(results) 

results <- run_3pg_waves(n_waves = 1, target_levels = 0.25)
display_results(results) 


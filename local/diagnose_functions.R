## ============================================================
## SDA / EnKF Diagnostic Pipeline
## Author: Yang Gu workflow helper
## Purpose:
##   1. Read SDA output files from a directory.
##   2. Build final_summary_all, diag_summary_all, and obs_long_all.
##   3. Compute R2/RMSE/bias diagnostics.
##   4. Generate scatter plots, NIS figures, R2 visualizations,
##      and optional filter diagnostics.
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

## ============================================================
## 0. Utility functions
## ============================================================

calc_R2_sse <- function(obs, pred) {
  ## SSE-based R2:
  ## R2 = 1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2)
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  
  if (length(obs) < 2) return(NA_real_)
  
  sst <- sum((obs - mean(obs))^2)
  if (!is.finite(sst) || sst == 0) return(NA_real_)
  
  1 - sum((pred - obs)^2) / sst
}

safe_log10 <- function(x) {
  dplyr::if_else(is.finite(x) & x > 0, log10(x), NA_real_)
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", as.character(x))
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_plot_if_needed <- function(plot, filename, save_plot = TRUE, width = 12, height = 8, dpi = 300) {
  if (isTRUE(save_plot)) {
    ensure_dir(dirname(filename))
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  invisible(plot)
}

## ============================================================
## 1. Summarize SDA outputs
## ============================================================

pick_block_date_by_year <- function(block_all, year, date_key = NULL) {
  ## Select the block.list.all date matching the target year.
  non_null_dates <- names(block_all)[!vapply(block_all, is.null, logical(1))]
  
  if (length(non_null_dates) == 0) {
    stop("No non-null dates found in block.list.all.")
  }
  
  if (!is.null(date_key)) {
    if (!date_key %in% non_null_dates) {
      stop(
        "date_key not found in block.list.all: ", date_key,
        "\nAvailable dates: ", paste(non_null_dates, collapse = ", ")
      )
    }
    return(date_key)
  }
  
  block_dates <- suppressWarnings(as.Date(non_null_dates))
  
  ## Fallback when names include time stamps or extra characters.
  if (all(is.na(block_dates))) {
    block_dates <- suppressWarnings(as.Date(substr(non_null_dates, 1, 10)))
  }
  
  matched_dates <- non_null_dates[
    !is.na(block_dates) &
      as.integer(format(block_dates, "%Y")) == as.integer(year)
  ]
  
  if (length(matched_dates) == 0) {
    stop(
      "No block.list.all date found for year = ", year,
      "\nAvailable dates: ", paste(non_null_dates, collapse = ", ")
    )
  }
  
  ## If multiple dates exist in one year, use the last one.
  matched_dates[length(matched_dates)]
}

extract_site_diag_from_block <- function(blk, site_index, site_id, state_vars, year, date_key) {
  ## Extract diagonal Pf, Q, R for one local block/site.
  ## Notes:
  ## - blk$data$r is assumed to be observation precision, so R = solve(r).
  ## - q.type == 4: Wishart precision q; approximate Q = aq / bq following the existing workflow.
  ## - q.type == 3: vector precision q; E[q] = aq / bq, so Q = diag(1 / E[q]).
  
  H_mat <- blk$H
  
  obs_idx <- blk$constant$H
  if (is.null(obs_idx)) {
    obs_idx <- which(colSums(H_mat) != 0)
  }
  
  obs_vars <- state_vars[obs_idx]
  
  ## Observation covariance R.
  R <- tryCatch(
    solve(as.matrix(blk$data$r)),
    error = function(e) {
      warning("Failed to invert R precision for site ", site_id, ": ", conditionMessage(e))
      matrix(NA_real_, nrow = length(obs_vars), ncol = length(obs_vars))
    }
  )
  rownames(R) <- colnames(R) <- obs_vars
  
  ## Forecast covariance Pf.
  Pf <- as.matrix(blk$data$pf)
  rownames(Pf) <- colnames(Pf) <- state_vars
  
  ## Use posterior-updated Q hyperparameters when available.
  aq <- blk$update$aq
  bq <- blk$update$bq
  if (is.null(aq) || is.null(bq)) {
    aq <- blk$data$aq
    bq <- blk$data$bq
  }
  
  q_type <- blk$constant$q.type
  
  if (q_type == 4) {
    ## Wishart q case.
    Q <- as.matrix(aq) / as.numeric(bq)
    rownames(Q) <- colnames(Q) <- obs_vars
    
  } else if (q_type == 3) {
    ## Vector q case.
    q_mean <- as.numeric(aq) / as.numeric(bq)
    Q <- diag(1 / q_mean)
    rownames(Q) <- colnames(Q) <- obs_vars
    
  } else {
    stop("Unsupported q.type: ", q_type)
  }
  
  Pf_diag_all <- data.frame(
    year = year,
    date_key = date_key,
    site_index = site_index,
    site_id = as.character(site_id),
    variable = state_vars,
    Pf = as.numeric(diag(Pf)),
    row.names = NULL
  )
  
  QR_diag <- data.frame(
    year = year,
    date_key = date_key,
    site_index = site_index,
    site_id = as.character(site_id),
    variable = obs_vars,
    Q = as.numeric(diag(Q)),
    R = as.numeric(diag(R)),
    row.names = NULL
  ) %>%
    mutate(
      Q_over_R = Q / R
    )
  
  diag_summary <- Pf_diag_all %>%
    left_join(
      QR_diag,
      by = c("year", "date_key", "site_index", "site_id", "variable")
    ) %>%
    mutate(
      Pf_over_R = Pf / R
    )
  
  list(
    site_id = as.character(site_id),
    obs_vars = obs_vars,
    H = H_mat,
    Q = Q,
    R = R,
    Pf = Pf,
    diag_summary = diag_summary
  )
}

summarize_sda_output <- function(
    sda.outputs,
    year,
    date_key = NULL,
    state_vars = c("AbvGrndWood", "NEE", "Qle", "LAI", "SoilMoistFrac", "TotSoilCarb")
) {
  ## Summarize one SDA output object into long-format diagnostics.
  
  forecast <- as.data.frame(sda.outputs[["forecast"]])
  analysis <- as.data.frame(sda.outputs[["analysis"]])
  obs_mean <- sda.outputs[["obs.mean"]]
  block_all <- sda.outputs[["enkf.params"]][["block.list.all"]]
  
  date_key <- pick_block_date_by_year(
    block_all = block_all,
    year = year,
    date_key = date_key
  )
  
  blk_list <- block_all[[date_key]]
  
  nvar <- length(state_vars)
  nens <- nrow(forecast)
  nsite <- ncol(forecast) / nvar
  
  stopifnot(ncol(forecast) %% nvar == 0)
  stopifnot(all(dim(forecast) == dim(analysis)))
  stopifnot(length(blk_list) == nsite)
  
  site_ids <- vapply(
    blk_list,
    function(x) paste(as.character(x$site.ids), collapse = "_"),
    character(1)
  )
  
  ## Map each forecast/analysis column to year-site-variable.
  col_map <- data.frame(
    year = year,
    date_key = date_key,
    site_index = rep(seq_len(nsite), each = nvar),
    site_id = rep(site_ids, each = nvar),
    variable = rep(state_vars, times = nsite),
    col_index = seq_len(ncol(forecast)),
    row.names = NULL
  )
  
  state_mean <- col_map %>%
    mutate(
      forecast_mean = as.numeric(colMeans(forecast, na.rm = TRUE)),
      analysis_mean = as.numeric(colMeans(analysis, na.rm = TRUE)),
      analysis_minus_forecast = analysis_mean - forecast_mean
    )
  
  ## Convert obs.mean list to long format.
  obs_long <- bind_rows(lapply(names(obs_mean), function(sid) {
    x <- obs_mean[[sid]]
    
    if (is.null(x) || ncol(x) == 0) {
      return(data.frame())
    }
    
    data.frame(
      year = year,
      date_key = date_key,
      site_id = as.character(sid),
      variable = colnames(x),
      obs_value = as.numeric(x[1, ]),
      row.names = NULL
    )
  }))
  
  matrices_by_site <- lapply(seq_along(blk_list), function(i) {
    extract_site_diag_from_block(
      blk = blk_list[[i]],
      site_index = i,
      site_id = site_ids[i],
      state_vars = state_vars,
      year = year,
      date_key = date_key
    )
  })
  names(matrices_by_site) <- site_ids
  
  diag_summary <- bind_rows(lapply(matrices_by_site, `[[`, "diag_summary"))
  
  final_summary <- state_mean %>%
    left_join(
      obs_long,
      by = c("year", "date_key", "site_id", "variable")
    ) %>%
    left_join(
      diag_summary,
      by = c("year", "date_key", "site_index", "site_id", "variable")
    ) %>%
    mutate(
      forecast_minus_obs = forecast_mean - obs_value,
      analysis_minus_obs = analysis_mean - obs_value,
      abs_forecast_error = abs(forecast_minus_obs),
      abs_analysis_error = abs(analysis_minus_obs),
      increment = analysis_mean - forecast_mean,
      error_improvement = abs_forecast_error - abs_analysis_error,
      worsened = abs_analysis_error > abs_forecast_error
    )
  
  list(
    year = year,
    date_key = date_key,
    nens = nens,
    nsite = nsite,
    final_summary = final_summary,
    diag_summary = diag_summary,
    obs_long = obs_long,
    matrices_by_site = matrices_by_site
  )
}

summarize_sda_directory <- function(
    sda_dir,
    file_index = 1:13,
    start_year = 2012,
    state_vars = c("AbvGrndWood", "NEE", "Qle", "LAI", "SoilMoistFrac", "TotSoilCarb"),
    save_csv = TRUE,
    out_dir = sda_dir,
    prefix = "SDA_all_years"
) {
  ## Read multiple sda.output*.Rdata files and combine diagnostics.
  
  all_results <- list()
  
  for (i in file_index) {
    year_i <- start_year + i - 1
    file_i <- file.path(sda_dir, paste0("sda.output", i, ".Rdata"))
    
    if (!file.exists(file_i)) {
      warning("File not found: ", file_i)
      next
    }
    
    message("Reading: ", file_i, " | year = ", year_i)
    
    env <- new.env(parent = emptyenv())
    load(file_i, envir = env)
    
    if (!exists("sda.outputs", envir = env)) {
      warning("No object named sda.outputs in: ", file_i)
      next
    }
    
    res_i <- summarize_sda_output(
      sda.outputs = env$sda.outputs,
      year = year_i,
      date_key = NULL,
      state_vars = state_vars
    )
    
    all_results[[as.character(year_i)]] <- res_i
  }
  
  if (length(all_results) == 0) {
    stop("No valid SDA output files were summarized.")
  }
  
  final_summary_all <- bind_rows(lapply(all_results, `[[`, "final_summary"))
  diag_summary_all <- bind_rows(lapply(all_results, `[[`, "diag_summary"))
  obs_long_all <- bind_rows(lapply(all_results, `[[`, "obs_long"))
  
  if (isTRUE(save_csv)) {
    ensure_dir(out_dir)
    
    readr::write_csv(
      final_summary_all,
      file.path(out_dir, paste0(prefix, "_final_summary.csv"))
    )
    
    readr::write_csv(
      diag_summary_all,
      file.path(out_dir, paste0(prefix, "_diag_summary.csv"))
    )
    
    readr::write_csv(
      obs_long_all,
      file.path(out_dir, paste0(prefix, "_obs_long.csv"))
    )
  }
  
  list(
    final_summary_all = final_summary_all,
    diag_summary_all = diag_summary_all,
    obs_long_all = obs_long_all,
    yearly_results = all_results
  )
}

## ============================================================
## 2. Prepare plot-scale data
## ============================================================

prepare_diagnostic_data <- function(final_summary_all, soilmoist_div10 = FALSE) {
  ## Create plotting-scale forecast/analysis/Q/R/Pf columns.
  ## If soilmoist_div10 = TRUE, SoilMoistFrac means are divided by 10,
  ## and variance terms are divided by 100.
  
  final_summary_all %>%
    mutate(
      forecast_plot = if_else(
        soilmoist_div10 & variable == "SoilMoistFrac",
        forecast_mean / 10,
        forecast_mean
      ),
      analysis_plot = if_else(
        soilmoist_div10 & variable == "SoilMoistFrac",
        analysis_mean / 10,
        analysis_mean
      ),
      obs_plot = obs_value,
      Pf_plot = if_else(
        soilmoist_div10 & variable == "SoilMoistFrac",
        Pf / 100,
        Pf
      ),
      Q_plot = if_else(
        soilmoist_div10 & variable == "SoilMoistFrac",
        Q / 100,
        Q
      ),
      R_plot = if_else(
        soilmoist_div10 & variable == "SoilMoistFrac",
        R / 100,
        R
      ),
      Q_over_R_plot = Q_plot / R_plot,
      Pf_over_R_plot = Pf_plot / R_plot,
      forecast_resid = forecast_plot - obs_plot,
      analysis_resid = analysis_plot - obs_plot,
      increment_plot = analysis_plot - forecast_plot,
      abs_forecast_error = abs(forecast_resid),
      abs_analysis_error = abs(analysis_resid),
      error_improvement = abs_forecast_error - abs_analysis_error,
      worsened = abs_analysis_error > abs_forecast_error
    )
}

## ============================================================
## 3. Performance tables
## ============================================================

compute_performance_tables <- function(final_summary_all, soilmoist_div10 = FALSE) {
  ## Compute R2, RMSE, and bias for forecast and analysis.
  
  dat <- prepare_diagnostic_data(final_summary_all, soilmoist_div10 = soilmoist_div10)
  
  pred_long <- dat %>%
    filter(is.finite(obs_plot)) %>%
    select(year, site_id, variable, obs_plot, forecast_plot, analysis_plot) %>%
    pivot_longer(
      cols = c(forecast_plot, analysis_plot),
      names_to = "type",
      values_to = "pred"
    ) %>%
    mutate(
      type = recode(
        type,
        forecast_plot = "forecast",
        analysis_plot = "analysis"
      ),
      residual = pred - obs_plot,
      abs_residual = abs(residual)
    )
  
  performance_by_variable <- pred_long %>%
    group_by(variable, type) %>%
    summarise(
      n = n(),
      R2_sse = calc_R2_sse(obs_plot, pred),
      RMSE = sqrt(mean((pred - obs_plot)^2, na.rm = TRUE)),
      bias = mean(pred - obs_plot, na.rm = TRUE),
      MAE = mean(abs(pred - obs_plot), na.rm = TRUE),
      .groups = "drop"
    )
  
  performance_by_year_variable <- pred_long %>%
    group_by(year, variable, type) %>%
    summarise(
      n = n(),
      R2_sse = calc_R2_sse(obs_plot, pred),
      RMSE = sqrt(mean((pred - obs_plot)^2, na.rm = TRUE)),
      bias = mean(pred - obs_plot, na.rm = TRUE),
      MAE = mean(abs(pred - obs_plot), na.rm = TRUE),
      .groups = "drop"
    )
  
  performance_overall <- pred_long %>%
    group_by(type) %>%
    summarise(
      n = n(),
      R2_sse = calc_R2_sse(obs_plot, pred),
      RMSE = sqrt(mean((pred - obs_plot)^2, na.rm = TRUE)),
      bias = mean(pred - obs_plot, na.rm = TRUE),
      MAE = mean(abs(pred - obs_plot), na.rm = TRUE),
      .groups = "drop"
    )
  
  r2_wide_by_variable <- performance_by_variable %>%
    select(variable, type, n, R2_sse, RMSE, bias, MAE) %>%
    pivot_wider(
      names_from = type,
      values_from = c(n, R2_sse, RMSE, bias, MAE)
    ) %>%
    mutate(
      R2_improvement = R2_sse_analysis - R2_sse_forecast,
      RMSE_reduction = RMSE_forecast - RMSE_analysis,
      MAE_reduction = MAE_forecast - MAE_analysis
    )
  
  r2_wide_by_year_variable <- performance_by_year_variable %>%
    select(year, variable, type, n, R2_sse, RMSE, bias, MAE) %>%
    pivot_wider(
      names_from = type,
      values_from = c(n, R2_sse, RMSE, bias, MAE)
    ) %>%
    mutate(
      R2_improvement = R2_sse_analysis - R2_sse_forecast,
      RMSE_reduction = RMSE_forecast - RMSE_analysis,
      MAE_reduction = MAE_forecast - MAE_analysis
    )
  
  list(
    pred_long = pred_long,
    performance_by_variable = performance_by_variable,
    performance_by_year_variable = performance_by_year_variable,
    performance_overall = performance_overall,
    r2_wide_by_variable = r2_wide_by_variable,
    r2_wide_by_year_variable = r2_wide_by_year_variable
  )
}

write_performance_tables <- function(performance, out_dir) {
  ensure_dir(out_dir)
  
  readr::write_csv(performance$performance_by_variable, file.path(out_dir, "performance_by_variable.csv"))
  readr::write_csv(performance$performance_by_year_variable, file.path(out_dir, "performance_by_year_variable.csv"))
  readr::write_csv(performance$performance_overall, file.path(out_dir, "performance_overall.csv"))
  readr::write_csv(performance$r2_wide_by_variable, file.path(out_dir, "R2_wide_by_variable.csv"))
  readr::write_csv(performance$r2_wide_by_year_variable, file.path(out_dir, "R2_wide_by_year_variable.csv"))
  
  invisible(performance)
}

## ============================================================
## 4. Plot functions
## ============================================================

plot_obs_forecast_analysis <- function(
    final_summary_all,
    soilmoist_div10 = FALSE,
    save_plot = FALSE,
    out_file = "obs_forecast_analysis_scatter.png",
    point_size = 2.2,
    point_alpha = 0.70
) {
  ## Scatter plot: Observation vs forecast/analysis ensemble mean.
  
  performance <- compute_performance_tables(final_summary_all, soilmoist_div10 = soilmoist_div10)
  plot_dat <- performance$pred_long
  
  anno_pos <- plot_dat %>%
    group_by(variable) %>%
    summarise(
      x = min(obs_plot, na.rm = TRUE),
      y = max(pred, na.rm = TRUE),
      .groups = "drop"
    )
  
  anno_lab <- performance$performance_by_variable %>%
    select(variable, type, n, R2_sse) %>%
    pivot_wider(names_from = type, values_from = c(n, R2_sse)) %>%
    left_join(anno_pos, by = "variable") %>%
    mutate(
      label = paste0(
        "Forecast R² = ", sprintf("%.2f", R2_sse_forecast),
        "\nAnalysis R² = ", sprintf("%.2f", R2_sse_analysis),
        "\nN = ", n_forecast
      )
    )
  
  p <- ggplot(plot_dat, aes(x = obs_plot, y = pred, color = type)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.7,
      color = "gray40"
    ) +
    geom_point(size = point_size, alpha = point_alpha) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
    geom_text(
      data = anno_lab,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.5,
      color = "black"
    ) +
    facet_wrap(~ variable, scales = "free") +
    scale_color_manual(
      values = c("forecast" = "#2C7FB8", "analysis" = "#D7301F")
    ) +
    labs(
      x = "Observation",
      y = "Prediction",
      color = NULL,
      title = "Observed vs Forecast / Analysis",
      subtitle = ifelse(
        soilmoist_div10,
        "SoilMoistFrac was divided by 10 for plotting",
        "Raw scale"
      )
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 12, height = 8)
  
  list(
    plot = p,
    plot_data = plot_dat,
    r2_table = performance$performance_by_variable,
    annotation_table = anno_lab
  )
}

plot_R2_by_variable <- function(
    performance,
    save_plot = FALSE,
    out_file = "R2_by_variable.png"
) {
  ## Visualize forecast and analysis R2 by variable.
  
  plot_dat <- performance$performance_by_variable %>%
    mutate(
      type = factor(type, levels = c("forecast", "analysis"))
    )
  
  p <- ggplot(plot_dat, aes(x = variable, y = R2_sse, fill = type)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    coord_flip() +
    scale_fill_manual(values = c("forecast" = "#2C7FB8", "analysis" = "#D7301F")) +
    labs(
      x = NULL,
      y = "SSE-based R²",
      fill = NULL,
      title = "Forecast vs Analysis R² by Variable",
      subtitle = "Higher R² means better agreement with observations"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 10, height = 6)
  
  list(plot = p, data = plot_dat)
}

plot_R2_improvement_by_variable <- function(
    performance,
    save_plot = FALSE,
    out_file = "R2_improvement_by_variable.png"
) {
  ## Visualize R2_analysis - R2_forecast by variable.
  
  plot_dat <- performance$r2_wide_by_variable %>%
    mutate(
      improved = R2_improvement > 0
    )
  
  p <- ggplot(plot_dat, aes(x = variable, y = R2_improvement, fill = improved)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_col(width = 0.7, alpha = 0.9) +
    coord_flip() +
    scale_fill_manual(
      values = c("FALSE" = "#D7301F", "TRUE" = "#2C7FB8"),
      labels = c("FALSE" = "analysis worse", "TRUE" = "analysis better")
    ) +
    labs(
      x = NULL,
      y = "R² improvement = R²_analysis - R²_forecast",
      fill = NULL,
      title = "R² Improvement After Assimilation by Variable"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 10, height = 6)
  
  list(plot = p, data = plot_dat)
}

plot_R2_by_year_variable <- function(
    performance,
    save_plot = FALSE,
    out_file = "R2_by_year_variable.png"
) {
  ## Visualize R2 through time for each variable.
  
  plot_dat <- performance$performance_by_year_variable %>%
    mutate(
      type = factor(type, levels = c("forecast", "analysis"))
    )
  
  p <- ggplot(plot_dat, aes(x = year, y = R2_sse, color = type, linetype = type)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 2.2, na.rm = TRUE) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_color_manual(values = c("forecast" = "#2C7FB8", "analysis" = "#D7301F")) +
    labs(
      x = "Year",
      y = "SSE-based R²",
      color = NULL,
      linetype = NULL,
      title = "R² by Year and Variable"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 12, height = 8)
  
  list(plot = p, data = plot_dat)
}

plot_R2_improvement_by_year_variable <- function(
    performance,
    save_plot = FALSE,
    out_file = "R2_improvement_by_year_variable.png"
) {
  ## Visualize R2 improvement through time for each variable.
  
  plot_dat <- performance$r2_wide_by_year_variable %>%
    mutate(
      improved = R2_improvement > 0
    )
  
  p <- ggplot(plot_dat, aes(x = year, y = R2_improvement, fill = improved)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_col(alpha = 0.9) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_fill_manual(
      values = c("FALSE" = "#D7301F", "TRUE" = "#2C7FB8"),
      labels = c("FALSE" = "analysis worse", "TRUE" = "analysis better")
    ) +
    labs(
      x = "Year",
      y = "R² improvement",
      fill = NULL,
      title = "R² Improvement by Year and Variable",
      subtitle = "Positive values mean analysis improved over forecast"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 12, height = 8)
  
  list(plot = p, data = plot_dat)
}

## ============================================================
## 5. NIS diagnostics
## ============================================================

calc_scalar_nis <- function(final_summary_all, soilmoist_div10 = FALSE) {
  ## Compute scalar NIS = innovation^2 / (Pf + Q + R).
  ## This is valid as a scalar per year-site-variable diagnostic.
  
  prepare_diagnostic_data(final_summary_all, soilmoist_div10 = soilmoist_div10) %>%
    filter(
      is.finite(obs_plot),
      is.finite(forecast_plot),
      is.finite(Pf_plot),
      is.finite(Q_plot),
      is.finite(R_plot),
      Pf_plot >= 0,
      Q_plot >= 0,
      R_plot > 0,
      (Pf_plot + Q_plot + R_plot) > 0
    ) %>%
    mutate(
      innovation = forecast_plot - obs_plot,
      S = Pf_plot + Q_plot + R_plot,
      expected_sd = sqrt(S),
      NIS = innovation^2 / S,
      abs_normalized_innovation = sqrt(NIS),
      exceed_95 = NIS > stats::qchisq(0.95, df = 1),
      exceed_99 = NIS > stats::qchisq(0.99, df = 1)
    )
}

plot_NIS_by_year <- function(
    nis_dat,
    save_plot = FALSE,
    out_file = "NIS_by_year.png",
    y_log10 = FALSE
) {
  ## Boxplot of NIS by year and variable.
  
  p <- ggplot(nis_dat, aes(x = factor(year), y = NIS)) +
    geom_hline(
      yintercept = stats::qchisq(0.95, df = 1),
      linetype = "dashed",
      color = "orange",
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = stats::qchisq(0.99, df = 1),
      linetype = "dashed",
      color = "red",
      linewidth = 0.8
    ) +
    geom_boxplot(outlier.alpha = 0.5) +
    facet_wrap(~ variable, scales = "free_y") +
    labs(
      x = "Year",
      y = "NIS = innovation² / (Pf + Q + R)",
      title = "Normalized Innovation Squared (NIS) by Year",
      subtitle = "Orange = 95% chi-square threshold; Red = 99% threshold"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
  
  if (isTRUE(y_log10)) {
    p <- p +
      scale_y_continuous(trans = "log10") +
      labs(y = "NIS, log10 scale")
  }
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 12, height = 8)
  
  list(plot = p, data = nis_dat)
}

plot_NIS_exceedance_rate <- function(
    nis_dat,
    save_plot = FALSE,
    out_file = "NIS_exceedance_rate.png"
) {
  ## Visualize the fraction of observations with NIS > 95% and 99% thresholds.
  
  plot_dat <- nis_dat %>%
    group_by(variable, year) %>%
    summarise(
      n = n(),
      exceed_95_rate = mean(exceed_95, na.rm = TRUE),
      exceed_99_rate = mean(exceed_99, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(exceed_95_rate, exceed_99_rate),
      names_to = "threshold",
      values_to = "exceedance_rate"
    ) %>%
    mutate(
      threshold = recode(
        threshold,
        exceed_95_rate = "NIS > 95% threshold",
        exceed_99_rate = "NIS > 99% threshold"
      )
    )
  
  p <- ggplot(plot_dat, aes(x = year, y = exceedance_rate, color = threshold)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Year",
      y = "Exceedance rate",
      color = NULL,
      title = "NIS Exceedance Rate by Year and Variable",
      subtitle = "High exceedance indicates forecast/filter overconfidence"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 12, height = 8)
  
  list(plot = p, data = plot_dat)
}

plot_NIS_each_site <- function(
    nis_dat,
    out_dir = "NIS_each_site",
    save_plot = TRUE,
    y_log10 = FALSE
) {
  ## Save one NIS time-series plot per site.
  
  if (isTRUE(save_plot)) ensure_dir(out_dir)
  
  site_ids <- unique(as.character(nis_dat$site_id))
  plot_list <- vector("list", length(site_ids))
  names(plot_list) <- site_ids
  
  for (sid in site_ids) {
    dat_site <- nis_dat %>%
      filter(as.character(site_id) == sid)
    
    if (nrow(dat_site) == 0) next
    
    p <- ggplot(dat_site, aes(x = year, y = NIS)) +
      geom_hline(
        yintercept = stats::qchisq(0.95, df = 1),
        linetype = "dashed",
        color = "orange",
        linewidth = 0.8
      ) +
      geom_hline(
        yintercept = stats::qchisq(0.99, df = 1),
        linetype = "dashed",
        color = "red",
        linewidth = 0.8
      ) +
      geom_line(aes(group = variable), color = "gray40", linewidth = 0.7, alpha = 0.7) +
      geom_point(aes(color = exceed_99), size = 2.8, alpha = 0.85) +
      facet_wrap(~ variable, scales = "free_y") +
      scale_color_manual(
        values = c("FALSE" = "#2C7FB8", "TRUE" = "#D7301F"),
        labels = c("FALSE" = "NIS <= 99% threshold", "TRUE" = "NIS > 99% threshold")
      ) +
      labs(
        x = "Year",
        y = "NIS",
        color = NULL,
        title = paste0("NIS by Year - Site ", sid),
        subtitle = "Orange = 95% chi-square threshold; Red = 99% threshold"
      ) +
      theme_bw(base_size = 13) +
      theme(
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    if (isTRUE(y_log10)) {
      p <- p +
        scale_y_continuous(trans = "log10") +
        labs(y = "NIS, log10 scale")
    }
    
    plot_list[[sid]] <- p
    
    if (isTRUE(save_plot)) {
      ggplot2::ggsave(
        filename = file.path(out_dir, paste0("NIS_site_", safe_filename(sid), ".png")),
        plot = p,
        width = 11,
        height = 7,
        dpi = 300
      )
    }
  }
  
  plot_list
}

## ============================================================
## 6. Optional filter diagnostic plots
## ============================================================

plot_update_direction <- function(
    final_summary_all,
    soilmoist_div10 = FALSE,
    save_plot = FALSE,
    out_file = "update_direction.png"
) {
  ## Diagnose whether analysis moves forecast toward observations.
  
  plot_dat <- prepare_diagnostic_data(final_summary_all, soilmoist_div10 = soilmoist_div10) %>%
    filter(
      is.finite(obs_plot),
      is.finite(forecast_plot),
      is.finite(analysis_plot)
    )
  
  p <- ggplot(plot_dat, aes(x = forecast_resid, y = increment_plot, color = worsened)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
    geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black") +
    geom_point(alpha = 0.75, size = 2.4) +
    facet_wrap(~ variable, scales = "free") +
    scale_color_manual(
      values = c("FALSE" = "#2C7FB8", "TRUE" = "#D7301F"),
      labels = c("FALSE" = "improved / not worse", "TRUE" = "worse")
    ) +
    labs(
      x = "Forecast residual = Forecast - Observation",
      y = "Analysis increment = Analysis - Forecast",
      color = NULL,
      title = "Filter Update Direction",
      subtitle = "Dashed line is ideal: increment = -forecast residual"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 12, height = 8)
  
  list(plot = p, data = plot_dat)
}

plot_expected_uncertainty_vs_error <- function(
    final_summary_all,
    soilmoist_div10 = FALSE,
    save_plot = FALSE,
    out_file = "forecast_error_vs_expected_uncertainty.png"
) {
  ## Compare observed forecast error with sqrt(Pf + Q + R).
  
  plot_dat <- calc_scalar_nis(final_summary_all, soilmoist_div10 = soilmoist_div10)
  
  p <- ggplot(plot_dat, aes(x = expected_sd, y = abs_forecast_error)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_abline(slope = 2, intercept = 0, linetype = "dotted", color = "orange", linewidth = 0.8) +
    geom_abline(slope = 3, intercept = 0, linetype = "dotted", color = "red", linewidth = 0.8) +
    geom_point(aes(color = exceed_99, shape = worsened), alpha = 0.75, size = 2.5) +
    facet_wrap(~ variable, scales = "free") +
    scale_color_manual(
      values = c("FALSE" = "#2C7FB8", "TRUE" = "#D7301F"),
      labels = c("FALSE" = "NIS <= 99%", "TRUE" = "NIS > 99%")
    ) +
    labs(
      x = expression("Expected uncertainty: " * sqrt(P[f] + Q + R)),
      y = expression("|Forecast - Observation|"),
      color = NULL,
      shape = "Analysis worse?",
      title = "Forecast Error vs Filter Expected Uncertainty",
      subtitle = "Dashed = 1 sd, orange dotted = 2 sd, red dotted = 3 sd"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 12, height = 8)
  
  list(plot = p, data = plot_dat)
}

plot_QR_vs_abs_forecast_resid <- function(
    final_summary_all,
    soilmoist_div10 = FALSE,
    save_plot = FALSE,
    out_file = "QR_vs_abs_forecast_resid.png"
) {
  ## Plot log10(Q/R) against absolute forecast residual.
  
  plot_dat <- prepare_diagnostic_data(final_summary_all, soilmoist_div10 = soilmoist_div10) %>%
    filter(
      is.finite(obs_plot),
      is.finite(forecast_plot),
      is.finite(Q_over_R_plot),
      Q_over_R_plot > 0
    ) %>%
    mutate(log_Q_over_R = log10(Q_over_R_plot))
  
  p <- ggplot(plot_dat, aes(x = log_Q_over_R, y = abs_forecast_error)) +
    geom_point(alpha = 0.7, size = 2.3, color = "#2C7FB8") +
    geom_smooth(method = "lm", se = FALSE, color = "#D7301F", linewidth = 0.9) +
    facet_wrap(~ variable, scales = "free_y") +
    labs(
      x = expression(log[10](Q / R)),
      y = expression("|Forecast - Observation|"),
      title = "Q/R vs Absolute Forecast Residual",
      subtitle = ifelse(
        soilmoist_div10,
        "SoilMoistFrac forecast divided by 10 before computing residual",
        "Using raw forecast values"
      )
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 11, height = 7)
  
  list(plot = p, data = plot_dat)
}

plot_PfR_vs_abs_forecast_resid <- function(
    final_summary_all,
    soilmoist_div10 = FALSE,
    save_plot = FALSE,
    out_file = "PfR_vs_abs_forecast_resid.png"
) {
  ## Plot log10(Pf/R) against absolute forecast residual.
  
  plot_dat <- prepare_diagnostic_data(final_summary_all, soilmoist_div10 = soilmoist_div10) %>%
    filter(
      is.finite(obs_plot),
      is.finite(forecast_plot),
      is.finite(Pf_over_R_plot),
      Pf_over_R_plot > 0
    ) %>%
    mutate(log_Pf_over_R = log10(Pf_over_R_plot))
  
  p <- ggplot(plot_dat, aes(x = log_Pf_over_R, y = abs_forecast_error)) +
    geom_point(alpha = 0.7, size = 2.3, color = "#2C7FB8") +
    geom_smooth(method = "lm", se = FALSE, color = "#D7301F", linewidth = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~ variable, scales = "free_y") +
    labs(
      x = expression(log[10](P[f] / R)),
      y = expression("|Forecast - Observation|"),
      title = "Pf/R vs Absolute Forecast Residual",
      subtitle = ifelse(
        soilmoist_div10,
        "SoilMoistFrac forecast divided by 10 before computing residual",
        "Using raw forecast values"
      )
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 11, height = 7)
  
  list(plot = p, data = plot_dat)
}

## ============================================================
## 7. Main wrapper
## ============================================================

run_sda_diagnostic_pipeline <- function(
    sda_dir,
    out_dir,
    file_index = 1:13,
    start_year = 2012,
    state_vars = c("AbvGrndWood", "NEE", "Qle", "LAI", "SoilMoistFrac", "TotSoilCarb"),
    soilmoist_div10 = FALSE,
    prefix = "SDA_all_years",
    save_outputs = TRUE,
    make_site_nis_plots = TRUE
) {
  ## Main pipeline:
  ## 1. Summarize SDA output files.
  ## 2. Compute R2/RMSE/bias tables.
  ## 3. Compute NIS diagnostics.
  ## 4. Save tables and figures.
  ## 5. Return all objects as a list.
  
  tables_dir <- file.path(out_dir, "tables")
  figures_dir <- file.path(out_dir, "figures")
  nis_site_dir <- file.path(figures_dir, "NIS_each_site")
  
  ensure_dir(out_dir)
  ensure_dir(tables_dir)
  ensure_dir(figures_dir)
  
  ## Step 1: summarize raw SDA files.
  summary_res <- summarize_sda_directory(
    sda_dir = sda_dir,
    file_index = file_index,
    start_year = start_year,
    state_vars = state_vars,
    save_csv = save_outputs,
    out_dir = tables_dir,
    prefix = prefix
  )
  
  final_summary_all <- summary_res$final_summary_all
  diag_summary_all <- summary_res$diag_summary_all
  obs_long_all <- summary_res$obs_long_all
  
  ## Step 2: performance metrics.
  performance <- compute_performance_tables(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10
  )
  
  if (isTRUE(save_outputs)) {
    write_performance_tables(performance, tables_dir)
  }
  
  ## Step 3: NIS diagnostics.
  nis_dat <- calc_scalar_nis(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10
  )
  
  if (isTRUE(save_outputs)) {
    readr::write_csv(nis_dat, file.path(tables_dir, "NIS_diagnostics.csv"))
  }
  
  ## Step 4: figures.
  plots <- list()
  
  plots$scatter <- plot_obs_forecast_analysis(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "01_obs_forecast_analysis_scatter.png")
  )
  
  plots$R2_by_variable <- plot_R2_by_variable(
    performance = performance,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "02_R2_by_variable.png")
  )
  
  plots$R2_improvement_by_variable <- plot_R2_improvement_by_variable(
    performance = performance,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "03_R2_improvement_by_variable.png")
  )
  
  plots$R2_by_year_variable <- plot_R2_by_year_variable(
    performance = performance,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "04_R2_by_year_variable.png")
  )
  
  plots$R2_improvement_by_year_variable <- plot_R2_improvement_by_year_variable(
    performance = performance,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "05_R2_improvement_by_year_variable.png")
  )
  
  plots$NIS_by_year <- plot_NIS_by_year(
    nis_dat = nis_dat,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "06_NIS_by_year.png")
  )
  
  plots$NIS_exceedance_rate <- plot_NIS_exceedance_rate(
    nis_dat = nis_dat,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "07_NIS_exceedance_rate.png")
  )
  
  plots$update_direction <- plot_update_direction(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "08_update_direction.png")
  )
  
  plots$expected_uncertainty_vs_error <- plot_expected_uncertainty_vs_error(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "09_forecast_error_vs_expected_uncertainty.png")
  )
  
  plots$QR_vs_abs_forecast_resid <- plot_QR_vs_abs_forecast_resid(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "10_QR_vs_abs_forecast_resid.png")
  )
  
  plots$PfR_vs_abs_forecast_resid <- plot_PfR_vs_abs_forecast_resid(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "11_PfR_vs_abs_forecast_resid.png")
  )
  
  if (isTRUE(make_site_nis_plots)) {
    plots$NIS_each_site <- plot_NIS_each_site(
      nis_dat = nis_dat,
      out_dir = nis_site_dir,
      save_plot = save_outputs,
      y_log10 = FALSE
    )
  }
  
  ## Step 5: return all objects.
  list(
    final_summary_all = final_summary_all,
    diag_summary_all = diag_summary_all,
    obs_long_all = obs_long_all,
    nis_dat = nis_dat,
    performance = performance,
    plots = plots,
    yearly_results = summary_res$yearly_results,
    output_paths = list(
      out_dir = out_dir,
      tables_dir = tables_dir,
      figures_dir = figures_dir,
      nis_site_dir = nis_site_dir
    )
  )
}

## ============================================================
## 8. Example usage
## ============================================================

## Example:
## res_diag <- run_sda_diagnostic_pipeline(
##   sda_dir = "/projectnb/dietzelab/guYANG/pecan/output/param_uncalibrated",
##   out_dir = "/projectnb/dietzelab/guYANG/pecan/output/adjust/SDA_diagnostics",
##   file_index = 1:13,
##   start_year = 2012,
##   soilmoist_div10 = FALSE,
##   prefix = "param_uncalibrated",
##   save_outputs = TRUE,
##   make_site_nis_plots = TRUE
## )
##
## final_summary_all <- res_diag$final_summary_all
## diag_summary_all <- res_diag$diag_summary_all
## obs_long_all <- res_diag$obs_long_all
## nis_dat <- res_diag$nis_dat
##
## res_diag$performance$performance_by_variable
## res_diag$performance$r2_wide_by_variable
## res_diag$plots$scatter$plot
## res_diag$plots$NIS_by_year$plot
## res_diag$plots$R2_by_variable$plot

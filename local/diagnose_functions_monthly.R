## ============================================================
## SDA / EnKF Monthly Diagnostic Pipeline
## q.type == 3 vector-Q robust version
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
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  
  if (length(obs) < 2) return(NA_real_)
  
  sst <- sum((obs - mean(obs))^2)
  if (!is.finite(sst) || sst == 0) return(NA_real_)
  
  1 - sum((pred - obs)^2) / sst
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", as.character(x))
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

save_plot_if_needed <- function(plot, filename, save_plot = TRUE,
                                width = 12, height = 8, dpi = 300) {
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

make_monthly_assim_dates <- function(
    start_date = "2012-07-15",
    end_date = "2024-07-15"
) {
  seq(as.Date(start_date), as.Date(end_date), by = "1 month")
}

parse_block_dates <- function(non_null_dates) {
  block_dates <- suppressWarnings(as.Date(non_null_dates))
  
  if (all(is.na(block_dates))) {
    block_dates <- suppressWarnings(as.Date(substr(non_null_dates, 1, 10)))
  }
  
  block_dates
}

pick_block_date_by_assim_date <- function(
    block_all,
    assim_date,
    date_key = NULL,
    allow_same_month_fallback = TRUE,
    allow_single_date_fallback = TRUE
) {
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
  
  assim_date <- as.Date(assim_date)
  block_dates <- parse_block_dates(non_null_dates)
  
  matched_dates <- non_null_dates[
    !is.na(block_dates) & block_dates == assim_date
  ]
  
  if (length(matched_dates) > 0) {
    return(matched_dates[1])
  }
  
  if (isTRUE(allow_same_month_fallback)) {
    target_ym <- format(assim_date, "%Y-%m")
    
    matched_dates <- non_null_dates[
      !is.na(block_dates) & format(block_dates, "%Y-%m") == target_ym
    ]
    
    if (length(matched_dates) > 0) {
      return(matched_dates[length(matched_dates)])
    }
  }
  
  if (isTRUE(allow_single_date_fallback) && length(non_null_dates) == 1) {
    return(non_null_dates[1])
  }
  
  stop(
    "No block.list.all date found for assim_date = ", assim_date,
    "\nAvailable dates: ", paste(non_null_dates, collapse = ", ")
  )
}

make_diag_table_from_values_monthly <- function(
    values,
    value_name,
    state_vars,
    obs_vars = NULL,
    time_index,
    assim_date,
    date_key,
    site_index,
    site_id,
    prefer_obs = FALSE
) {
  ## Build diagnostic table from a vector.
  ## This avoids dimnames errors when monthly observed variables differ.
  
  assim_date <- as.Date(assim_date)
  year <- as.integer(format(assim_date, "%Y"))
  month <- as.integer(format(assim_date, "%m"))
  year_month <- format(assim_date, "%Y-%m")
  
  values <- as.numeric(values)
  n_val <- length(values)
  
  n_state <- length(state_vars)
  n_obs <- if (!is.null(obs_vars)) length(obs_vars) else NA_integer_
  
  if (n_val == 0) {
    return(data.frame())
  }
  
  if (prefer_obs && !is.null(obs_vars) && n_val == n_obs) {
    vars <- obs_vars
    
  } else if (n_val == n_state) {
    vars <- state_vars
    
  } else if (!is.null(obs_vars) && n_val == n_obs) {
    vars <- obs_vars
    
  } else if (n_val > n_state) {
    warning(
      value_name, " length = ", n_val,
      " is larger than number of state_vars = ", n_state,
      " at site ", site_id, ", date ", assim_date,
      ". Using first ", n_state, " values."
    )
    
    values <- values[seq_len(n_state)]
    vars <- state_vars
    
  } else {
    warning(
      value_name, " length = ", n_val,
      " does not match state_vars or obs_vars at site ", site_id,
      ", date ", assim_date,
      ". Using first ", n_val, " state variable names."
    )
    
    vars <- state_vars[seq_len(n_val)]
  }
  
  out <- data.frame(
    time_index = time_index,
    assim_date = assim_date,
    year = year,
    month = month,
    year_month = year_month,
    date_key = date_key,
    site_index = site_index,
    site_id = as.character(site_id),
    variable = vars,
    value = values,
    row.names = NULL
  )
  
  names(out)[names(out) == "value"] <- value_name
  out
}

## ============================================================
## 1. Extract Pf / Q / R from one block
## ============================================================

extract_site_diag_from_block_monthly <- function(
    blk,
    site_index,
    site_id,
    state_vars,
    time_index,
    assim_date,
    date_key
) {
  ## Robust monthly version.
  ##
  ## Important:
  ## - Pf is usually a full state covariance matrix.
  ## - R is observation covariance, usually only observed variables.
  ## - q.type == 3 means vector precision q, not Q covariance matrix.
  ##   Therefore:
  ##       q_mean = aq / bq
  ##       Q variance = 1 / q_mean
  
  assim_date <- as.Date(assim_date)
  n_state <- length(state_vars)
  
  H_mat <- blk$H
  
  obs_idx <- blk$constant$H
  if (is.null(obs_idx)) {
    obs_idx <- which(colSums(H_mat) != 0)
  }
  
  obs_vars <- state_vars[((obs_idx - 1) %% n_state) + 1]
  obs_vars <- obs_vars[!is.na(obs_vars)]
  obs_vars <- unique(obs_vars)
  
  ## -----------------------------
  ## 1. Pf
  ## -----------------------------
  Pf_mat <- as.matrix(blk$data$pf)
  Pf_values <- as.numeric(diag(Pf_mat))
  
  Pf_diag <- make_diag_table_from_values_monthly(
    values = Pf_values,
    value_name = "Pf",
    state_vars = state_vars,
    obs_vars = obs_vars,
    time_index = time_index,
    assim_date = assim_date,
    date_key = date_key,
    site_index = site_index,
    site_id = site_id,
    prefer_obs = FALSE
  )
  
  ## -----------------------------
  ## 2. R
  ## -----------------------------
  R_mat <- tryCatch(
    solve(as.matrix(blk$data$r)),
    error = function(e) {
      warning(
        "Failed to invert R precision for site ", site_id,
        " at date ", assim_date, ": ", conditionMessage(e)
      )
      
      n_r <- nrow(as.matrix(blk$data$r))
      matrix(NA_real_, nrow = n_r, ncol = n_r)
    }
  )
  
  R_values <- as.numeric(diag(R_mat))
  
  R_diag <- make_diag_table_from_values_monthly(
    values = R_values,
    value_name = "R",
    state_vars = state_vars,
    obs_vars = obs_vars,
    time_index = time_index,
    assim_date = assim_date,
    date_key = date_key,
    site_index = site_index,
    site_id = site_id,
    prefer_obs = TRUE
  )
  
  ## -----------------------------
  ## 3. Q
  ## -----------------------------
  aq <- blk$update$aq
  bq <- blk$update$bq
  
  if (is.null(aq) || is.null(bq)) {
    aq <- blk$data$aq
    bq <- blk$data$bq
  }
  
  q_type <- blk$constant$q.type
  
  if (q_type == 3) {
    ## Vector precision q.
    aq_vec <- as.numeric(aq)
    bq_vec <- as.numeric(bq)
    
    if (length(bq_vec) == 1) {
      q_mean <- aq_vec / bq_vec
    } else if (length(aq_vec) == length(bq_vec)) {
      q_mean <- aq_vec / bq_vec
    } else {
      stop(
        "q.type == 3 but aq and bq lengths do not match at site ", site_id,
        ", date ", assim_date,
        ". length(aq) = ", length(aq_vec),
        ", length(bq) = ", length(bq_vec)
      )
    }
    
    Q_values <- 1 / q_mean
    
    Q_diag <- make_diag_table_from_values_monthly(
      values = Q_values,
      value_name = "Q",
      state_vars = state_vars,
      obs_vars = obs_vars,
      time_index = time_index,
      assim_date = assim_date,
      date_key = date_key,
      site_index = site_index,
      site_id = site_id,
      prefer_obs = TRUE
    )
    
  } else if (q_type == 4) {
    ## Wishart compatibility.
    Q_mat <- as.matrix(aq) / as.numeric(bq)
    Q_values <- as.numeric(diag(Q_mat))
    
    Q_diag <- make_diag_table_from_values_monthly(
      values = Q_values,
      value_name = "Q",
      state_vars = state_vars,
      obs_vars = obs_vars,
      time_index = time_index,
      assim_date = assim_date,
      date_key = date_key,
      site_index = site_index,
      site_id = site_id,
      prefer_obs = FALSE
    )
    
  } else {
    stop("Unsupported q.type: ", q_type)
  }
  
  ## -----------------------------
  ## 4. Merge
  ## -----------------------------
  by_cols <- c(
    "time_index", "assim_date", "year", "month", "year_month",
    "date_key", "site_index", "site_id", "variable"
  )
  
  diag_summary <- Pf_diag %>%
    left_join(Q_diag, by = by_cols) %>%
    left_join(R_diag, by = by_cols) %>%
    mutate(
      Q_over_R = Q / R,
      Pf_over_R = Pf / R
    )
  
  list(
    site_id = as.character(site_id),
    obs_vars = obs_vars,
    H = H_mat,
    Q = Q_diag,
    R = R_diag,
    Pf = Pf_diag,
    diag_summary = diag_summary
  )
}

## ============================================================
## 2. Summarize one SDA output
## ============================================================

summarize_sda_output_monthly <- function(
    sda.outputs,
    time_index,
    assim_date,
    date_key = NULL,
    state_vars = c("AbvGrndWood", "NEE", "Qle", "LAI", "SoilMoistFrac", "TotSoilCarb")
) {
  assim_date <- as.Date(assim_date)
  year <- as.integer(format(assim_date, "%Y"))
  month <- as.integer(format(assim_date, "%m"))
  year_month <- format(assim_date, "%Y-%m")
  
  forecast <- as.data.frame(sda.outputs[["forecast"]])
  analysis <- as.data.frame(sda.outputs[["analysis"]])
  obs_mean <- sda.outputs[["obs.mean"]]
  block_all <- sda.outputs[["enkf.params"]][["block.list.all"]]
  
  date_key <- pick_block_date_by_assim_date(
    block_all = block_all,
    assim_date = assim_date,
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
  
  col_map <- data.frame(
    time_index = time_index,
    assim_date = assim_date,
    year = year,
    month = month,
    year_month = year_month,
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
  
  obs_long <- bind_rows(lapply(names(obs_mean), function(sid) {
    x <- obs_mean[[sid]]
    
    if (is.null(x) || ncol(x) == 0) {
      return(data.frame())
    }
    
    data.frame(
      time_index = time_index,
      assim_date = assim_date,
      year = year,
      month = month,
      year_month = year_month,
      date_key = date_key,
      site_id = as.character(sid),
      variable = colnames(x),
      obs_value = as.numeric(x[1, ]),
      row.names = NULL
    )
  }))
  
  matrices_by_site <- lapply(seq_along(blk_list), function(i) {
    extract_site_diag_from_block_monthly(
      blk = blk_list[[i]],
      site_index = i,
      site_id = site_ids[i],
      state_vars = state_vars,
      time_index = time_index,
      assim_date = assim_date,
      date_key = date_key
    )
  })
  
  names(matrices_by_site) <- site_ids
  
  diag_summary <- bind_rows(lapply(matrices_by_site, `[[`, "diag_summary"))
  
  final_summary <- state_mean %>%
    left_join(
      obs_long,
      by = c(
        "time_index", "assim_date", "year", "month", "year_month",
        "date_key", "site_id", "variable"
      )
    ) %>%
    left_join(
      diag_summary,
      by = c(
        "time_index", "assim_date", "year", "month", "year_month",
        "date_key", "site_index", "site_id", "variable"
      )
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
    time_index = time_index,
    assim_date = assim_date,
    year = year,
    month = month,
    year_month = year_month,
    date_key = date_key,
    nens = nens,
    nsite = nsite,
    final_summary = final_summary,
    diag_summary = diag_summary,
    obs_long = obs_long,
    matrices_by_site = matrices_by_site
  )
}

## ============================================================
## 3. Summarize directory
## ============================================================

summarize_sda_directory_monthly <- function(
    sda_dir,
    file_index,
    assim_dates,
    state_vars = c("AbvGrndWood", "NEE", "Qle", "LAI", "SoilMoistFrac", "TotSoilCarb"),
    save_csv = TRUE,
    out_dir = sda_dir,
    prefix = "SDA_all_months"
) {
  assim_dates <- as.Date(assim_dates)
  
  if (length(assim_dates) < max(file_index)) {
    stop(
      "Length of assim_dates is smaller than max(file_index). ",
      "length(assim_dates) = ", length(assim_dates),
      ", max(file_index) = ", max(file_index)
    )
  }
  
  all_results <- list()
  
  for (i in file_index) {
    assim_date_i <- assim_dates[i]
    file_i <- file.path(sda_dir, paste0("sda.output", i, ".Rdata"))
    
    if (!file.exists(file_i)) {
      warning("File not found: ", file_i)
      next
    }
    
    message("Reading: ", file_i, " | assim_date = ", assim_date_i)
    
    env <- new.env(parent = emptyenv())
    load(file_i, envir = env)
    
    if (!exists("sda.outputs", envir = env)) {
      warning("No object named sda.outputs in: ", file_i)
      next
    }
    
    res_i <- summarize_sda_output_monthly(
      sda.outputs = env$sda.outputs,
      time_index = i,
      assim_date = assim_date_i,
      date_key = NULL,
      state_vars = state_vars
    )
    
    all_results[[as.character(assim_date_i)]] <- res_i
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
    monthly_results = all_results
  )
}

## ============================================================
## 4. Prepare diagnostic data
## ============================================================

prepare_diagnostic_data_monthly <- function(final_summary_all, soilmoist_div10 = FALSE) {
  final_summary_all %>%
    mutate(
      assim_date = as.Date(assim_date),
      
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
## 5. Performance tables
## ============================================================

compute_performance_tables_monthly <- function(final_summary_all, soilmoist_div10 = FALSE) {
  dat <- prepare_diagnostic_data_monthly(
    final_summary_all,
    soilmoist_div10 = soilmoist_div10
  )
  
  pred_long <- dat %>%
    filter(is.finite(obs_plot)) %>%
    select(
      time_index, assim_date, year, month, year_month,
      site_id, variable,
      obs_plot, forecast_plot, analysis_plot
    ) %>%
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
      type = factor(type, levels = c("forecast", "analysis")),
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
  
  performance_by_date_variable <- pred_long %>%
    group_by(time_index, assim_date, year, month, year_month, variable, type) %>%
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
  
  performance_by_calendar_month_variable <- pred_long %>%
    group_by(month, variable, type) %>%
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
  
  r2_wide_by_date_variable <- performance_by_date_variable %>%
    select(time_index, assim_date, year, month, year_month,
           variable, type, n, R2_sse, RMSE, bias, MAE) %>%
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
  
  r2_wide_by_calendar_month_variable <- performance_by_calendar_month_variable %>%
    select(month, variable, type, n, R2_sse, RMSE, bias, MAE) %>%
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
    performance_by_date_variable = performance_by_date_variable,
    performance_by_year_variable = performance_by_year_variable,
    performance_by_calendar_month_variable = performance_by_calendar_month_variable,
    performance_overall = performance_overall,
    r2_wide_by_variable = r2_wide_by_variable,
    r2_wide_by_date_variable = r2_wide_by_date_variable,
    r2_wide_by_year_variable = r2_wide_by_year_variable,
    r2_wide_by_calendar_month_variable = r2_wide_by_calendar_month_variable
  )
}

write_performance_tables_monthly <- function(performance, out_dir) {
  ensure_dir(out_dir)
  
  readr::write_csv(performance$performance_by_variable,
                   file.path(out_dir, "performance_by_variable.csv"))
  
  readr::write_csv(performance$performance_by_date_variable,
                   file.path(out_dir, "performance_by_date_variable.csv"))
  
  readr::write_csv(performance$performance_by_year_variable,
                   file.path(out_dir, "performance_by_year_variable.csv"))
  
  readr::write_csv(performance$performance_by_calendar_month_variable,
                   file.path(out_dir, "performance_by_calendar_month_variable.csv"))
  
  readr::write_csv(performance$performance_overall,
                   file.path(out_dir, "performance_overall.csv"))
  
  readr::write_csv(performance$r2_wide_by_variable,
                   file.path(out_dir, "R2_wide_by_variable.csv"))
  
  readr::write_csv(performance$r2_wide_by_date_variable,
                   file.path(out_dir, "R2_wide_by_date_variable.csv"))
  
  readr::write_csv(performance$r2_wide_by_year_variable,
                   file.path(out_dir, "R2_wide_by_year_variable.csv"))
  
  readr::write_csv(performance$r2_wide_by_calendar_month_variable,
                   file.path(out_dir, "R2_wide_by_calendar_month_variable.csv"))
  
  invisible(performance)
}

## ============================================================
## 6. NIS
## ============================================================

calc_scalar_nis_monthly <- function(final_summary_all, soilmoist_div10 = FALSE) {
  prepare_diagnostic_data_monthly(
    final_summary_all,
    soilmoist_div10 = soilmoist_div10
  ) %>%
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
      assim_date = as.Date(assim_date),
      innovation = forecast_plot - obs_plot,
      S = Pf_plot + Q_plot + R_plot,
      expected_sd = sqrt(S),
      NIS = innovation^2 / S,
      abs_normalized_innovation = sqrt(NIS),
      exceed_95 = NIS > stats::qchisq(0.95, df = 1),
      exceed_99 = NIS > stats::qchisq(0.99, df = 1)
    )
}

## ============================================================
## 7. Plot functions
## ============================================================

plot_obs_forecast_analysis_monthly <- function(
    final_summary_all,
    soilmoist_div10 = FALSE,
    save_plot = FALSE,
    out_file = "obs_forecast_analysis_scatter_monthly.png"
) {
  performance <- compute_performance_tables_monthly(
    final_summary_all,
    soilmoist_div10 = soilmoist_div10
  )
  
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
    geom_point(size = 2.0, alpha = 0.55) +
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
    scale_color_manual(values = c("forecast" = "#2C7FB8", "analysis" = "#D7301F")) +
    labs(
      x = "Observation",
      y = "Prediction",
      color = NULL,
      title = "Observed vs Forecast / Analysis",
      subtitle = "Monthly diagnostics"
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

plot_R2_by_variable_monthly <- function(
    performance,
    save_plot = FALSE,
    out_file = "R2_by_variable_monthly.png"
) {
  plot_dat <- performance$performance_by_variable %>%
    mutate(type = factor(type, levels = c("forecast", "analysis")))
  
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
      subtitle = "All monthly assimilation dates pooled together"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 10, height = 6)
  
  list(plot = p, data = plot_dat)
}

plot_R2_improvement_by_variable_monthly <- function(
    performance,
    save_plot = FALSE,
    out_file = "R2_improvement_by_variable_monthly.png"
) {
  plot_dat <- performance$r2_wide_by_variable %>%
    mutate(improved = R2_improvement > 0)
  
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

plot_R2_by_date_variable_monthly <- function(
    performance,
    save_plot = FALSE,
    out_file = "R2_by_date_variable_monthly.png"
) {
  plot_dat <- performance$performance_by_date_variable %>%
    mutate(
      assim_date = as.Date(assim_date),
      type = factor(type, levels = c("forecast", "analysis"))
    )
  
  p <- ggplot(plot_dat, aes(x = assim_date, y = R2_sse, color = type, linetype = type)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    geom_line(linewidth = 0.75, na.rm = TRUE) +
    geom_point(size = 1.7, na.rm = TRUE) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_color_manual(values = c("forecast" = "#2C7FB8", "analysis" = "#D7301F")) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      x = "Assimilation date",
      y = "SSE-based R²",
      color = NULL,
      linetype = NULL,
      title = "R² by Monthly Assimilation Date and Variable"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 13, height = 8)
  
  list(plot = p, data = plot_dat)
}

plot_R2_improvement_by_date_variable_monthly <- function(
    performance,
    save_plot = FALSE,
    out_file = "R2_improvement_by_date_variable_monthly.png"
) {
  plot_dat <- performance$r2_wide_by_date_variable %>%
    mutate(
      assim_date = as.Date(assim_date),
      improved = R2_improvement > 0
    )
  
  p <- ggplot(plot_dat, aes(x = assim_date, y = R2_improvement, fill = improved)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_col(alpha = 0.9, width = 25) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_fill_manual(
      values = c("FALSE" = "#D7301F", "TRUE" = "#2C7FB8"),
      labels = c("FALSE" = "analysis worse", "TRUE" = "analysis better")
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      x = "Assimilation date",
      y = "R² improvement",
      fill = NULL,
      title = "R² Improvement by Monthly Assimilation Date and Variable",
      subtitle = "Positive values mean analysis improved over forecast"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 13, height = 8)
  
  list(plot = p, data = plot_dat)
}

plot_R2_by_calendar_month_variable_monthly <- function(
    performance,
    save_plot = FALSE,
    out_file = "R2_by_calendar_month_variable.png"
) {
  plot_dat <- performance$performance_by_calendar_month_variable %>%
    mutate(type = factor(type, levels = c("forecast", "analysis")))
  
  p <- ggplot(plot_dat, aes(x = month, y = R2_sse, color = type, linetype = type)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 2.2, na.rm = TRUE) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_color_manual(values = c("forecast" = "#2C7FB8", "analysis" = "#D7301F")) +
    scale_x_continuous(breaks = 1:12) +
    labs(
      x = "Calendar month",
      y = "SSE-based R²",
      color = NULL,
      linetype = NULL,
      title = "Seasonal R² by Calendar Month and Variable"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 13, height = 8)
  
  list(plot = p, data = plot_dat)
}

plot_NIS_by_date_monthly <- function(
    nis_dat,
    save_plot = FALSE,
    out_file = "NIS_by_date_monthly.png",
    y_log10 = FALSE
) {
  plot_dat <- nis_dat %>%
    mutate(assim_date = as.Date(assim_date)) %>%
    group_by(assim_date, year, month, year_month, variable) %>%
    summarise(
      n = n(),
      NIS_median = median(NIS, na.rm = TRUE),
      NIS_p25 = quantile(NIS, 0.25, na.rm = TRUE),
      NIS_p75 = quantile(NIS, 0.75, na.rm = TRUE),
      NIS_mean = mean(NIS, na.rm = TRUE),
      exceed_95_rate = mean(exceed_95, na.rm = TRUE),
      exceed_99_rate = mean(exceed_99, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(plot_dat, aes(x = assim_date, y = NIS_median)) +
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
    geom_ribbon(aes(ymin = NIS_p25, ymax = NIS_p75), alpha = 0.25) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      x = "Assimilation date",
      y = "Median NIS with IQR",
      title = "Normalized Innovation Squared by Monthly Assimilation Date",
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
      labs(y = "Median NIS with IQR, log10 scale")
  }
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 13, height = 8)
  
  list(plot = p, data = plot_dat)
}

plot_NIS_exceedance_rate_monthly <- function(
    nis_dat,
    save_plot = FALSE,
    out_file = "NIS_exceedance_rate_monthly.png"
) {
  plot_dat <- nis_dat %>%
    mutate(assim_date = as.Date(assim_date)) %>%
    group_by(assim_date, year, month, year_month, variable) %>%
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
  
  p <- ggplot(plot_dat, aes(x = assim_date, y = exceedance_rate, color = threshold)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.7) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      x = "Assimilation date",
      y = "Exceedance rate",
      color = NULL,
      title = "NIS Exceedance Rate by Monthly Assimilation Date and Variable",
      subtitle = "High exceedance indicates forecast/filter overconfidence"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 13, height = 8)
  
  list(plot = p, data = plot_dat)
}

plot_NIS_each_site_monthly <- function(
    nis_dat,
    out_dir = "NIS_each_site_monthly",
    save_plot = TRUE,
    y_log10 = FALSE
) {
  if (isTRUE(save_plot)) ensure_dir(out_dir)
  
  site_ids <- unique(as.character(nis_dat$site_id))
  plot_list <- vector("list", length(site_ids))
  names(plot_list) <- site_ids
  
  for (sid in site_ids) {
    dat_site <- nis_dat %>%
      mutate(assim_date = as.Date(assim_date)) %>%
      filter(as.character(site_id) == sid)
    
    if (nrow(dat_site) == 0) next
    
    p <- ggplot(dat_site, aes(x = assim_date, y = NIS)) +
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
      geom_point(aes(color = exceed_99), size = 2.2, alpha = 0.85) +
      facet_wrap(~ variable, scales = "free_y") +
      scale_color_manual(
        values = c("FALSE" = "#2C7FB8", "TRUE" = "#D7301F"),
        labels = c("FALSE" = "NIS <= 99% threshold", "TRUE" = "NIS > 99% threshold")
      ) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      labs(
        x = "Assimilation date",
        y = "NIS",
        color = NULL,
        title = paste0("NIS by Monthly Assimilation Date - Site ", sid),
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

plot_update_direction_monthly <- function(
    final_summary_all,
    soilmoist_div10 = FALSE,
    save_plot = FALSE,
    out_file = "update_direction_monthly.png"
) {
  plot_dat <- prepare_diagnostic_data_monthly(
    final_summary_all,
    soilmoist_div10 = soilmoist_div10
  ) %>%
    filter(
      is.finite(obs_plot),
      is.finite(forecast_plot),
      is.finite(analysis_plot)
    )
  
  p <- ggplot(plot_dat, aes(x = forecast_resid, y = increment_plot, color = worsened)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
    geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black") +
    geom_point(alpha = 0.65, size = 2.1) +
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

plot_expected_uncertainty_vs_error_monthly <- function(
    final_summary_all,
    soilmoist_div10 = FALSE,
    save_plot = FALSE,
    out_file = "forecast_error_vs_expected_uncertainty_monthly.png"
) {
  plot_dat <- calc_scalar_nis_monthly(
    final_summary_all,
    soilmoist_div10 = soilmoist_div10
  )
  
  p <- ggplot(plot_dat, aes(x = expected_sd, y = abs_forecast_error)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "black", linewidth = 0.8) +
    geom_abline(slope = 2, intercept = 0, linetype = "dotted",
                color = "orange", linewidth = 0.8) +
    geom_abline(slope = 3, intercept = 0, linetype = "dotted",
                color = "red", linewidth = 0.8) +
    geom_point(aes(color = exceed_99, shape = worsened), alpha = 0.65, size = 2.2) +
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
      title = "Forecast Error vs Filter Expected Uncertainty"
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

plot_QR_vs_abs_forecast_resid_monthly <- function(
    final_summary_all,
    soilmoist_div10 = FALSE,
    save_plot = FALSE,
    out_file = "QR_vs_abs_forecast_resid_monthly.png"
) {
  plot_dat <- prepare_diagnostic_data_monthly(
    final_summary_all,
    soilmoist_div10 = soilmoist_div10
  ) %>%
    filter(
      is.finite(obs_plot),
      is.finite(forecast_plot),
      is.finite(Q_over_R_plot),
      Q_over_R_plot > 0
    ) %>%
    mutate(log_Q_over_R = log10(Q_over_R_plot))
  
  p <- ggplot(plot_dat, aes(x = log_Q_over_R, y = abs_forecast_error)) +
    geom_point(alpha = 0.65, size = 2.1, color = "#2C7FB8") +
    geom_smooth(method = "lm", se = FALSE, color = "#D7301F", linewidth = 0.9) +
    facet_wrap(~ variable, scales = "free_y") +
    labs(
      x = expression(log[10](Q / R)),
      y = expression("|Forecast - Observation|"),
      title = "Q/R vs Absolute Forecast Residual"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  save_plot_if_needed(p, out_file, save_plot = save_plot, width = 11, height = 7)
  
  list(plot = p, data = plot_dat)
}

plot_PfR_vs_abs_forecast_resid_monthly <- function(
    final_summary_all,
    soilmoist_div10 = FALSE,
    save_plot = FALSE,
    out_file = "PfR_vs_abs_forecast_resid_monthly.png"
) {
  plot_dat <- prepare_diagnostic_data_monthly(
    final_summary_all,
    soilmoist_div10 = soilmoist_div10
  ) %>%
    filter(
      is.finite(obs_plot),
      is.finite(forecast_plot),
      is.finite(Pf_over_R_plot),
      Pf_over_R_plot > 0
    ) %>%
    mutate(log_Pf_over_R = log10(Pf_over_R_plot))
  
  p <- ggplot(plot_dat, aes(x = log_Pf_over_R, y = abs_forecast_error)) +
    geom_point(alpha = 0.65, size = 2.1, color = "#2C7FB8") +
    geom_smooth(method = "lm", se = FALSE, color = "#D7301F", linewidth = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~ variable, scales = "free_y") +
    labs(
      x = expression(log[10](P[f] / R)),
      y = expression("|Forecast - Observation|"),
      title = "Pf/R vs Absolute Forecast Residual"
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
## 8. Main wrapper
## ============================================================

run_sda_diagnostic_pipeline_monthly <- function(
    sda_dir,
    out_dir,
    file_index = NULL,
    assim_dates = NULL,
    start_date = "2012-07-15",
    end_date = "2024-07-15",
    state_vars = c("AbvGrndWood", "NEE", "Qle", "LAI", "SoilMoistFrac", "TotSoilCarb"),
    soilmoist_div10 = FALSE,
    prefix = "SDA_all_months",
    save_outputs = TRUE,
    make_site_nis_plots = TRUE
) {
  if (is.null(assim_dates)) {
    assim_dates <- make_monthly_assim_dates(
      start_date = start_date,
      end_date = end_date
    )
  }
  
  assim_dates <- as.Date(assim_dates)
  
  if (is.null(file_index)) {
    file_index <- seq_along(assim_dates)
  }
  
  tables_dir <- file.path(out_dir, "tables")
  figures_dir <- file.path(out_dir, "figures")
  nis_site_dir <- file.path(figures_dir, "NIS_each_site")
  
  ensure_dir(out_dir)
  ensure_dir(tables_dir)
  ensure_dir(figures_dir)
  
  summary_res <- summarize_sda_directory_monthly(
    sda_dir = sda_dir,
    file_index = file_index,
    assim_dates = assim_dates,
    state_vars = state_vars,
    save_csv = save_outputs,
    out_dir = tables_dir,
    prefix = prefix
  )
  
  final_summary_all <- summary_res$final_summary_all
  diag_summary_all <- summary_res$diag_summary_all
  obs_long_all <- summary_res$obs_long_all
  
  performance <- compute_performance_tables_monthly(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10
  )
  
  if (isTRUE(save_outputs)) {
    write_performance_tables_monthly(performance, tables_dir)
  }
  
  nis_dat <- calc_scalar_nis_monthly(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10
  )
  
  if (isTRUE(save_outputs)) {
    readr::write_csv(nis_dat, file.path(tables_dir, "NIS_diagnostics.csv"))
  }
  
  plots <- list()
  
  plots$scatter <- plot_obs_forecast_analysis_monthly(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "01_obs_forecast_analysis_scatter.png")
  )
  
  plots$R2_by_variable <- plot_R2_by_variable_monthly(
    performance = performance,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "02_R2_by_variable.png")
  )
  
  plots$R2_improvement_by_variable <- plot_R2_improvement_by_variable_monthly(
    performance = performance,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "03_R2_improvement_by_variable.png")
  )
  
  plots$R2_by_date_variable <- plot_R2_by_date_variable_monthly(
    performance = performance,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "04_R2_by_date_variable.png")
  )
  
  plots$R2_improvement_by_date_variable <- plot_R2_improvement_by_date_variable_monthly(
    performance = performance,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "05_R2_improvement_by_date_variable.png")
  )
  
  plots$R2_by_calendar_month_variable <- plot_R2_by_calendar_month_variable_monthly(
    performance = performance,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "06_R2_by_calendar_month_variable.png")
  )
  
  plots$NIS_by_date <- plot_NIS_by_date_monthly(
    nis_dat = nis_dat,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "07_NIS_by_date.png")
  )
  
  plots$NIS_exceedance_rate <- plot_NIS_exceedance_rate_monthly(
    nis_dat = nis_dat,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "08_NIS_exceedance_rate.png")
  )
  
  plots$update_direction <- plot_update_direction_monthly(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "09_update_direction.png")
  )
  
  plots$expected_uncertainty_vs_error <- plot_expected_uncertainty_vs_error_monthly(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "10_forecast_error_vs_expected_uncertainty.png")
  )
  
  plots$QR_vs_abs_forecast_resid <- plot_QR_vs_abs_forecast_resid_monthly(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "11_QR_vs_abs_forecast_resid.png")
  )
  
  plots$PfR_vs_abs_forecast_resid <- plot_PfR_vs_abs_forecast_resid_monthly(
    final_summary_all = final_summary_all,
    soilmoist_div10 = soilmoist_div10,
    save_plot = save_outputs,
    out_file = file.path(figures_dir, "12_PfR_vs_abs_forecast_resid.png")
  )
  
  if (isTRUE(make_site_nis_plots)) {
    plots$NIS_each_site <- plot_NIS_each_site_monthly(
      nis_dat = nis_dat,
      out_dir = nis_site_dir,
      save_plot = save_outputs,
      y_log10 = FALSE
    )
  }
  
  list(
    final_summary_all = final_summary_all,
    diag_summary_all = diag_summary_all,
    obs_long_all = obs_long_all,
    nis_dat = nis_dat,
    performance = performance,
    plots = plots,
    monthly_results = summary_res$monthly_results,
    output_paths = list(
      out_dir = out_dir,
      tables_dir = tables_dir,
      figures_dir = figures_dir,
      nis_site_dir = nis_site_dir
    )
  )
}
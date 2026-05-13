##' @title Read restart function for SDA with SIPNET
##' 
##' @author Ann Raiho \email{araiho@@nd.edu}
##' 
##' @inheritParams PEcAn.ModelName::read_restart.ModelName
##' 
##' @description Read Restart for SIPNET
##' 
##' @return X.vec      vector of forecasts
##' @export
read_restart.SIPNET <- function(outdir, runid, stop.time, settings, var.names, params) {
  
  prior.sla <- params[[which(!names(params) %in% c("soil", "soil_SDA", "restart"))[1]]]$SLA
  
  forecast <- list()
  params$restart <- c()
  
  state.vars <- c(
    "SWE",
    "SoilMoist",
    "SoilMoistFrac",
    "AbvGrndWood",
    "NEE",
    "Qle",
    "TotSoilCarb",
    "LAI",
    "litter_carbon_content",
    "fine_root_carbon_content",
    "coarse_root_carbon_content",
    "litter_mass_content_of_water"
  )
  
  params$restart <- rep(NA, length(setdiff(state.vars, var.names)))
  names(params$restart) <- setdiff(state.vars, var.names)
  
  # ------------------------------------------------------------------
  # 1. Find and open the interval NetCDF.
  #    Example:
  #    20120716T000000_to_20120815T235959.nc
  # ------------------------------------------------------------------
  
  nc_path <- find_sipnet_interval_nc(
    outdir = outdir,
    runid = runid,
    stop.time = stop.time
  )
  
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  # ------------------------------------------------------------------
  # 2. Read variables from interval NetCDF.
  # ------------------------------------------------------------------
  
  vars_to_read <- unique(c(
    state.vars,
    "time_bounds",
    "GWBI",
    "TotLivBiom",
    "leaf_carbon_content"
  ))
  
  ens <- list()
  
  for (v in vars_to_read) {
    ens[[v]] <- get_nc_var_or_null(nc, v)
  }
  
  # Keep behavior close to read.output: ensure missing vars exist as NULL.
  for (v in vars_to_read) {
    if (!v %in% names(ens)) {
      ens[[v]] <- NULL
    }
  }
  
  # ------------------------------------------------------------------
  # 3. Check key variables needed for restart.
  # ------------------------------------------------------------------
  
  must_have <- c(
    "AbvGrndWood",
    "fine_root_carbon_content",
    "coarse_root_carbon_content",
    "LAI",
    "litter_carbon_content",
    "SoilMoist",
    "SoilMoistFrac",
    "SWE",
    "TotSoilCarb"
  )
  
  missing <- must_have[sapply(ens[must_have], is.null)]
  
  if (length(missing) > 0) {
    available <- names(nc$var)
    
    stop(
      "Missing variables in interval NetCDF for runid ", runid, ": ",
      paste(missing, collapse = ", "),
      "\nNetCDF file: ", nc_path,
      "\nAvailable variables: ", paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  
  # ------------------------------------------------------------------
  # 4. Determine final time index from NetCDF time, not from AGB.
  # ------------------------------------------------------------------
  
  last <- get_last_time_index_matching_stop(nc, stop.time)
  
  # Optional sanity check.
  if (is.na(last) || last <= 0) {
    stop("Invalid final time index `last`: ", last, call. = FALSE)
  }
  
  # ------------------------------------------------------------------
  # 5. Pool/state variables use final timestep.
  #    Flux variables use interval mean.
  # ------------------------------------------------------------------
  
  AbvGrndWood_last <- get_state_at_time(ens, "AbvGrndWood", last)
  fine_root_last   <- get_state_at_time(ens, "fine_root_carbon_content", last)
  coarse_root_last <- get_state_at_time(ens, "coarse_root_carbon_content", last)
  
  wood_total_C <- AbvGrndWood_last + fine_root_last + coarse_root_last
  
  if (is.na(wood_total_C) || wood_total_C <= 0) {
    wood_total_C <- 1e-04
  }
  
  if ("AbvGrndWood" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- PEcAn.utils::ud_convert(
      AbvGrndWood_last,
      "kg/m^2",
      "Mg/ha"
    )
    names(forecast[[length(forecast)]]) <- c("AbvGrndWood")
    
    params$restart["abvGrndWoodFrac"] <- AbvGrndWood_last / wood_total_C
    params$restart["coarseRootFrac"]  <- coarse_root_last / wood_total_C
    params$restart["fineRootFrac"]    <- fine_root_last / wood_total_C
    
  } else {
    
    params$restart["AbvGrndWood"] <- PEcAn.utils::ud_convert(
      AbvGrndWood_last,
      "kg/m^2",
      "g/m^2"
    )
    
    params$restart["abvGrndWoodFrac"] <- AbvGrndWood_last / wood_total_C
    params$restart["coarseRootFrac"]  <- coarse_root_last / wood_total_C
    params$restart["fineRootFrac"]    <- fine_root_last / wood_total_C
  }
  
  # GWBI is a flux-like variable, so use interval mean.
  if ("GWBI" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- PEcAn.utils::ud_convert(
      get_interval_mean(ens, "GWBI"),
      "kg/m^2/s",
      "Mg/ha/yr"
    )
    
    names(forecast[[length(forecast)]]) <- c("GWBI")
  }
  
  # NEE is flux-like, so use interval mean.
  if ("NEE" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- nee_model_to_obs(
      get_interval_mean(ens, "NEE")
    )
    
    names(forecast[[length(forecast)]]) <- c("NEE")
  }
  
  # Qle / LE is flux-like, so use interval mean.
  if ("Qle" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- get_interval_mean(ens, "Qle")
    
    names(forecast[[length(forecast)]]) <- c("Qle")
  }
  
  if ("leaf_carbon_content" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- get_state_at_time(
      ens,
      "leaf_carbon_content",
      last
    )
    
    names(forecast[[length(forecast)]]) <- c("LeafC")
  }
  
  if ("LAI" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- get_state_at_time(ens, "LAI", last)
    
    names(forecast[[length(forecast)]]) <- c("LAI")
    
  } else {
    
    params$restart["LAI"] <- get_state_at_time(ens, "LAI", last)
  }
  
  if ("litter_carbon_content" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- get_state_at_time(
      ens,
      "litter_carbon_content",
      last
    )
    
    names(forecast[[length(forecast)]]) <- c("litter_carbon_content")
    
  } else {
    
    params$restart["litter_carbon_content"] <- PEcAn.utils::ud_convert(
      get_state_at_time(ens, "litter_carbon_content", last),
      "kg m-2",
      "g m-2"
    )
  }
  
  if ("litter_mass_content_of_water" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- get_state_at_time(
      ens,
      "litter_mass_content_of_water",
      last
    )
    
    names(forecast[[length(forecast)]]) <- c("litter_mass_content_of_water")
    
  } else {
    
    if (!is.null(ens$litter_mass_content_of_water)) {
      params$restart["litter_mass_content_of_water"] <- get_state_at_time(
        ens,
        "litter_mass_content_of_water",
        last
      )
    }
  }
  
  if ("SoilMoist" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- get_state_at_time(
      ens,
      "SoilMoist",
      last
    )
    
    names(forecast[[length(forecast)]]) <- c("SoilMoist")
    
  } else {
    
    params$restart["SoilMoist"] <- get_state_at_time(ens, "SoilMoist", last)
  }
  
  if ("SoilMoistFrac" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- get_state_at_time(
      ens,
      "SoilMoistFrac",
      last
    ) * 100
    
    names(forecast[[length(forecast)]]) <- c("SoilMoistFrac")
    
  } else {
    
    params$restart["SoilMoistFrac"] <- get_state_at_time(
      ens,
      "SoilMoistFrac",
      last
    )
  }
  
  if ("SWE" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- get_state_at_time(ens, "SWE", last)
    
    names(forecast[[length(forecast)]]) <- c("SWE")
    
  } else {
    
    params$restart["SWE"] <- get_state_at_time(ens, "SWE", last) / 10
  }
  
  if ("TotLivBiom" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- PEcAn.utils::ud_convert(
      get_state_at_time(ens, "TotLivBiom", last),
      "kg/m^2",
      "Mg/ha"
    )
    
    names(forecast[[length(forecast)]]) <- c("TotLivBiom")
  }
  
  if ("TotSoilCarb" %in% var.names) {
    
    forecast[[length(forecast) + 1]] <- get_state_at_time(
      ens,
      "TotSoilCarb",
      last
    )
    
    names(forecast[[length(forecast)]]) <- c("TotSoilCarb")
    
  } else {
    
    params$restart["TotSoilCarb"] <- PEcAn.utils::ud_convert(
      get_state_at_time(ens, "TotSoilCarb", last),
      "kg m-2",
      "g m-2"
    )
  }
  
  params$restart <- stats::na.omit(params$restart)
  
  X_tmp <- list(
    X = unlist(forecast),
    params = params
  )
  
  return(X_tmp)
}

parse_sipnet_datetime <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(as.POSIXct(x, tz = "UTC"))
  }
  
  x <- as.character(x)
  
  out <- lubridate::ymd_hms(x, truncated = 3, tz = "UTC", quiet = TRUE)
  
  if (is.na(out)) {
    out <- lubridate::ymd(x, tz = "UTC", quiet = TRUE)
  }
  
  if (is.na(out)) {
    stop("Cannot parse datetime: ", x, call. = FALSE)
  }
  
  as.POSIXct(out, tz = "UTC")
}


sipnet_interval_tag <- function(x) {
  x <- parse_sipnet_datetime(x)
  format(x, "%Y%m%dT%H%M%S", tz = "UTC")
}


find_sipnet_interval_nc <- function(outdir, runid, stop.time) {
  
  run_outdir <- file.path(outdir, runid)
  
  if (!dir.exists(run_outdir)) {
    stop("Run output directory does not exist: ", run_outdir, call. = FALSE)
  }
  
  stop_tag <- sipnet_interval_tag(stop.time)
  
  files <- list.files(
    run_outdir,
    pattern = paste0("_to_", stop_tag, "\\.nc$"),
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    available_nc <- list.files(run_outdir, pattern = "\\.nc$", full.names = FALSE)
    
    stop(
      "Interval NetCDF not found for runid ", runid,
      "\nstop.time = ", as.character(stop.time),
      "\nExpected pattern: *_to_", stop_tag, ".nc",
      "\nSearch folder: ", run_outdir,
      "\nAvailable nc files: ", paste(available_nc, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (length(files) > 1) {
    stop(
      "Multiple interval NetCDF files found for runid ", runid,
      " and stop.time = ", as.character(stop.time),
      ":\n", paste(files, collapse = "\n"),
      call. = FALSE
    )
  }
  
  files[[1]]
}


get_nc_time_values <- function(nc) {
  
  if ("time" %in% names(nc$dim)) {
    time_vals <- nc$dim$time$vals
    time_units <- nc$dim$time$units
    
    if (!is.null(time_vals) && length(time_vals) > 0) {
      return(list(vals = time_vals, units = time_units))
    }
  }
  
  if ("time" %in% names(nc$var)) {
    time_vals <- ncdf4::ncvar_get(nc, "time")
    time_units <- ncdf4::ncatt_get(nc, "time", "units")$value
    
    if (!is.null(time_vals) && length(time_vals) > 0) {
      return(list(vals = time_vals, units = time_units))
    }
  }
  
  stop("Cannot find valid time dimension or time variable in NetCDF.", call. = FALSE)
}


parse_nc_time_origin <- function(time_units) {
  
  if (is.null(time_units) || !grepl("^days since ", time_units)) {
    stop("Unsupported NetCDF time units: ", time_units, call. = FALSE)
  }
  
  origin_str <- sub("^days since ", "", time_units)
  
  origin_time <- lubridate::ymd_hms(origin_str, tz = "UTC", quiet = TRUE)
  
  if (is.na(origin_time)) {
    origin_time <- lubridate::ymd(origin_str, tz = "UTC", quiet = TRUE)
  }
  
  if (is.na(origin_time)) {
    stop("Cannot parse NetCDF time origin: ", origin_str, call. = FALSE)
  }
  
  as.POSIXct(origin_time, tz = "UTC")
}


get_last_time_index_matching_stop <- function(nc, stop.time, tolerance_seconds = 3600) {
  
  stop.time <- parse_sipnet_datetime(stop.time)
  
  time_info <- get_nc_time_values(nc)
  time_vals <- as.numeric(time_info$vals)
  time_units <- time_info$units
  
  origin_time <- parse_nc_time_origin(time_units)
  
  # IMPORTANT:
  # NetCDF time values are often fractional days, e.g. 0, 0.125, 0.25.
  # Do NOT use lubridate::days(time_vals), because Period requires integer values.
  # Convert fractional days to seconds instead.
  real_time <- origin_time + time_vals * 86400
  
  # Prefer exact / near stop.time match.
  diff_seconds <- abs(as.numeric(difftime(real_time, stop.time, units = "secs")))
  near_idx <- which(diff_seconds <= tolerance_seconds)
  
  if (length(near_idx) > 0) {
    return(near_idx[length(near_idx)])
  }
  
  # Fallback: last timestep not after stop.time.
  idx <- which(real_time <= stop.time)
  
  if (length(idx) > 0) {
    return(idx[length(idx)])
  }
  
  stop(
    "No NetCDF timestep matches or precedes stop.time.",
    "\nstop.time = ", as.character(stop.time),
    "\nNetCDF time range = ",
    as.character(min(real_time, na.rm = TRUE)),
    " to ",
    as.character(max(real_time, na.rm = TRUE)),
    call. = FALSE
  )
}



get_nc_var_or_null <- function(nc, v) {
  if (!(v %in% names(nc$var))) {
    return(NULL)
  }
  ncdf4::ncvar_get(nc, v)
}


get_state_at_time <- function(ens, v, last) {
  
  x <- ens[[v]]
  
  if (is.null(x)) {
    stop("Variable `", v, "` is missing from NetCDF.", call. = FALSE)
  }
  
  x <- as.numeric(x)
  
  if (length(x) < last || is.na(last)) {
    stop(
      "Invalid time index for variable `", v, "`.",
      " last = ", last,
      "; length = ", length(x),
      call. = FALSE
    )
  }
  
  x[last]
}


get_interval_mean <- function(ens, v) {
  
  x <- ens[[v]]
  
  if (is.null(x)) {
    stop("Variable `", v, "` is missing from NetCDF.", call. = FALSE)
  }
  
  x <- as.numeric(x)
  mean(x, na.rm = TRUE)
}


nee_model_to_obs <- function(x) {
  x * 1e8 / 1.157407
}


nee_obs_to_model <- function(x) {
  x * 1.157407 / 1e8
}

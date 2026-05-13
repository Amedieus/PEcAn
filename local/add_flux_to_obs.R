library(data.table)
library(lubridate)

# ============================================================
# 0. User settings
# ============================================================

obs_rdata_dir <- "/projectnb/dietzelab/guYANG/pecan/output/obs/Rdata"

nee_file <- "/projectnb/dietzelab/guYANG/Gap_fill/results/ens_ec_3h.csv"
qle_file <- "/projectnb/dietzelab/guYANG/Gap_fill/results/le_ens_ec_3h.csv"
match_file <- "/projectnb/dietzelab/guYANG/Validation/validation/matched_within_1km.csv"

out_obs_mean_file <- file.path(obs_rdata_dir, "obs.mean.monthly.6obs.Rdata")
out_obs_cov_file  <- file.path(obs_rdata_dir, "obs.cov.monthly.6obs.Rdata")

# Your SDA / XML state order.
# If SMAP has been removed from obs, SoilMoistFrac can still be omitted here.
state_vars_xml <- c(
  "AbvGrndWood",
  "NEE",
  "Qle",
  "LAI",
  "TotSoilCarb"
)

# Inflation for NEE/Qle observation variances.
flux_var_inflation <- 1.1

# ============================================================
# 2. Helper: aggregate 3-hourly ensemble EC data to annual SDA time
# ============================================================

aggregate_flux_to_annual <- function(flux_file,
                                     match_file,
                                     flux_name = c("NEE", "Qle")) {
  flux_name <- match.arg(flux_name)
  
  flux_df <- fread(flux_file)
  close_points_df <- fread(match_file)
  
  setDT(flux_df)
  setDT(close_points_df)
  
  # Map Site_ID to the index used by your obs.mean list names.
  idx_map <- close_points_df[, .(Site_ID, index)]
  setkey(idx_map, Site_ID)
  
  flux_df <- idx_map[flux_df, on = "Site_ID"]
  dt <- as.data.table(flux_df)
  
  # Clean time
  dt <- dt[!is.na(utc)]
  dt[, utc := as.POSIXct(utc, tz = "UTC")]
  
  # Define assimilation year anchored on July 15.
  # <= July 15 goes to current year July 15.
  # > July 15 goes to next year July 15.
  dt[, assim_year := ifelse(
    month(utc) > 7 | (month(utc) == 7 & day(utc) > 15),
    year(utc) + 1L,
    year(utc)
  )]
  
  dt[, time := as.Date(paste0(assim_year, "-07-15"))]
  
  # Ensemble columns: ens01, ens02, ...
  ens_cols <- grep("^ens\\d{2}$", names(dt), value = TRUE)
  ens_cols <- sort(ens_cols)
  stopifnot(length(ens_cols) >= 1)
  stopifnot("ens_mean" %in% names(dt))
  
  # Annual mean for each index x time.
  annual_df <- dt[
    ,
    c(
      lapply(.SD, mean, na.rm = TRUE),
      list(ens_mean = mean(ens_mean, na.rm = TRUE))
    ),
    by = .(index, time),
    .SDcols = ens_cols
  ]
  
  # Remove rows where all ensemble values and ens_mean are NA.
  cols_check <- c(ens_cols, "ens_mean")
  annual_df <- annual_df[
    rowSums(!is.na(annual_df[, ..cols_check])) > 0
  ]
  
  # Remove 2025 if you do not want it.
  annual_df <- annual_df[time != as.Date("2025-07-15")]
  
  # Rename ensemble columns.
  old_names <- c(ens_cols, "ens_mean")
  new_names <- c(paste0(flux_name, "_", ens_cols), flux_name)
  setnames(annual_df, old = old_names, new = new_names)
  
  flux_ens_cols <- paste0(flux_name, "_", ens_cols)
  
  # Ensemble spread based on annual ensemble means.
  annual_df[
    ,
    paste0(flux_name, "_sd") := apply(.SD, 1, sd, na.rm = TRUE),
    .SDcols = flux_ens_cols
  ]
  
  annual_df[
    ,
    paste0(flux_name, "_var") := apply(.SD, 1, var, na.rm = TRUE),
    .SDcols = flux_ens_cols
  ]
  
  # Keep only what we need.
  keep_cols <- c("index", "time", flux_name,
                 paste0(flux_name, "_sd"),
                 paste0(flux_name, "_var"))
  
  annual_df <- annual_df[, ..keep_cols]
  annual_df[, index := as.character(index)]
  annual_df[, time := as.Date(time)]
  
  setorder(annual_df, index, time)
  
  annual_df
}


# ============================================================
# 3. Build annual flux observation table: NEE + Qle
# ============================================================

nee_annual_df <- aggregate_flux_to_annual(
  flux_file = nee_file,
  match_file = match_file,
  flux_name = "NEE"
)

qle_annual_df <- aggregate_flux_to_annual(
  flux_file = qle_file,
  match_file = match_file,
  flux_name = "Qle"
)

# Merge NEE and Qle.
# all = TRUE means if one exists but the other is missing, still keep the available one.
annual_flux_df <- merge(
  nee_annual_df,
  qle_annual_df,
  by = c("index", "time"),
  all = TRUE
)

annual_flux_df[, index := as.character(index)]
annual_flux_df[, time := as.Date(time)]

setorder(annual_flux_df, index, time)

message("Annual flux table:")
print(dim(annual_flux_df))
print(head(annual_flux_df))


# ============================================================
# 4. Diagnostics before insertion
# ============================================================

obs_dates <- names(obs.mean)
flux_dates <- as.character(unique(annual_flux_df$time))

message("Date overlap between obs.mean and flux table:")
print(length(intersect(obs_dates, flux_dates)))
print(intersect(obs_dates, flux_dates))

first_date <- names(obs.mean)[1]
obs_site_names <- names(obs.mean[[first_date]])

message("Example obs.mean site names:")
print(head(obs_site_names))

message("Example annual_flux_df index values:")
print(head(unique(annual_flux_df$index)))

message("Site-name overlap:")
print(length(intersect(obs_site_names, unique(annual_flux_df$index))))


# ============================================================
# 5. Add NEE and Qle to obs.mean by site name
# ============================================================

add_flux_to_obsmean_by_name <- function(obs.mean,
                                        annual_flux_df,
                                        state_vars_xml) {
  obs.mean_new <- obs.mean
  
  flux <- as.data.table(copy(annual_flux_df))
  flux[, time_chr := as.character(as.Date(time))]
  flux[, site_chr := as.character(index)]
  
  n_add_nee <- 0L
  n_add_qle <- 0L
  n_skip_date <- 0L
  n_skip_site <- 0L
  n_skip_no_flux <- 0L
  
  for (r in seq_len(nrow(flux))) {
    d <- flux$time_chr[r]
    sid <- flux$site_chr[r]
    
    if (!(d %in% names(obs.mean_new))) {
      n_skip_date <- n_skip_date + 1L
      next
    }
    
    if (!(sid %in% names(obs.mean_new[[d]]))) {
      n_skip_site <- n_skip_site + 1L
      next
    }
    
    nee <- flux$NEE[r]
    qle <- flux$Qle[r]
    
    if (is.na(nee) && is.na(qle)) {
      n_skip_no_flux <- n_skip_no_flux + 1L
      next
    }
    
    df0 <- obs.mean_new[[d]][[sid]]
    
    if (is.null(df0)) {
      df0 <- data.frame()
    }
    
    if (!is.data.frame(df0)) {
      df0 <- as.data.frame(df0)
    }
    
    if (!is.na(nee)) {
      df0$NEE <- nee
      n_add_nee <- n_add_nee + 1L
    }
    
    if (!is.na(qle)) {
      df0$Qle <- qle
      n_add_qle <- n_add_qle + 1L
    }
    
    # Reorder variables according to SDA/XML order.
    vars_present <- names(df0)
    vars_ordered <- state_vars_xml[state_vars_xml %in% vars_present]
    extra_vars <- setdiff(vars_present, vars_ordered)
    
    df0 <- df0[, c(vars_ordered, extra_vars), drop = FALSE]
    
    obs.mean_new[[d]][[sid]] <- df0
  }
  
  message("Added NEE values: ", n_add_nee)
  message("Added Qle values: ", n_add_qle)
  message("Skipped rows because date not found: ", n_skip_date)
  message("Skipped rows because site not found: ", n_skip_site)
  message("Skipped rows because both NEE and Qle are NA: ", n_skip_no_flux)
  
  obs.mean_new
}

obs.mean_new <- add_flux_to_obsmean_by_name(
  obs.mean = obs.mean,
  annual_flux_df = annual_flux_df,
  state_vars_xml = state_vars_xml
)


# ============================================================
# 6. Add NEE and Qle variances to obs.cov by site name
# ============================================================

add_flux_to_obscov_by_name <- function(obs.mean_old,
                                       obs.cov_old,
                                       obs.mean_new,
                                       annual_flux_df,
                                       state_vars_xml,
                                       inflation = 1.1) {
  obs.cov_new <- obs.cov_old
  
  flux <- as.data.table(copy(annual_flux_df))
  flux[, time_chr := as.character(as.Date(time))]
  flux[, site_chr := as.character(index)]
  setkey(flux, time_chr, site_chr)
  
  n_add_nee_var <- 0L
  n_add_qle_var <- 0L
  n_drop_missing_var <- 0L
  
  for (d in names(obs.mean_new)) {
    for (sid in names(obs.mean_new[[d]])) {
      
      y_new <- obs.mean_new[[d]][[sid]]
      
      if (is.null(y_new) || !is.data.frame(y_new) || ncol(y_new) == 0) {
        obs.cov_new[[d]][[sid]] <- matrix(0, 0, 0)
        next
      }
      
      vars_new <- names(y_new)
      
      # Old observation mean/cov for old variables.
      y_old <- NULL
      R_old <- NULL
      
      if (d %in% names(obs.mean_old) &&
          sid %in% names(obs.mean_old[[d]])) {
        y_old <- obs.mean_old[[d]][[sid]]
      }
      
      if (d %in% names(obs.cov_old) &&
          sid %in% names(obs.cov_old[[d]])) {
        R_old <- obs.cov_old[[d]][[sid]]
      }
      
      old_var_diag <- numeric(0)
      
      if (!is.null(y_old) &&
          is.data.frame(y_old) &&
          ncol(y_old) > 0 &&
          !is.null(R_old)) {
        
        R_old <- as.matrix(R_old)
        old_vars <- names(y_old)
        
        if (nrow(R_old) == length(old_vars) &&
            ncol(R_old) == length(old_vars)) {
          old_var_diag <- diag(R_old)
          names(old_var_diag) <- old_vars
        }
      }
      
      # New variance vector in the same order as obs.mean_new.
      vv <- rep(NA_real_, length(vars_new))
      names(vv) <- vars_new
      
      # Fill old variables from old obs.cov.
      common_old <- intersect(names(old_var_diag), vars_new)
      vv[common_old] <- old_var_diag[common_old]
      
      # Fill NEE/Qle from annual flux ensemble spread.
      row_flux <- flux[J(d, sid), nomatch = 0]
      
      if (nrow(row_flux) > 0) {
        row_flux <- row_flux[1]
        
        if ("NEE" %in% vars_new && !is.na(row_flux$NEE_var)) {
          vv["NEE"] <- row_flux$NEE_var * inflation
          n_add_nee_var <- n_add_nee_var + 1L
        }
        
        if ("Qle" %in% vars_new && !is.na(row_flux$Qle_var)) {
          vv["Qle"] <- row_flux$Qle_var * inflation
          n_add_qle_var <- n_add_qle_var + 1L
        }
      }
      
      # If a variable has mean but no variance, drop it to avoid dimension mismatch.
      ok <- !is.na(vv)
      
      if (!all(ok)) {
        vars_keep <- names(vv)[ok]
        vars_drop <- names(vv)[!ok]
        
        if (length(vars_drop) > 0) {
          n_drop_missing_var <- n_drop_missing_var + length(vars_drop)
        }
        
        y_new <- y_new[, vars_keep, drop = FALSE]
        obs.mean_new[[d]][[sid]] <- y_new
        
        vv <- vv[ok]
        vars_new <- names(vv)
      }
      
      if (length(vv) == 0) {
        obs.cov_new[[d]][[sid]] <- matrix(0, 0, 0)
      } else {
        R_new <- diag(as.numeric(vv), nrow = length(vv), ncol = length(vv))
        dimnames(R_new) <- list(vars_new, vars_new)
        obs.cov_new[[d]][[sid]] <- R_new
      }
    }
  }
  
  message("Added NEE variances: ", n_add_nee_var)
  message("Added Qle variances: ", n_add_qle_var)
  message("Dropped variables because variance was missing: ", n_drop_missing_var)
  
  list(
    obs.mean = obs.mean_new,
    obs.cov = obs.cov_new
  )
}

tmp <- add_flux_to_obscov_by_name(
  obs.mean_old = obs.mean,
  obs.cov_old = obs.cov,
  obs.mean_new = obs.mean_new,
  annual_flux_df = annual_flux_df,
  state_vars_xml = state_vars_xml,
  inflation = flux_var_inflation
)

obs.mean_new <- tmp$obs.mean
obs.cov_new  <- tmp$obs.cov


# ============================================================
# 7. Final format checks
# ============================================================

check_obs_pair <- function(obs.mean, obs.cov) {
  bad <- list()
  k <- 1L
  
  for (d in names(obs.mean)) {
    if (!(d %in% names(obs.cov))) {
      bad[[k]] <- data.frame(time = d, site = NA, issue = "date missing in obs.cov")
      k <- k + 1L
      next
    }
    
    for (sid in names(obs.mean[[d]])) {
      y <- obs.mean[[d]][[sid]]
      
      if (!(sid %in% names(obs.cov[[d]]))) {
        bad[[k]] <- data.frame(time = d, site = sid, issue = "site missing in obs.cov")
        k <- k + 1L
        next
      }
      
      R <- obs.cov[[d]][[sid]]
      
      if (is.null(y) || !is.data.frame(y) || ncol(y) == 0) next
      
      if (is.null(R)) {
        bad[[k]] <- data.frame(time = d, site = sid, issue = "R is NULL")
        k <- k + 1L
        next
      }
      
      R <- as.matrix(R)
      
      if (ncol(y) != nrow(R) || ncol(y) != ncol(R)) {
        bad[[k]] <- data.frame(
          time = d,
          site = sid,
          issue = paste0(
            "dim mismatch: y=",
            ncol(y),
            ", R=",
            paste(dim(R), collapse = "x")
          )
        )
        k <- k + 1L
        next
      }
      
      if (!identical(names(y), colnames(R))) {
        bad[[k]] <- data.frame(
          time = d,
          site = sid,
          issue = paste0(
            "name mismatch: y=",
            paste(names(y), collapse = ","),
            " | R=",
            paste(colnames(R), collapse = ",")
          )
        )
        k <- k + 1L
      }
    }
  }
  
  if (length(bad) == 0) {
    message("obs.mean and obs.cov are aligned.")
    return(invisible(NULL))
  }
  
  rbindlist(bad, fill = TRUE)
}

check_result <- check_obs_pair(obs.mean_new, obs.cov_new)

if (!is.null(check_result)) {
  print(check_result)
  stop("obs.mean and obs.cov are not aligned. Fix before saving.")
}


# ============================================================
# 8. Confirm NEE and Qle are present
# ============================================================

vars_all <- unique(unlist(lapply(obs.mean_new, function(x) {
  unlist(lapply(x, names))
})))

message("Variables found in obs.mean_new:")
print(vars_all)

message("Number of site-time entries containing NEE:")
print(sum(unlist(lapply(obs.mean_new, function(date_list) {
  lapply(date_list, function(df) {
    is.data.frame(df) && "NEE" %in% names(df)
  })
}))))

message("Number of site-time entries containing Qle:")
print(sum(unlist(lapply(obs.mean_new, function(date_list) {
  lapply(date_list, function(df) {
    is.data.frame(df) && "Qle" %in% names(df)
  })
}))))

# ============================================================
# 9. Save using object names expected by SDA
# ============================================================

obs.mean <- obs.mean_new
obs.cov  <- obs.cov_new

save(obs.mean, file = out_obs_mean_file)
save(obs.cov,  file = out_obs_cov_file)

message("Saved obs.mean to: ", out_obs_mean_file)
message("Saved obs.cov to: ", out_obs_cov_file)